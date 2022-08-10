/*
 * IndexerPlain.cpp
 *
 * SimpleMonochromaticDiffractionPatternPrediction.h
 *
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019      Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
 *
 * This file is part of XGANDALF.
 *
 * XGANDALF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * XGANDALF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with XGANDALF.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <xgandalf/IndexerPlain.h>
#include <algorithm>
#include <sstream>
#include <vector>

using namespace Eigen;
using namespace std;

namespace xgandalf
{
    IndexerPlain::IndexerPlain(const ExperimentSettings& experimentSettings)
        : IndexerBase(experimentSettings)
    {
        precompute();
    }

    void IndexerPlain::precompute()
    {
        maxCloseToPointDeviation = 0.15;
        maxPeaksToUseForIndexing = 250;

        setSamplingPitch(SamplingPitch::standard);
        setGradientDescentIterationsCount(GradientDescentIterationsCount::standard);

        float minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.2;
        float maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.2;

        try
        {
            sparsePeakFinder.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);
        }
        catch (const BadInputException& e)
        {
            cerr << "\nThe difference between largest and smallest possible lattice vector is very big. Reducing some parameters to nevertheless make fitting "
                    "possible.\n This will reduce the fitting performance. Consider reducing the quotient between largest and smallest possible lattice "
                    "vector.\n";
            minSpacingBetweenPeaks = experimentSettings.getDifferentRealLatticeVectorLengths_A().minCoeff() * 0.4;
            maxPossiblePointNorm = experimentSettings.getDifferentRealLatticeVectorLengths_A().maxCoeff() * 1.1;

            try
            {
                sparsePeakFinder.precompute(minSpacingBetweenPeaks, maxPossiblePointNorm);
            }
            catch (const BadInputException& e)
            {
                stringstream errStream;
                errStream
                    << "\nThe quotient between largest and smallest possible lattice vector is too big. To be able to run, make it smaller or contact the "
                       "developers of xgandalf.\n";
                throw BadInputException(errStream.str());
            }
        }

        inverseSpaceTransform = InverseSpaceTransform(maxCloseToPointDeviation);
        inverseSpaceTransform.setFunctionSelection(9);
        inverseSpaceTransform.setOptionalFunctionArgument(8);
        inverseSpaceTransform.setLocalTransformFlag();
        inverseSpaceTransform.clearRadialWeightingFlag();

        accuracyConstants_LatticeAssembler.maxCountGlobalPassingWeightFilter = 500;
        accuracyConstants_LatticeAssembler.maxCountLocalPassingWeightFilter = 15;
        accuracyConstants_LatticeAssembler.maxCountPassingRelativeDefectFilter = 50;
        accuracyConstants_LatticeAssembler.minPointsOnLattice = 5;
        accuracyConstants_LatticeAssembler.maxCloseToPointDeviation = maxCloseToPointDeviation;
        accuracyConstants_LatticeAssembler.refineWithExactLattice = false;
        latticeAssembler.setAccuracyConstants(accuracyConstants_LatticeAssembler);


        if (experimentSettings.isLatticeParametersKnown())
        {
            //    latticeAssembler.setDeterminantRange(experimentSettings.getMinRealLatticeDeterminant_A3(),
            //    experimentSettings.getMaxRealLatticeDeterminant_A3());
            latticeAssembler.setDeterminantRange(experimentSettings.getRealLatticeDeterminant_A3() * 0.8,
                                                 experimentSettings.getRealLatticeDeterminant_A3() * 1.2); // debug
            latticeAssembler.setKnownLatticeParameters(experimentSettings.getSampleRealLattice_A(), experimentSettings.getTolerance());
        }
        else
        {
            latticeAssembler.setDeterminantRange(experimentSettings.getMinRealLatticeDeterminant_A3(), experimentSettings.getMaxRealLatticeDeterminant_A3());
        }
    }

    void IndexerPlain::setSamplingPitch(SamplingPitch samplingPitch)
    {
        float unitPitch = 0;
        bool coverSecondaryMillerIndices = false;

        switch (samplingPitch)
        {
            case SamplingPitch::extremelyLoose:
                unitPitch = 0.10;
                break;
            case SamplingPitch::loose:
                unitPitch = 0.075;
                break;
            case SamplingPitch::standard:
                unitPitch = 0.05;
                break;
            case SamplingPitch::dense:
                unitPitch = 0.025;
                break;
            case SamplingPitch::extremelyDense:
                unitPitch = 0.01;
                break;

            case SamplingPitch::standardWithSeondaryMillerIndices:
                unitPitch = 0.05;
                coverSecondaryMillerIndices = true;
                break;
            case SamplingPitch::denseWithSeondaryMillerIndices:
                unitPitch = 0.025;
                coverSecondaryMillerIndices = true;
                break;
            case SamplingPitch::extremelyDenseWithSeondaryMillerIndices:
                unitPitch = 0.01;
                coverSecondaryMillerIndices = true;
                break;
        }

        setSamplingPitch(unitPitch, coverSecondaryMillerIndices);
    }

    void IndexerPlain::setSamplingPitch(float unitPitch, bool coverSecondaryMillerIndices)
    {
        if (experimentSettings.isLatticeParametersKnown())
        {
            // float tolerance = max(unitPitch, experimentSettings.getTolerance());
            float tolerance = experimentSettings.getTolerance();

            if (!coverSecondaryMillerIndices)
            {
                samplePointsGenerator.getTightGrid(precomputedSamplePoints, unitPitch, tolerance, experimentSettings.getDifferentRealLatticeVectorLengths_A());
            }
            else
            {
                const Matrix3f& sampleRealBasis = experimentSettings.getSampleRealLattice_A().getBasis();

                vector<float> radii(6);
                radii[0] = sampleRealBasis.col(0).norm();
                radii[1] = sampleRealBasis.col(1).norm();
                radii[2] = sampleRealBasis.col(2).norm();
                radii[3] = (sampleRealBasis.col(0) + sampleRealBasis.col(1)).norm();
                radii[4] = (sampleRealBasis.col(1) + sampleRealBasis.col(2)).norm();
                radii[5] = (sampleRealBasis.col(2) + sampleRealBasis.col(0)).norm();

                sort(radii.begin(), radii.end());
                auto it = std::unique(radii.begin(), radii.end());
                radii.resize(std::distance(radii.begin(), it));

                ArrayXf radii_array = Eigen::Map<ArrayXf>(radii.data(), radii.size(), 1);

                samplePointsGenerator.getTightGrid(precomputedSamplePoints, unitPitch, tolerance, radii_array);
            }
        }
        else
        {
            if (!coverSecondaryMillerIndices)
            {
                float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
                float maxRadius = experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

                samplePointsGenerator.getDenseGrid(precomputedSamplePoints, unitPitch, minRadius, maxRadius);
            }
            else
            {
                float minRadius = experimentSettings.getMinRealLatticeVectorLength_A() * 0.98;
                float maxRadius = 2 * experimentSettings.getMaxRealLatticeVectorLength_A() * 1.02;

                samplePointsGenerator.getDenseGrid(precomputedSamplePoints, unitPitch, minRadius, maxRadius);
            }
        }
    }

    void IndexerPlain::setRefineWithExactLattice(bool flag)
    {
        LatticeAssembler::accuracyConstants_t accuracyConstants = latticeAssembler.getAccuracyConstants();

        accuracyConstants.refineWithExactLattice = flag;

        latticeAssembler.setAccuracyConstants(accuracyConstants);
    }

    void IndexerPlain::setMaxPeaksToUseForIndexing(int maxPeaksToUseForIndexing)
    {
        this->maxPeaksToUseForIndexing = maxPeaksToUseForIndexing;
    }

    void IndexerPlain::index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A)
    {
        vector<int> peakCountOnLattices;
        index(assembledLattices, reciprocalPeaks_1_per_A, peakCountOnLattices);
    }

    void IndexerPlain::index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A, std::vector<int>& peakCountOnLattices)
    {
        if (precomputedSamplePoints.size() == 0)
        {
            precompute();
        }

        Matrix3Xf samplePoints = precomputedSamplePoints;

        Matrix3Xf reciprocalPeaksReduced_1_per_A = reciprocalPeaks_1_per_A;
        reducePeakCount(reciprocalPeaksReduced_1_per_A);

        // global hill climbing
        hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_global);
        hillClimbingOptimizer.performOptimization(reciprocalPeaksReduced_1_per_A, samplePoints);
        RowVectorXf globalHillClimbingPointEvaluation = hillClimbingOptimizer.getLastInverseTransformEvaluation();
        Matrix3Xf globalHillClimbingSamplePoints = samplePoints;

        // additional global hill climbing
        hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_additionalGlobal);
        hillClimbingOptimizer.performOptimization(reciprocalPeaksReduced_1_per_A, samplePoints);
        RowVectorXf& additionalGlobalHillClimbingPointEvaluation = hillClimbingOptimizer.getLastInverseTransformEvaluation();
        Matrix3Xf& additionalGlobalHillClimbingSamplePoints = samplePoints;

        // find peaks
        uint32_t maxGlobalPeaksToTakeCount = 50;
        sparsePeakFinder.findPeaks_fast(globalHillClimbingSamplePoints, globalHillClimbingPointEvaluation);
        keepSamplePointsWithHighestEvaluation(globalHillClimbingSamplePoints, globalHillClimbingPointEvaluation, maxGlobalPeaksToTakeCount);

        uint32_t maxAdditionalGlobalPeaksToTakeCount = 50;
        sparsePeakFinder.findPeaks_fast(additionalGlobalHillClimbingSamplePoints, additionalGlobalHillClimbingPointEvaluation);
        keepSamplePointsWithHighestEvaluation(additionalGlobalHillClimbingSamplePoints, additionalGlobalHillClimbingPointEvaluation,
                                              maxAdditionalGlobalPeaksToTakeCount);

        Matrix3Xf peakSamplePoints(3, globalHillClimbingSamplePoints.cols() + additionalGlobalHillClimbingSamplePoints.cols());
        peakSamplePoints << globalHillClimbingSamplePoints, additionalGlobalHillClimbingSamplePoints;

        // peaks hill climbing
        hillClimbingOptimizer.setHillClimbingAccuracyConstants(hillClimbing_accuracyConstants_peaks);
        hillClimbingOptimizer.performOptimization(reciprocalPeaksReduced_1_per_A, peakSamplePoints);

        // final peaks extra evaluation
        inverseSpaceTransform.setPointsToTransform(reciprocalPeaks_1_per_A);
        inverseSpaceTransform.performTransform(peakSamplePoints);

        // find peaks , TODO: check, whether better performance without peak finding here
        // sparsePeakFinder.findPeaks_fast(peakSamplePoints, inverseSpaceTransform.getInverseTransformEvaluation());

        // assemble lattices
        vector<LatticeAssembler::assembledLatticeStatistics_t> assembledLatticesStatistics;
        Matrix3Xf& candidateVectors = peakSamplePoints;
        RowVectorXf& candidateVectorWeights = inverseSpaceTransform.getInverseTransformEvaluation();
        vector<vector<uint16_t>>& pointIndicesOnVector = inverseSpaceTransform.getPointsCloseToEvaluationPositions_indices();
        Matrix3Xf reciprocalPeaksCopy_1_per_A = reciprocalPeaks_1_per_A;
        latticeAssembler.assembleLattices(assembledLattices, assembledLatticesStatistics, candidateVectors, candidateVectorWeights, pointIndicesOnVector,
                                          reciprocalPeaksCopy_1_per_A);

        peakCountOnLattices.clear();
        peakCountOnLattices.reserve(assembledLatticesStatistics.size());
        for (auto assembledLatticeStatistics = assembledLatticesStatistics.cbegin(); assembledLatticeStatistics != assembledLatticesStatistics.cend();
             ++assembledLatticeStatistics)
        {
            peakCountOnLattices.push_back(assembledLatticeStatistics->occupiedLatticePointsCount);
        }


        //    cout << assembledLatticesStatistics[0].meanDefect << " " << assembledLatticesStatistics[0].meanRelativeDefect << " "
        //            << assembledLatticesStatistics[0].occupiedLatticePointsCount << " " << assembledLatticesStatistics.size() << endl <<
        //            assembledLattices[0].det()
        //            << endl;

        //    ofstream ofs("workfolder/samplePoints", ofstream::out);
        //    ofs << samplePoints.transpose().eval();
    }


    void IndexerPlain::setGradientDescentIterationsCount(GradientDescentIterationsCount gradientDescentIterationsCount)
    {
        HillClimbingOptimizer::hillClimbingAccuracyConstants_t& global = hillClimbing_accuracyConstants_global;
        HillClimbingOptimizer::hillClimbingAccuracyConstants_t& additionalGlobal = hillClimbing_accuracyConstants_additionalGlobal;
        HillClimbingOptimizer::hillClimbingAccuracyConstants_t& peaks = hillClimbing_accuracyConstants_peaks;

        float meanRealLatticeVectorLength = experimentSettings.getDifferentRealLatticeVectorLengths_A().mean();

        switch (gradientDescentIterationsCount)
        {
            case GradientDescentIterationsCount::exremelyFew:
                global.initialIterationCount = 10;
                global.calmDownIterationCount = 3;
                global.calmDownFactor = 0.65;
                global.localFitIterationCount = 4;
                global.localCalmDownIterationCount = 3;
                global.localCalmDownFactor = 0.65;

                global.stepComputationAccuracyConstants.gamma = 0.75;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 15;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 150;
                global.stepComputationAccuracyConstants.directionChangeFactor = 2;
                break;
            case GradientDescentIterationsCount::few:
                global.initialIterationCount = 20;
                global.calmDownIterationCount = 4;
                global.calmDownFactor = 0.73;
                global.localFitIterationCount = 6;
                global.localCalmDownIterationCount = 5;
                global.localCalmDownFactor = 0.73;

                global.stepComputationAccuracyConstants.gamma = 0.7;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 18;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 180;
                global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

                break;
            case GradientDescentIterationsCount::standard:
                global.initialIterationCount = 40;
                global.calmDownIterationCount = 5;
                global.calmDownFactor = 0.8;
                global.localFitIterationCount = 8;
                global.localCalmDownIterationCount = 6;
                global.localCalmDownFactor = 0.8;

                global.stepComputationAccuracyConstants.gamma = 0.65;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 20;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 200;
                global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

                break;
            case GradientDescentIterationsCount::many:
                global.initialIterationCount = 60;
                global.calmDownIterationCount = 8;
                global.calmDownFactor = 0.83;
                global.localFitIterationCount = 10;
                global.localCalmDownIterationCount = 9;
                global.localCalmDownFactor = 0.83;

                global.stepComputationAccuracyConstants.gamma = 0.65;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 20;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 200;
                global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

                break;
            case GradientDescentIterationsCount::manyMany:
                global.initialIterationCount = 90;
                global.calmDownIterationCount = 10;
                global.calmDownFactor = 0.85;
                global.localFitIterationCount = 15;
                global.localCalmDownIterationCount = 15;
                global.localCalmDownFactor = 0.9;

                global.stepComputationAccuracyConstants.gamma = 0.60;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 23;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 230;
                global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

                break;
            case GradientDescentIterationsCount::extremelyMany:
                global.initialIterationCount = 200;
                global.calmDownIterationCount = 15;
                global.calmDownFactor = 0.9;
                global.localFitIterationCount = 15;
                global.localCalmDownIterationCount = 15;
                global.localCalmDownFactor = 0.9;

                global.stepComputationAccuracyConstants.gamma = 0.60;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 23;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 230;
                global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;

                break;
            case GradientDescentIterationsCount::custom:
                global.initialIterationCount = 40;
                global.calmDownIterationCount = 5;
                global.calmDownFactor = 0.8;
                global.localFitIterationCount = 8;
                global.localCalmDownIterationCount = 6;
                global.localCalmDownFactor = 0.8;

                global.stepComputationAccuracyConstants.gamma = 0.65;
                global.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 20;
                global.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 200;
                global.stepComputationAccuracyConstants.directionChangeFactor = 1.5;
                break;
        }

        global.functionSelection = 1;
        global.optionalFunctionArgument = 1;
        global.maxCloseToPointDeviation = maxCloseToPointDeviation;

        additionalGlobal.functionSelection = 9;
        additionalGlobal.optionalFunctionArgument = 4;
        additionalGlobal.maxCloseToPointDeviation = maxCloseToPointDeviation;

        additionalGlobal.initialIterationCount = 0;
        additionalGlobal.calmDownIterationCount = 0;
        additionalGlobal.calmDownFactor = 0;
        additionalGlobal.localFitIterationCount = 4;
        additionalGlobal.localCalmDownIterationCount = 3;
        additionalGlobal.localCalmDownFactor = 0.7;

        additionalGlobal.stepComputationAccuracyConstants.gamma = 0.65;
        additionalGlobal.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 50;
        additionalGlobal.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 20000;
        additionalGlobal.stepComputationAccuracyConstants.directionChangeFactor = 2.5;

        peaks.functionSelection = 9;
        peaks.optionalFunctionArgument = 4;
        peaks.maxCloseToPointDeviation = maxCloseToPointDeviation;

        peaks.initialIterationCount = 0;
        peaks.calmDownIterationCount = 0;
        peaks.calmDownFactor = 0;
        peaks.localFitIterationCount = 20;
        peaks.localCalmDownIterationCount = 300;
        peaks.localCalmDownFactor = 0.98;

        peaks.stepComputationAccuracyConstants.gamma = 0.1;
        peaks.stepComputationAccuracyConstants.maxStep = meanRealLatticeVectorLength / 300;
        peaks.stepComputationAccuracyConstants.minStep = meanRealLatticeVectorLength / 20000;
        peaks.stepComputationAccuracyConstants.directionChangeFactor = 2.5;
    }

    void IndexerPlain::reducePeakCount(Matrix3Xf& reciprocalPeaks_1_per_A)
    {
        if (maxPeaksToUseForIndexing >= reciprocalPeaks_1_per_A.cols())
            return;

        RowVectorXf norms = reciprocalPeaks_1_per_A.colwise().norm();
        RowVectorXf norms_sorted = norms;
        sort((float*)norms_sorted.data(), (float*)norms_sorted.data() + norms_sorted.size());
        float maxNorm = norms_sorted[maxPeaksToUseForIndexing - 1];

        int peaksKept_norms = 0;
        for (int i = 0; i < norms.cols(); i++)
        {
            if (norms[i] <= maxNorm)
            {
                reciprocalPeaks_1_per_A.col(peaksKept_norms) = reciprocalPeaks_1_per_A.col(i);
                peaksKept_norms++;
            }
        }

        reciprocalPeaks_1_per_A.conservativeResize(NoChange, maxPeaksToUseForIndexing);
    }
} // namespace xgandalf
