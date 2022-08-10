/*
 * LatticeAssembler.cpp
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

#include <xgandalf/refinement.h>
#include <xgandalf/LatticeAssembler.h>

#include <algorithm>
#include <ctype.h>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <iostream>

namespace xgandalf
{
    using namespace Eigen;
    using namespace std;

    LatticeAssembler::LatticeAssembler()
    {
        setStandardValues();
    }
    LatticeAssembler::LatticeAssembler(const Vector2f& determinantRange)
    {
        setStandardValues();
        setDeterminantRange(determinantRange);
    }
    LatticeAssembler::LatticeAssembler(const Vector2f& determinantRange, const Lattice& sampleRealLattice_A, float knownLatticeTolerance)
    {
        setStandardValues();
        setDeterminantRange(determinantRange);
        setKnownLatticeParameters(sampleRealLattice_A, knownLatticeTolerance);
    }
    LatticeAssembler::LatticeAssembler(const Vector2f& determinantRange, const accuracyConstants_t& accuracyConstants)
    {
        setStandardValues();
        setDeterminantRange(determinantRange);
        setAccuracyConstants(accuracyConstants);
    }
    LatticeAssembler::LatticeAssembler(const Vector2f& determinantRange, const Lattice& sampleRealLattice_A, float knownLatticeTolerance,
                                       const accuracyConstants_t& accuracyConstants)
    {
        setDeterminantRange(determinantRange);
        setKnownLatticeParameters(sampleRealLattice_A, knownLatticeTolerance);
        setAccuracyConstants(accuracyConstants);
    }

    void LatticeAssembler::setStandardValues()
    {
        determinantRange << 0, numeric_limits<float>::max();

        latticeParametersKnown = false;
        knownLatticeParametersTolerance = 0;

        accuracyConstants.maxCountGlobalPassingWeightFilter = 500;
        accuracyConstants.maxCountLocalPassingWeightFilter = 15;
        accuracyConstants.maxCountPassingRelativeDefectFilter = 50;

        accuracyConstants.minPointsOnLattice = 5;

        accuracyConstants.maxCloseToPointDeviation = 0.15;
    }

    void LatticeAssembler::setKnownLatticeParameters(const Lattice& sampleRealLattice_A, float tolerance)
    {
        latticeParametersKnown = true;
        knownSampleRealLattice_A = sampleRealLattice_A;
        knownLatticeParameters << sampleRealLattice_A.getBasisVectorNorms(), sampleRealLattice_A.getBasisVectorAnglesNormalized_deg();
        knownLatticeParametersInverse = 1.0f / knownLatticeParameters;
        knownLatticeParametersTolerance = tolerance;
    }

    LatticeAssembler::accuracyConstants_t LatticeAssembler::getAccuracyConstants()
    {
        return accuracyConstants;
    }

    void LatticeAssembler::assembleLattices(vector<Lattice>& assembledLattices, Matrix3Xf& candidateVectors, RowVectorXf& candidateVectorWeights,
                                            vector<vector<uint16_t>>& pointIndicesOnVector, Matrix3Xf& pointsToFitInReciprocalSpace)
    {
        vector<assembledLatticeStatistics_t> assembledLatticesStatistics;
        assembleLattices(assembledLattices, assembledLatticesStatistics, candidateVectors, candidateVectorWeights, pointIndicesOnVector,
                         pointsToFitInReciprocalSpace);
    }

    void LatticeAssembler::assembleLattices(vector<Lattice>& assembledLattices, vector<assembledLatticeStatistics_t>& assembledLatticesStatistics,
                                            Matrix3Xf& candidateVectors, RowVectorXf& candidateVectorWeights, vector<vector<uint16_t>>& pointIndicesOnVector,
                                            Matrix3Xf& pointsToFitInReciprocalSpace)
    {
        reset();

        computeCandidateLattices(candidateVectors, candidateVectorWeights, pointIndicesOnVector);

        list<candidateLattice_t> finalCandidateLattices;

        filterCandidateLatticesByWeight(accuracyConstants.maxCountGlobalPassingWeightFilter);

        for (auto candidateLattice = candidateLattices.begin(); candidateLattice != candidateLattices.end(); ++candidateLattice)
        {
            candidateLattice->realSpaceLattice.minimize();
            candidateLattice->det = abs(candidateLattice->realSpaceLattice.det());
            computeAssembledLatticeStatistics(*candidateLattice, pointsToFitInReciprocalSpace);
        };

        // assume that candidateVectors is sorted descending for weight!
        finalCandidateLattices.insert(finalCandidateLattices.end(), candidateLattices.begin(),
                                      candidateLattices.begin() + min((uint32_t)candidateLattices.size(), accuracyConstants.maxCountLocalPassingWeightFilter));

        filterCandidateBasesByMeanRelativeDefect(accuracyConstants.maxCountPassingRelativeDefectFilter);

        finalCandidateLattices.insert(finalCandidateLattices.end(), candidateLattices.begin(), candidateLattices.end());

        selectBestLattices(assembledLattices, assembledLatticesStatistics, finalCandidateLattices);

        for (auto lattice = assembledLattices.begin(); lattice != assembledLattices.end(); ++lattice)
        {
            // cout << latticeParametersKnown << accuracyConstants.refineWithExactLattice;
            if (latticeParametersKnown)
            {
                if (accuracyConstants.refineWithExactLattice)
                {
                    lattice->reorder(knownSampleRealLattice_A);
                    refineLattice_peaksAndAngle_fixedBasisParameters(*lattice, pointsToFitInReciprocalSpace);
                }
                else
                {
                    refineLattice_peaksAndAngle(*lattice, pointsToFitInReciprocalSpace);
                }
            }
            else
            {
                refineLattice_peaksAndAngle(*lattice, pointsToFitInReciprocalSpace);

                lattice->normalizeAngles();
            }
        }
    }

    // clang-format off
void LatticeAssembler::selectBestLattices(vector< Lattice >& assembledLattices, vector< assembledLatticeStatistics_t >& assembledLatticesStatistics,
        list< candidateLattice_t >& finalCandidateLattices)
{
    assembledLattices.clear();
    if (finalCandidateLattices.size() == 0) {
        return;
    }

    finalCandidateLattices.sort([&](candidateLattice_t i, candidateLattice_t j) {return i.pointOnLatticeIndices.size() > j.pointOnLatticeIndices.size();}); //descending

    float significantDetReductionFactor = 0.75f;
    float significantPointCountReductionFactor = 0.85f;
    float significantMeanDefectReductionFactor = 0.7f;
    float significantMeanRelativeDefectReductionFactor = 0.8f;

    vector< uint16_t > combinedPointIndicesOnSelectedLattices;
    vector< uint16_t > tmp_indices;
    tmp_indices.reserve(4000); //just for speed
    combinedPointIndicesOnSelectedLattices.reserve(1000); //just for speed

    auto bestCandidateLattice = finalCandidateLattices.begin();
    auto nextCandidateLattice = finalCandidateLattices.begin();
    while (bestCandidateLattice != prev(finalCandidateLattices.end()) && bestCandidateLattice != finalCandidateLattices.end()) {
        tmp_indices.clear();
        set_union(combinedPointIndicesOnSelectedLattices.begin(), combinedPointIndicesOnSelectedLattices.end(),
                bestCandidateLattice->pointOnLatticeIndices.begin(), bestCandidateLattice->pointOnLatticeIndices.end(),
                back_inserter(tmp_indices));   // first resizing to take all elements and then rezizing to right size may be slightly faster than back-inserter
        combinedPointIndicesOnSelectedLattices.swap(tmp_indices);

        nextCandidateLattice = next(bestCandidateLattice);
        while (nextCandidateLattice != finalCandidateLattices.end()) {
            tmp_indices.clear();
            set_difference(nextCandidateLattice->pointOnLatticeIndices.begin(), nextCandidateLattice->pointOnLatticeIndices.end(),
                    combinedPointIndicesOnSelectedLattices.begin(), combinedPointIndicesOnSelectedLattices.end(),
                    back_inserter(tmp_indices));
            uint16_t uniquelyReachedNodesCount = tmp_indices.size();
            if (uniquelyReachedNodesCount >= accuracyConstants.minPointsOnLattice) { //enough new points on lattice => lattice cannot be rejected
                ++nextCandidateLattice;
            } else if (nextCandidateLattice->pointOnLatticeIndices.size()
                    > bestCandidateLattice->pointOnLatticeIndices.size() * significantPointCountReductionFactor) { //subset of previous lattice + zero to few points; Not significantly fewer points than current lattice
                if (
                (nextCandidateLattice->det <= bestCandidateLattice->det * significantDetReductionFactor
                        && nextCandidateLattice->assembledLatticeStatistics.meanDefect
                                * min(significantMeanDefectReductionFactor, 1.5f / (bestCandidateLattice->det / nextCandidateLattice->det))
                                < bestCandidateLattice->assembledLatticeStatistics.meanDefect
                        && nextCandidateLattice->assembledLatticeStatistics.meanRelativeDefect * significantMeanRelativeDefectReductionFactor
                                < bestCandidateLattice->assembledLatticeStatistics.meanRelativeDefect)
                        ||
                        (nextCandidateLattice->det * significantDetReductionFactor <= bestCandidateLattice->det
                                && nextCandidateLattice->assembledLatticeStatistics.meanDefect < bestCandidateLattice->assembledLatticeStatistics.meanDefect)
                        ||
                        (nextCandidateLattice->assembledLatticeStatistics.meanDefect
                                < bestCandidateLattice->assembledLatticeStatistics.meanDefect
                                        * min(significantMeanDefectReductionFactor, 1.5f / (nextCandidateLattice->det / bestCandidateLattice->det))
                                && nextCandidateLattice->assembledLatticeStatistics.meanRelativeDefect
                                        < bestCandidateLattice->assembledLatticeStatistics.meanRelativeDefect * significantMeanRelativeDefectReductionFactor)
                        ) {

                    auto temp = next(nextCandidateLattice);
                    finalCandidateLattices.splice(next(bestCandidateLattice), finalCandidateLattices, nextCandidateLattice);
                    bestCandidateLattice = finalCandidateLattices.erase(bestCandidateLattice);
                    nextCandidateLattice = temp;
                } else {
                    nextCandidateLattice = finalCandidateLattices.erase(nextCandidateLattice);
                }
            } else {
                nextCandidateLattice = finalCandidateLattices.erase(nextCandidateLattice);
            }
        }
        ++bestCandidateLattice;
    }

    assembledLattices.reserve(finalCandidateLattices.size());
    for (auto selectedCandidateLattice = finalCandidateLattices.begin(); selectedCandidateLattice != finalCandidateLattices.end(); 
		++selectedCandidateLattice) {
        assembledLattices.push_back(selectedCandidateLattice->realSpaceLattice);
        assembledLatticesStatistics.push_back(selectedCandidateLattice->assembledLatticeStatistics);
    }

}
    // clang-format on

    uint16_t LatticeAssembler::countUniqueColumns(const Matrix3Xf& millerIndices)
    {
        uint16_t count = 0;

        sortIndices.resize(millerIndices.cols());
        iota(sortIndices.begin(), sortIndices.end(), 0);
        sort(sortIndices.begin(), sortIndices.end(), [&millerIndices](uint16_t i, uint16_t j) { return millerIndices(0, i) < millerIndices(0, j); });

        Matrix3Xf millerIndices_sorted(3, millerIndices.cols() + 1);
        for (int i = 0; i < millerIndices.cols(); ++i)
        {
            millerIndices_sorted.col(i) = millerIndices.col(sortIndices[i]);
        }
        millerIndices_sorted(0, millerIndices.cols()) = 0.5f; // not a valid miller index! This saves an out of bounds check in a later loop

        for (int i = 0; i < millerIndices.cols(); ++i)
        {
            bool isUnique = true;
            for (int j = i + 1; millerIndices_sorted(0, i) == millerIndices_sorted(0, j); j++)
            {
                if ((millerIndices_sorted(1, i) == millerIndices_sorted(1, j)) & (millerIndices_sorted(2, i) == millerIndices_sorted(2, j)))
                {
                    isUnique = false;
                    break;
                }
            }
            if (isUnique)
            {
                count++;
            }
        }

        return count;
    }

    void LatticeAssembler::computeCandidateLattices(Matrix3Xf& candidateVectors, RowVectorXf& candidateVectorWeights,
                                                    vector<vector<uint16_t>>& pointIndicesOnVector)
    {
        // hand-crafted remove-if for three variables.
        for (int i = candidateVectors.cols() - 1; i >= 0; i--)
        {
            if (pointIndicesOnVector[i].size() < accuracyConstants.minPointsOnLattice)
            {
                pointIndicesOnVector[i] = pointIndicesOnVector.back();
                pointIndicesOnVector.pop_back();

                candidateVectors.col(i) = candidateVectors.col(candidateVectors.cols() - 1);
                candidateVectors.conservativeResize(NoChange,
                                                    candidateVectors.cols() - 1); // possibly slow, because of call to realloc. can be avoided if needed

                candidateVectorWeights[i] = candidateVectorWeights[candidateVectorWeights.size() - 1];
                candidateVectorWeights.conservativeResize(candidateVectorWeights.size() - 1);
            }
        }

        // needed for easier intersection computation later
        const auto end = pointIndicesOnVector.end();
        for (auto it = pointIndicesOnVector.begin(); it < end; ++it)
        {
            sort(it->begin(), it->end());
        }

        candidateLattices.reserve(10000);
        int candidateVectorsCount = candidateVectors.cols();
        for (uint16_t i = 0; i < candidateVectorsCount - 2; ++i)
        {
            for (uint16_t j = (i + 1); j < candidateVectorsCount - 1; ++j)
            {
                for (uint16_t k = (j + 1); k < candidateVectorsCount; ++k)
                {
                    Lattice latticeToCheck(candidateVectors.col(i), candidateVectors.col(j), candidateVectors.col(k));

                    float absDet = abs(latticeToCheck.det());
                    if ((absDet < determinantRange[0]) | (absDet > determinantRange[1]))
                    {
                        continue;
                    }

                    vector<uint16_t> pointIndicesOnTwoVectors(max(pointIndicesOnVector[i].size(), pointIndicesOnVector[j].size()));
                    auto it = set_intersection(pointIndicesOnVector[i].begin(), pointIndicesOnVector[i].end(), pointIndicesOnVector[j].begin(),
                                               pointIndicesOnVector[j].end(), pointIndicesOnTwoVectors.begin());
                    uint16_t pointsOnBothVectorsCount = it - pointIndicesOnTwoVectors.begin();
                    if (pointsOnBothVectorsCount < accuracyConstants.minPointsOnLattice)
                    {
                        continue;
                    }
                    pointIndicesOnTwoVectors.resize(pointsOnBothVectorsCount);

                    vector<uint16_t> pointIndicesOnLatticeToCheck(max(pointIndicesOnTwoVectors.size(), pointIndicesOnVector[k].size()));
                    it = set_intersection(pointIndicesOnTwoVectors.begin(), pointIndicesOnTwoVectors.end(), pointIndicesOnVector[k].begin(),
                                          pointIndicesOnVector[k].end(), pointIndicesOnLatticeToCheck.begin());
                    uint16_t pointsOnLatticeToCheckCount = it - pointIndicesOnLatticeToCheck.begin();
                    if (pointsOnLatticeToCheckCount < accuracyConstants.minPointsOnLattice)
                    {
                        continue;
                    }
                    pointIndicesOnLatticeToCheck.resize(pointsOnLatticeToCheckCount);

                    if (latticeParametersKnown)
                    {
                        if (!checkLatticeParameters(latticeToCheck))
                        {
                            continue;
                        }
                    }

                    candidateLattices.resize(candidateLattices.size() + 1);
                    auto& newCandidateBasis = candidateLattices.back();
                    newCandidateBasis.realSpaceLattice = latticeToCheck;
                    newCandidateBasis.weight = candidateVectorWeights[i] + candidateVectorWeights[j] + candidateVectorWeights[k];
                    newCandidateBasis.pointOnLatticeIndices.swap(pointIndicesOnLatticeToCheck);
                    newCandidateBasis.vectorIndices = {i, j, k};
                }
            }
        }
    }

    // clang-format off
void LatticeAssembler::computeAssembledLatticeStatistics(candidateLattice_t& candidateLattice, const Matrix3Xf& pointsToFitInReciprocalSpace)
{
    auto& pointOnLatticeIndices = candidateLattice.pointOnLatticeIndices;
    Matrix3Xf currentPointsToFitInReciprocalSpace(3, pointOnLatticeIndices.size());
    for (uint32_t j = 0; j < pointOnLatticeIndices.size(); ++j) {
        currentPointsToFitInReciprocalSpace.col(j) = pointsToFitInReciprocalSpace.col(pointOnLatticeIndices[j]);
    }

//    auto reciprocalBasis = candidateLattice.realSpaceLattice.getBasis().transpose().inverse();
//    Matrix3Xf factorsToReachPoints = reciprocalBasis.inverse()*pointsToFitInInverseSpace;
    Matrix3Xf factorsToReachPoints = candidateLattice.realSpaceLattice.getBasis().transpose() * currentPointsToFitInReciprocalSpace; //realSpaceLattice is inverse of the transpose of the reciprocal basis. Inverse of the reciprocal basis is needed => transpose!

    Matrix3Xf millerIndices = factorsToReachPoints.array().round();

    candidateLattice.assembledLatticeStatistics.occupiedLatticePointsCount = countUniqueColumns(millerIndices);

    auto predictedPoints = candidateLattice.realSpaceLattice.getBasis().transpose().inverse() * millerIndices;

//    cout << predictedPoints << endl << endl;
//    cout << (predictedPoints - currentPointsToFitInReciprocalSpace).eval() << endl << endl;
//    cout << ((predictedPoints - currentPointsToFitInReciprocalSpace).colwise().norm()).eval() << endl << endl;

    candidateLattice.assembledLatticeStatistics.meanDefect = ((predictedPoints - currentPointsToFitInReciprocalSpace).colwise().norm()).mean();
    candidateLattice.assembledLatticeStatistics.meanRelativeDefect = ((factorsToReachPoints - millerIndices).colwise().norm()).mean();
}
    // clang-format on

    void LatticeAssembler::filterCandidateLatticesByWeight(uint32_t maxToTakeCount)
    {
        uint32_t toTakeCount = min(maxToTakeCount, (uint32_t)candidateLattices.size());

        sortIndices.resize(candidateLattices.size());
        iota(sortIndices.begin(), sortIndices.end(), 0);
        nth_element(sortIndices.begin(), sortIndices.begin() + toTakeCount, sortIndices.end(),
                    [&](uint32_t i, uint32_t j) { return candidateLattices[i].weight > candidateLattices[j].weight; }); // descending

        sortIndices.resize(toTakeCount);
        sort(sortIndices.begin(), sortIndices.end(),
             [&](uint32_t i, uint32_t j) { return candidateLattices[i].weight > candidateLattices[j].weight; }); // descending

        vector<candidateLattice_t> candidateLatticesFiltered;
        candidateLatticesFiltered.resize(toTakeCount);
        for (uint32_t i = 0; i < toTakeCount; ++i)
        {
            candidateLatticesFiltered[i] = candidateLattices[sortIndices[i]];
        }

        candidateLattices.swap(candidateLatticesFiltered);
    }


    bool LatticeAssembler::checkLatticeParameters(Lattice& lattice)
    {
        lattice.minimize();

        Vector3f n = lattice.getBasisVectorNorms();
        Vector3f a = lattice.getBasisVectorAnglesNormalized_deg();

        Array<float, 6, 6> allPermutations;

        // clang-format off
    allPermutations <<
            n[0], n[0], n[1], n[1], n[2], n[2],
            n[1], n[2], n[0], n[2], n[0], n[1],
            n[2], n[1], n[2], n[0], n[1], n[0],
            a[0], a[0], a[1], a[1], a[2], a[2],
            a[1], a[2], a[0], a[2], a[0], a[1],
            a[2], a[1], a[2], a[0], a[1], a[0];
        // clang-format on

        auto rasiduals = ((allPermutations.colwise() - knownLatticeParameters).colwise() * knownLatticeParametersInverse).abs(); // Array< bool, 6, 6 >

        auto parametersValid = rasiduals < knownLatticeParametersTolerance; // Array< bool, 6, 6 >

        auto permutationsValid = parametersValid.colwise().all(); // Array< bool, 1, 6 >

        bool latticeValid = permutationsValid.any();

        return latticeValid;
    }


    void LatticeAssembler::filterCandidateBasesByMeanRelativeDefect(uint32_t maxToTakeCount)
    {
        uint32_t toTakeCount = min(maxToTakeCount, (uint32_t)candidateLattices.size());

        sortIndices.resize(candidateLattices.size());
        iota(sortIndices.begin(), sortIndices.end(), 0);
        nth_element(sortIndices.begin(), sortIndices.begin() + toTakeCount, sortIndices.end(), [&](uint32_t i, uint32_t j) {
            return candidateLattices[i].assembledLatticeStatistics.meanRelativeDefect < candidateLattices[j].assembledLatticeStatistics.meanRelativeDefect;
        });

        sortIndices.resize(toTakeCount);
        sort(sortIndices.begin(), sortIndices.end(), [&](uint32_t i, uint32_t j) {
            return candidateLattices[i].assembledLatticeStatistics.meanRelativeDefect < candidateLattices[j].assembledLatticeStatistics.meanRelativeDefect;
        });

        vector<candidateLattice_t> candidateLatticesFiltered;
        candidateLatticesFiltered.resize(toTakeCount);
        for (uint32_t i = 0; i < toTakeCount; ++i)
        {
            candidateLatticesFiltered[i] = candidateLattices[sortIndices[i]];
        }

        candidateLattices.swap(candidateLatticesFiltered);
    }

    void LatticeAssembler::setAccuracyConstants(const accuracyConstants_t& accuracyConstants)
    {
        this->accuracyConstants = accuracyConstants;
    }

    void LatticeAssembler::setDeterminantRange(const Eigen::Vector2f& determinantRange)
    {
        this->determinantRange = determinantRange;
    }

    void LatticeAssembler::setDeterminantRange(float min, float max)
    {
        this->determinantRange = Vector2f(min, max);
    }

    void LatticeAssembler::reset()
    {
        candidateLattices.clear();
        validLattices.clear();
    }

//#define MEAN_SQUARED_DIST_REFINE
#define MEAN_DIST_REFINE
    static void keepGoodReciprocalPeaks(Matrix3Xf& keptPeaks, Matrix3Xf& keptMillerIndices, vector<uint16_t>& pointOnLatticeIndices,
                                        const Array<bool, 1, Dynamic>& goodPeaksFlags, const Matrix3Xf& allPeaks, const Matrix3Xf& allMillerIndices);
    void LatticeAssembler::refineLattice(Lattice& realSpaceLattice, vector<uint16_t>& pointOnLatticeIndices, const Matrix3Xf& pointsToFitInReciprocalSpace)
    {
        Lattice& bestLattice = realSpaceLattice;

        int maxIterationsCount = 5;
        int iterationCount = 0;

        Matrix3Xf factorsToReachNodes;
        Matrix3Xf millerIndices;
        Matrix3f gradient;
        Array<float, 1, Dynamic> maxRelativeDefects;
        Array<bool, 1, Dynamic> goodReciprocalPeaksFlags;
        Matrix3Xf reciprocalPeaksUsedForFitting_1_per_A;
        Matrix3Xf millerIndicesUsedForFitting;
        while (1)
        {
            factorsToReachNodes = bestLattice.getBasis().transpose() * pointsToFitInReciprocalSpace;
            millerIndices = factorsToReachNodes.array().round();

            maxRelativeDefects = (millerIndices - factorsToReachNodes).array().abs().colwise().maxCoeff();

            // tune!
            goodReciprocalPeaksFlags = maxRelativeDefects < accuracyConstants.maxCloseToPointDeviation;

            int goodReciprocalPeaksCount = goodReciprocalPeaksFlags.cast<uint16_t>().sum();
            if (goodReciprocalPeaksCount < 5) // if unstable
                break;

            reciprocalPeaksUsedForFitting_1_per_A.resize(3, goodReciprocalPeaksCount);
            millerIndicesUsedForFitting.resize(3, goodReciprocalPeaksCount);
            keepGoodReciprocalPeaks(reciprocalPeaksUsedForFitting_1_per_A, millerIndicesUsedForFitting, pointOnLatticeIndices, goodReciprocalPeaksFlags,
                                    pointsToFitInReciprocalSpace, millerIndices);

#ifdef MEAN_DIST_REFINE
            // gradient descent
            Matrix3f refinedReciprocalBasis = bestLattice.getReciprocalLattice().getBasis();
            float stepLength = refinedReciprocalBasis.maxCoeff() * 0.002;
            for (int i = 0; i < 50; i++)
            {
                getGradient_reciprocalPeakMatch_meanDist(gradient, refinedReciprocalBasis, millerIndicesUsedForFitting, pointsToFitInReciprocalSpace);
                float maxCoeff = gradient.cwiseAbs().maxCoeff();
                if (maxCoeff < 1e-20)
                {
                    break;
                }
                gradient /= maxCoeff;

                refinedReciprocalBasis = refinedReciprocalBasis - stepLength * gradient;

                if (i >= 25)
                {
                    stepLength *= 0.87;
                }
            }
#elif defined MEAN_SQUARED_DIST_REFINE
            Matrix3f refinedReciprocalBasis;
            refineReciprocalBasis_meanSquaredDist(refinedReciprocalBasis, factorsToReachNodes, pointsToFitInReciprocalSpace);
#endif

            Lattice refinedLattice = Lattice(refinedReciprocalBasis).getReciprocalLattice();
            refinedLattice.minimize();

            if ((refinedLattice.getBasis() - bestLattice.getBasis()).isZero(1e-3))
            {
                return;
            }
            else if (iterationCount++ > maxIterationsCount)
            {
                bestLattice = refinedLattice;
                return;
            }
            else
            {
                bestLattice = refinedLattice;
            }
        }
    }


    void LatticeAssembler::refineLattice_peaksAndAngle(Lattice& realSpaceLattice, const Matrix3Xf& pointsToFitInReciprocalSpace)
    {
        Lattice& bestLattice = realSpaceLattice;

        int maxIterationsCount = 5;
        int iterationCount = 0;

        Matrix3Xf factorsToReachNodes;
        Matrix3Xf millerIndices;
        Matrix3f gradient;
        Array<float, 1, Dynamic> maxRelativeDefects;
        Array<bool, 1, Dynamic> goodReciprocalPeaksFlags;
        Matrix3Xf reciprocalPeaksUsedForFitting_1_per_A;
        Matrix3Xf millerIndicesUsedForFitting;
        vector<uint16_t> pointOnLatticeIndices_junk;
        while (1)
        {
            factorsToReachNodes = bestLattice.getBasis().transpose() * pointsToFitInReciprocalSpace;
            millerIndices = factorsToReachNodes.array().round();

            maxRelativeDefects = (millerIndices - factorsToReachNodes).array().abs().colwise().maxCoeff();

            // tune!
            goodReciprocalPeaksFlags = maxRelativeDefects < accuracyConstants.maxCloseToPointDeviation;

            int goodReciprocalPeaksCount = goodReciprocalPeaksFlags.cast<uint16_t>().sum();
            if (goodReciprocalPeaksCount < 5) // if unstable
                break;

            reciprocalPeaksUsedForFitting_1_per_A.resize(3, goodReciprocalPeaksCount);
            millerIndicesUsedForFitting.resize(3, goodReciprocalPeaksCount);
            keepGoodReciprocalPeaks(reciprocalPeaksUsedForFitting_1_per_A, millerIndicesUsedForFitting, pointOnLatticeIndices_junk, goodReciprocalPeaksFlags,
                                    pointsToFitInReciprocalSpace, millerIndices);

            Matrix3f refinedReciprocalBasis = bestLattice.getReciprocalLattice().getBasis();
            refineReciprocalBasis_meanDist_peaksAndAngle(refinedReciprocalBasis, millerIndicesUsedForFitting, reciprocalPeaksUsedForFitting_1_per_A);

            Lattice refinedLattice = Lattice(refinedReciprocalBasis).getReciprocalLattice();
            refinedLattice.minimize();

            if ((refinedLattice.getBasis() - bestLattice.getBasis()).isZero(1e-3))
            {
                bestLattice.minimize();
                return;
            }
            else if (iterationCount++ > maxIterationsCount)
            {
                bestLattice = refinedLattice;
                bestLattice.minimize();
                return;
            }
            else
            {
                bestLattice = refinedLattice;
            }
        }
        bestLattice.minimize();
    }

    void LatticeAssembler::refineLattice_peaksAndAngle_fixedBasisParameters(Lattice& realSpaceLattice, const Matrix3Xf& pointsToFitInReciprocalSpace)
    {
        Lattice& bestLattice = realSpaceLattice;

        int maxIterationsCount = 5;
        int iterationCount = 0;

        Matrix3Xf factorsToReachNodes;
        Matrix3Xf millerIndices;
        Matrix3f gradient;
        Array<float, 1, Dynamic> maxRelativeDefects;
        Array<bool, 1, Dynamic> goodReciprocalPeaksFlags;
        Matrix3Xf reciprocalPeaksUsedForFitting_1_per_A;
        Matrix3Xf millerIndicesUsedForFitting;
        vector<uint16_t> pointOnLatticeIndices_junk;
        while (1)
        {
            factorsToReachNodes = bestLattice.getBasis().transpose() * pointsToFitInReciprocalSpace;
            millerIndices = factorsToReachNodes.array().round();

            maxRelativeDefects = (millerIndices - factorsToReachNodes).array().abs().colwise().maxCoeff();

            // tune!
            goodReciprocalPeaksFlags = maxRelativeDefects < accuracyConstants.maxCloseToPointDeviation;

            int goodReciprocalPeaksCount = goodReciprocalPeaksFlags.cast<uint16_t>().sum();
            if (goodReciprocalPeaksCount < 5) // if unstable
                break;

            reciprocalPeaksUsedForFitting_1_per_A.resize(3, goodReciprocalPeaksCount);
            millerIndicesUsedForFitting.resize(3, goodReciprocalPeaksCount);
            keepGoodReciprocalPeaks(reciprocalPeaksUsedForFitting_1_per_A, millerIndicesUsedForFitting, pointOnLatticeIndices_junk, goodReciprocalPeaksFlags,
                                    pointsToFitInReciprocalSpace, millerIndices);

            Matrix3f refinedReciprocalBasis = bestLattice.getReciprocalLattice().getBasis();
            refineReciprocalBasis_meanSquaredDist_fixedBasisParameters(refinedReciprocalBasis, millerIndicesUsedForFitting,
                                                                       reciprocalPeaksUsedForFitting_1_per_A,
                                                                       knownSampleRealLattice_A.getReciprocalLattice().getBasis());
            refineReciprocalBasis_meanDist_detectorAngleMatchFixedParameters(refinedReciprocalBasis, millerIndicesUsedForFitting,
                                                                             reciprocalPeaksUsedForFitting_1_per_A);

            Lattice refinedLattice = Lattice(refinedReciprocalBasis).getReciprocalLattice();

            if ((refinedLattice.getBasis() - bestLattice.getBasis()).isZero(1e-3))
            {
                return;
            }
            else if (iterationCount++ > maxIterationsCount)
            {
                bestLattice = refinedLattice;
                return;
            }
            else
            {
                bestLattice = refinedLattice;
            }
        }
    }

    static void keepGoodReciprocalPeaks(Matrix3Xf& keptPeaks, Matrix3Xf& keptMillerIndices, vector<uint16_t>& pointOnLatticeIndices,
                                        const Array<bool, 1, Dynamic>& goodPeaksFlags, const Matrix3Xf& allPeaks, const Matrix3Xf& allMillerIndices)
    {
        pointOnLatticeIndices.clear();

        int peeksKeptCount = 0;
        for (int i = 0; i < goodPeaksFlags.size(); i++)
        {
            if (goodPeaksFlags[i])
            {
                keptPeaks.col(peeksKeptCount) = allPeaks.col(i);
                keptMillerIndices.col(peeksKeptCount) = allMillerIndices.col(i);
                pointOnLatticeIndices.push_back(i);
                peeksKeptCount++;
            }
        }
    }
} // namespace xgandalf
