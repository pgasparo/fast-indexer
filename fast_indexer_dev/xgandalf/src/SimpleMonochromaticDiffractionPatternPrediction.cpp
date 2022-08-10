/*
 * samplePointsGenerator.cpp
 *
 * SimpleMonochromaticDiffractionPatternPrediction.cpp
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

#include <xgandalf/SimpleMonochromaticDiffractionPatternPrediction.h>
#include <xgandalf/eigenSTLContainers.h>

#include <cmath>

namespace xgandalf
{
    using namespace std;
    using namespace Eigen;

    SimpleMonochromaticDiffractionPatternPrediction::SimpleMonochromaticDiffractionPatternPrediction(const ExperimentSettings& experimentSettings)
        : reciprocalToRealProjection(experimentSettings)
    {
        maxResolutionAngle = experimentSettings.getMaxResolutionAngle_rad();
        reflectionRadius = experimentSettings.getReflectionRadius();
        reciprocalLambdaShort = experimentSettings.getReciprocalLambdaShort_1A();
        reciprocalLambdaLong = experimentSettings.getReciprocalLambdaLong_1A();
        reciprocalLambdaShort_extended_squared = (reciprocalLambdaShort + reflectionRadius) * (reciprocalLambdaShort + reflectionRadius);
        reciprocalLambdaLong_extended_squared = (reciprocalLambdaLong - reflectionRadius) * (reciprocalLambdaLong - reflectionRadius);
        detectorDistance = experimentSettings.getDetectorDistance_m();
    }

    void SimpleMonochromaticDiffractionPatternPrediction::predictPattern(Matrix2Xf& predictedPeaks, Matrix3Xi& millerIndices, Matrix3Xf& projectionDirections,
                                                                         const Lattice& lattice)
    {
        Matrix3Xf peaksOnEwaldSphere;
        getPeaksOnEwaldSphere(peaksOnEwaldSphere, millerIndices, lattice);

        reciprocalToRealProjection.project(predictedPeaks, peaksOnEwaldSphere);

        projectionDirections.resize(3, predictedPeaks.cols());
        projectionDirections.row(0).setConstant(detectorDistance);
        projectionDirections.bottomRows(2) = predictedPeaks;
        projectionDirections.colwise().normalize();
    }

    void SimpleMonochromaticDiffractionPatternPrediction::getPeaksOnEwaldSphere(Matrix3Xf& peaksOnEwaldSphere, Matrix3Xi& millerIndices, const Lattice& lattice)
    {
        EigenSTL::vector_Vector3f peaks;
        EigenSTL::vector_Vector3i millers;

        Matrix3f basis = lattice.getBasis();

        float maxYZ = sin(maxResolutionAngle) * reciprocalLambdaShort + reflectionRadius;
        float maxYZ_squared = maxYZ * maxYZ;
        float maxX = reciprocalLambdaShort * (1 - cos(maxResolutionAngle)) + reflectionRadius;
        float maxX_pos = maxX * 0.05;
        Matrix<float, 3, 8> maxPointsNeededToReach;
        // clang-format off
	maxPointsNeededToReach << maxX_pos, maxX_pos, maxX_pos, maxX_pos, -maxX,  -maxX,  -maxX,  -maxX,
								maxYZ,   -maxYZ,    maxYZ,   -maxYZ,   maxYZ, -maxYZ,  maxYZ, -maxYZ,
								maxYZ,    maxYZ,   -maxYZ,   -maxYZ,   maxYZ,  maxYZ, -maxYZ, -maxYZ;
        // clang-format on
        Matrix<float, 3, 8> millerIndicesNeededToReachMaxPoints = basis.inverse() * maxPointsNeededToReach;

        Array3f millerIndicesBounds_min, millerIndicesBounds_max;
        millerIndicesBounds_min = millerIndicesNeededToReachMaxPoints.array().rowwise().minCoeff().floor();
        millerIndicesBounds_max = millerIndicesNeededToReachMaxPoints.array().rowwise().maxCoeff().ceil();

        Vector3f peakH, peakHK, peakHKL;
        for (int h = millerIndicesBounds_min[0]; h < millerIndicesBounds_max[0]; h++)
        {
            peakH = basis.col(0) * h;
            for (int k = millerIndicesBounds_min[1]; k < millerIndicesBounds_max[1]; k++)
            {
                peakHK = peakH + basis.col(1) * k;
                for (int l = millerIndicesBounds_min[2]; l < millerIndicesBounds_max[2]; l++)
                {
                    peakHKL = peakHK + basis.col(2) * l;

                    // clang-format off
                if (peakHKL[0] < maxX_pos && 
					peakHKL[0] > -maxX && 
					peakHKL[1] < maxYZ && 
					peakHKL[1] > -maxYZ && 
					peakHKL[2] < maxYZ && 
					peakHKL[2] > -maxYZ &&
                    peakHKL.tail(2).squaredNorm() < maxYZ_squared
					) // clang-format on
                    {
                        Vector3f centerToBorder = peakHKL;
                        centerToBorder[0] += reciprocalLambdaLong;
                        float tailNorm = centerToBorder.tail(2).squaredNorm();
                        if (centerToBorder[0] * centerToBorder[0] + tailNorm > reciprocalLambdaLong_extended_squared)
                        {
                            centerToBorder[0] = peakHKL[0] + reciprocalLambdaShort;
                            if (centerToBorder[0] * centerToBorder[0] + tailNorm < reciprocalLambdaShort_extended_squared)
                            {
                                if (!peakHKL.isZero(0))
                                {
                                    peaks.push_back(peakHKL);
                                    millers.emplace_back(h, k, l);
                                }
                            }
                        }
                    }
                }
            }
        }

        peaksOnEwaldSphere = Map<Matrix3Xf>((float*)peaks.data(), 3, peaks.size());
        millerIndices = Map<Matrix3Xi>((int*)millers.data(), 3, millers.size());
    }
} // namespace xgandalf
