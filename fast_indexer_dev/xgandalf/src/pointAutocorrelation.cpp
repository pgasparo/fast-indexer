/*
 * pointAutocorrelation.cpp
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

#include <ctype.h>
#include <xgandalf/pointAutocorrelation.h>

namespace xgandalf
{
    using namespace std;
    using namespace Eigen;

    void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, const Matrix3Xf& points, float minNormInAutocorrelation, float maxNormInAutocorrelation)
    {
        uint32_t N = points.cols();
        uint32_t n = N - 1;
        uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

        float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;
        float minNormInAutocorrelation_squared = minNormInAutocorrelation * minNormInAutocorrelation;

        autocorrelationPoints.resize(3, maxAutocorrPointsCount);

        uint32_t autocorrelationPointsCount = 0;

        for (uint32_t i = 0; i < N; ++i)
        {
            uint32_t addedNodesInIteration = n - i;

            autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) =
                points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

            for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k)
            {
                if (autocorrelationPoints(2, j) < 0)
                {
                    autocorrelationPoints.col(j) *= -1;
                }
            }

            uint32_t endInThisIteration = autocorrelationPointsCount + addedNodesInIteration;
            for (uint32_t j = autocorrelationPointsCount, k = 0; j < endInThisIteration; ++j, ++k)
            {
                float squaredNorm = autocorrelationPoints.col(j).squaredNorm();
                if (squaredNorm < maxNormInAutocorrelation_squared && squaredNorm > minNormInAutocorrelation_squared)
                {
                    autocorrelationPoints.col(autocorrelationPointsCount) = autocorrelationPoints.col(j);
                    autocorrelationPointsCount++;
                }
            }
        }

        autocorrelationPoints.conservativeResize(3, autocorrelationPointsCount);
    }

    void getPointAutocorrelation(Matrix3Xf& autocorrelationPoints, VectorXi& centerPointIndices, VectorXi& shiftedPointIndices, const Matrix3Xf& points,
                                 float minNormInAutocorrelation, float maxNormInAutocorrelation)
    {
        uint32_t N = points.cols();
        uint32_t n = N - 1;
        uint32_t maxAutocorrPointsCount = n * (n + 1) / 2;

        float maxNormInAutocorrelation_squared = maxNormInAutocorrelation * maxNormInAutocorrelation;
        float minNormInAutocorrelation_squared = minNormInAutocorrelation * minNormInAutocorrelation;

        autocorrelationPoints.resize(3, maxAutocorrPointsCount);
        centerPointIndices.resize(maxAutocorrPointsCount);
        shiftedPointIndices.resize(maxAutocorrPointsCount);

        uint32_t autocorrelationPointsCount = 0;

        for (uint32_t i = 0; i < N; ++i)
        {
            uint32_t addedNodesInIteration = n - i;

            autocorrelationPoints.block(0, autocorrelationPointsCount, 3, addedNodesInIteration) =
                points.block(0, i + 1, 3, N - i - 1).colwise() - points.col(i);

            for (uint32_t j = autocorrelationPointsCount, k = i + 1; j < autocorrelationPointsCount + addedNodesInIteration; ++j, ++k)
            {
                if (autocorrelationPoints(2, j) >= 0)
                {
                    centerPointIndices[j] = i;
                    shiftedPointIndices[j] = k;
                }
                else
                {
                    autocorrelationPoints.col(j) *= -1;
                    centerPointIndices[j] = k;
                    shiftedPointIndices[j] = i;
                }
            }

            uint32_t endInThisIteration = autocorrelationPointsCount + addedNodesInIteration;
            for (uint32_t j = autocorrelationPointsCount, k = 0; j < endInThisIteration; ++j, ++k)
            {
                float squaredNorm = autocorrelationPoints.col(j).squaredNorm();
                if (squaredNorm < maxNormInAutocorrelation_squared && squaredNorm > minNormInAutocorrelation_squared)
                {
                    autocorrelationPoints.col(autocorrelationPointsCount) = autocorrelationPoints.col(j);
                    centerPointIndices[autocorrelationPointsCount] = centerPointIndices[j];
                    shiftedPointIndices[autocorrelationPointsCount] = shiftedPointIndices[j];
                    autocorrelationPointsCount++;
                }
            }
        }

        autocorrelationPoints.conservativeResize(3, autocorrelationPointsCount);
        centerPointIndices.conservativeResize(autocorrelationPointsCount);
        shiftedPointIndices.conservativeResize(autocorrelationPointsCount);
    }
} // namespace xgandalf
