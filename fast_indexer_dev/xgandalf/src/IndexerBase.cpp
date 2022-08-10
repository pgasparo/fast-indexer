/*
 * IndexerBase.cpp
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

#include <xgandalf/IndexerBase.h>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <vector>

using namespace Eigen;
using namespace std;

namespace xgandalf
{
    IndexerBase::IndexerBase(const ExperimentSettings& experimentSettings)
        : experimentSettings(experimentSettings)
    {
    }

    void IndexerBase::keepSamplePointsWithHighEvaluation(Matrix3Xf& samplePoints, RowVectorXf& samplePointsEvaluation, float minEvaluation)
    {
        uint32_t bigEvaluationSamplePointsCount = 0;
        for (int i = 0; i < samplePointsEvaluation.size(); ++i)
        {
            if (samplePointsEvaluation[i] >= minEvaluation)
            {
                samplePoints.col(bigEvaluationSamplePointsCount) = samplePoints.col(i);
                samplePointsEvaluation[bigEvaluationSamplePointsCount] = samplePointsEvaluation[i];
                bigEvaluationSamplePointsCount++;
            }
        }
        samplePoints.conservativeResize(3, bigEvaluationSamplePointsCount);
        samplePointsEvaluation.conservativeResize(bigEvaluationSamplePointsCount);
    }

    void IndexerBase::keepSamplePointsWithHighestEvaluation(Eigen::Matrix3Xf& samplePoints, RowVectorXf& samplePointsEvaluation, uint32_t maxToTakeCount)
    {
        uint32_t toTakeCount = min(maxToTakeCount, (uint32_t)samplePointsEvaluation.size());

        sortIndices.resize(samplePointsEvaluation.size());
        iota(sortIndices.begin(), sortIndices.end(), 0);
        nth_element(sortIndices.begin(), sortIndices.begin() + toTakeCount, sortIndices.end(),
                    [&](uint32_t i, uint32_t j) { return samplePointsEvaluation[i] > samplePointsEvaluation[j]; });

        sortIndices.resize(toTakeCount);
        sort(sortIndices.begin(), sortIndices.end(), [&](uint32_t i, uint32_t j) { return samplePointsEvaluation[i] > samplePointsEvaluation[j]; });

        Matrix3Xf samplePoints_filtered(3, toTakeCount);
        RowVectorXf samplePointsEvaluation_filtered(toTakeCount);
        for (uint32_t i = 0; i < toTakeCount; ++i)
        {
            samplePoints_filtered.col(i) = samplePoints.col(sortIndices[i]);
            samplePointsEvaluation_filtered[i] = samplePointsEvaluation[sortIndices[i]];
        }

        samplePoints.swap(samplePoints_filtered);
        samplePointsEvaluation.swap(samplePointsEvaluation_filtered);
    }
} // namespace xgandalf
