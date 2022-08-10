/*
 * IndexerAutocorrPrefit.h
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

#ifndef INDEXERAUTOCORRPREFIT_H_
#define INDEXERAUTOCORRPREFIT_H_

#include "HillClimbingOptimizer.h"
#include <xgandalf/IndexerBase.h>

namespace xgandalf
{
    class IndexerAutocorrPrefit : public IndexerBase
    {
      public:
        enum class SamplingPitch
        {
            extremelyLoose,
            loose,
            standard,
            dense,
            extremelyDense
        };

        IndexerAutocorrPrefit(const ExperimentSettings& experimentSettings);

        void index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A);

        void setSamplingPitch(SamplingPitch samplingPitch);
        void setSamplingPitch(float unitPitch);

      private:
        void precompute();

        void getGoodAutocorrelationPoints(Eigen::Matrix3Xf& goodAutocorrelationPoints, Eigen::RowVectorXf& goodAutocorrelationPointWeights,
                                          const Eigen::Matrix3Xf& points, uint32_t maxAutocorrelationPointsCount);
        void autocorrPrefit(const Eigen::Matrix3Xf& reciprocalPeaks_A, Eigen::Matrix3Xf& samplePoints,
                            HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_autocorr);

        Eigen::Matrix3Xf precomputedSamplePoints;

        HillClimbingOptimizer hillClimbingOptimizer;
        SparsePeakFinder sparsePeakFinder;
        InverseSpaceTransform inverseSpaceTransform;

        float maxCloseToPointDeviation;
        float maxNormInAutocorrelation;
        float minNormInAutocorrelation;
        float dbscanEpsilon;
        Dbscan dbscan;
    };
} // namespace xgandalf
#endif /* INDEXERAUTOCORRPREFIT_H_ */
