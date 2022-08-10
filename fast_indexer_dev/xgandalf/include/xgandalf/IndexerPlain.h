/*
 * IndexerPlain.h
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

#ifndef INDEXERPLAIN_H_
#define INDEXERPLAIN_H_

#include <xgandalf/HillClimbingOptimizer.h>
#include <xgandalf/IndexerBase.h>

namespace xgandalf
{
    class IndexerPlain : public IndexerBase
    {
      public:
        enum class SamplingPitch
        {
            extremelyLoose,
            loose,
            standard,
            dense,
            extremelyDense,

            standardWithSeondaryMillerIndices,
            denseWithSeondaryMillerIndices,
            extremelyDenseWithSeondaryMillerIndices
        };

        enum class GradientDescentIterationsCount
        {
            exremelyFew,
            few,
            standard,
            many,
            manyMany,
            extremelyMany,

            custom
        };

        IndexerPlain(const ExperimentSettings& experimentSettings);

        void index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A);
        void index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A, std::vector<int>& peakCountOnLattices);

        void setSamplingPitch(SamplingPitch samplingPitch);
        void setSamplingPitch(float unitPitch, bool coverSecondaryMillerIndices);
        void setRefineWithExactLattice(bool flag);
        void setMaxPeaksToUseForIndexing(int maxPeaksToUseForIndexing);

        void setGradientDescentIterationsCount(GradientDescentIterationsCount gradientDescentIterationsCount);

      private:
        void precompute();
        void reducePeakCount(Eigen::Matrix3Xf& reciprocalPeaks_1_per_A);

        Eigen::Matrix3Xf precomputedSamplePoints;

        HillClimbingOptimizer hillClimbingOptimizer;
        SparsePeakFinder sparsePeakFinder;
        InverseSpaceTransform inverseSpaceTransform;

        float maxCloseToPointDeviation;
        int maxPeaksToUseForIndexing;

        HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_global;
        HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_additionalGlobal;
        HillClimbingOptimizer::hillClimbingAccuracyConstants_t hillClimbing_accuracyConstants_peaks;
        LatticeAssembler::accuracyConstants_t accuracyConstants_LatticeAssembler;
    };
} // namespace xgandalf
#endif /* INDEXERPLAIN_H_ */
