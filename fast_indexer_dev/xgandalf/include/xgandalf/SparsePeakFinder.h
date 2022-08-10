/*
 * SparsePeakFinder.h
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

#ifndef SPARSEPEAKFINDER_H_
#define SPARSEPEAKFINDER_H_

#include "WrongUsageException.h"
#include <Eigen/Dense>
#include <algorithm>
#include <ctype.h>
#include <vector>

namespace xgandalf
{

    class SparsePeakFinder
    {
      public:
        SparsePeakFinder();
        SparsePeakFinder(float minDistanceBetweenRealPeaks, float maxPossiblePointNorm);
        void precompute(float minDistanceBetweenRealPeaks, float maxPossiblePointNorm);

        // the only warranty this function makes is, that peaks that are separated by more than minDistanceBetweenRealPeaks are found. Additionally some peaks
        // might be found that are not real peaks
        void findPeaks_fast(Eigen::Matrix3Xf& pointPositions, Eigen::RowVectorXf& pointValues);

      private:
        // fast, but insecure: if point lies out of scope, it gets relocated somewhere inside the scope. Maximum bins per direction is 200;
        inline uint32_t getIndex(const Eigen::Vector3f& position) const
        {
            return std::min((uint32_t)((position - bin1Position) * binWidth_reciprocal).array().floor().matrix().dot(strides),
                            binCountMinus1); // min could be avoided for the price of unsafety
        }

        typedef struct
        {
            int pointIndex;
            float value;
        } bin_t;

        float minDistanceBetweenRealPeaks;

        float minDistanceBetweenRealPeaks_squared;

        float binWidth, binWidth_reciprocal;
        int binsPerDimension;
        uint32_t binCountMinus1;
        Eigen::Vector3f bin1Position;
        Eigen::Vector3f strides;

        std::vector<bin_t> discretizationVolume;
        int32_t neighbourBinIndexOffsets[26];

        bool precomputed;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };


} // namespace xgandalf

#endif /* SPARSEPEAKFINDER_H_ */
