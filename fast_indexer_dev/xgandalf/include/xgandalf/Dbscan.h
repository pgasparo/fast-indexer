/*
 * Dbscan.h
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

#ifndef DBSCAN_H_
#define DBSCAN_H_

#include "WrongUsageException.h"
#include "eigenSTLContainers.h"
#include <Eigen/Dense>
#include <ctype.h>
#include <set>
#include <vector>

namespace xgandalf
{
    class Dbscan
    {
      public:
        typedef std::vector<uint32_t> cluster_t;

        // maxEpsilon defines performance. Best performance for maxEpsilon == epsilon. Too small maxEpsilon means very big discretizationVolume
        Dbscan(float maxEpsilon, float maxPossiblePointNorm);
        Dbscan();

        void init(float maxEpsilon, float maxPossiblePointNorm);

        void computeClusters(std::vector<cluster_t>& clusters, const Eigen::Matrix3Xf& points, uint16_t minPoints, float epsilon);

      private:
        typedef struct
        {
            Eigen::Vector4f point;
            uint32_t pointIndex;

            bool visited;
            bool isMemberOfCluster;
            bool visitedAndMemberOfCluster;

          public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        } binEntry_t;

        typedef std::vector<binEntry_t> bin_t;

        class Neighbour
        {
          public:
            Neighbour(bin_t* bin, binEntry_t* binEntry)
                : bin(bin)
                , binEntry(binEntry)
            {
            }

            bin_t* bin;
            binEntry_t* binEntry;

            inline bool isVisited() const
            {
                return (*binEntry).visited;
            }
            inline bool isMemberOfCluster() const
            {
                return (*binEntry).isMemberOfCluster;
            }
            inline uint32_t pointIndex() const
            {
                return (*binEntry).pointIndex;
            }
            inline Eigen::Vector4f point() const
            {
                return (*binEntry).point;
            }
            inline void markVisited()
            {
                (*binEntry).visited = true;
                if ((*binEntry).isMemberOfCluster)
                {
                    (*binEntry).visitedAndMemberOfCluster = true;
                }
            }
            inline void markIsMemberOfCluster()
            {
                (*binEntry).isMemberOfCluster = true;
                if ((*binEntry).visited)
                {
                    (*binEntry).visitedAndMemberOfCluster = true;
                }
            }
        };

        void expandCluster(cluster_t& cluster, std::vector<Neighbour>& neighbourhood, Neighbour& currentPoint, uint16_t minPoints);
        uint32_t regionQuery(std::vector<Neighbour>& nieghbourhood,
                             const Neighbour& currentPoint); // returns neighbours count. Visited neighbours that are members
                                                             // of a cluster are not included in neighbourhood, but counted

        float squaredEpsilon;

        // for less reallocation
        std::vector<Neighbour> mainNeighbourhood;
        std::vector<Neighbour> neighbourNeighbourhood;

        // stuff for dicretizationVolume
        float maxEpsilon;

        std::vector<bin_t> discretizationVolume;
        std::set<bin_t*> usedBins;

        int32_t neighbourBinIndexOffsets[27];

        float binWidth, binWidth_reciprocal;
        int binsPerDimension;
        uint32_t binCount, binCountMinus1;
        Eigen::Vector3f bin1Position;
        Eigen::Vector3d strides;

        void fillDiscretizationVolume(const Eigen::Matrix3Xf& points);
        void cleanUpDiscretizationVolume();

        // fast, but insecure: if point lies out of scope, it gets relocated somewhere inside the scope.
        inline uint32_t getIndex(const Eigen::Vector3f& position) const
        {
            Eigen::Vector3d temp = ((position - bin1Position) * binWidth_reciprocal).cast<double>().array().floor().matrix();
            return std::min((uint32_t)(temp[0] + temp.tail(2).dot(strides.tail(2))), binCountMinus1); // min could be avoided for the price of unsafety
        }

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
} // namespace xgandalf
#endif /* DBSCAN_H_ */
