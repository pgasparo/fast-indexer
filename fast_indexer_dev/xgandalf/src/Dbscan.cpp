/*
 * Dbscan.cpp
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

#include <xgandalf/Dbscan.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace Eigen;

namespace xgandalf
{
    Dbscan::Dbscan()
        : squaredEpsilon(0)
        , maxEpsilon(0)
        , binWidth(0)
        , binWidth_reciprocal(0)
        , binsPerDimension(0)
        , binCount(0)
        , binCountMinus1(0)
    {
    }

    Dbscan::Dbscan(float maxEpsilon, float maxPossiblePointNorm)
    {
        init(maxEpsilon, maxPossiblePointNorm);
    }

    void Dbscan::init(float maxEpsilon, float maxPossiblePointNorm)
    {
        this->maxEpsilon = maxEpsilon;

        if (maxEpsilon < maxPossiblePointNorm / 50)
        {
            cout << "dbscan histogram would take too much memory! Reducing performance in trade of memory!" << endl << endl;
            this->maxEpsilon = maxPossiblePointNorm / 50;
        }

        binWidth = this->maxEpsilon;

        binWidth_reciprocal = 1 / binWidth;

        binsPerDimension = 2 * ceil(maxPossiblePointNorm / binWidth) + 2 + 1; //+2 for one extra border bin, where nothing should be inside.
        binCount = binsPerDimension * binsPerDimension * binsPerDimension;
        binCountMinus1 = binCount - 1;
        bin1Position.setConstant(-1.0f * binsPerDimension / 2 * binWidth);
        strides << 1, binsPerDimension, binsPerDimension * binsPerDimension;

        discretizationVolume.resize(binCount);

        const uint32_t typicalMaxPointsPerBin = 50;
        for (auto bin = discretizationVolume.begin(); bin != discretizationVolume.end(); ++bin)
        {
            bin->reserve(typicalMaxPointsPerBin);
        }

        neighbourBinIndexOffsets[0] = 0;
        int neighboursIndex = 1;
        for (int x = -1; x <= 1; x++)
        {
            for (int y = -1; y <= 1; y++)
            {
                for (int z = -1; z <= 1; z++)
                {
                    if (!(x == 0 && y == 0 && z == 0))
                    {
                        neighbourBinIndexOffsets[neighboursIndex] = x + y * strides.y() + z * strides.z();
                        neighboursIndex++;
                    }
                }
            }
        }

        mainNeighbourhood.reserve(27 * typicalMaxPointsPerBin);
        neighbourNeighbourhood.reserve(27 * typicalMaxPointsPerBin);
    }

    void Dbscan::computeClusters(vector<cluster_t>& clusters, const Matrix3Xf& points, uint16_t minPoints, float epsilon)
    {
        if (epsilon > maxEpsilon)
        {
            stringstream errStream;
            errStream << "epsilon must be smaller than maxEpsilon" << endl;
            errStream << "epsilon = " << epsilon << endl;
            errStream << "maxEpsilon = " << maxEpsilon << endl;
            throw WrongUsageException(errStream.str());
        }

        squaredEpsilon = epsilon * epsilon;

        clusters.reserve(50); // just for performance

        fillDiscretizationVolume(points);

        const auto end = usedBins.end();
        for (auto bin_p = usedBins.begin(); bin_p != end; ++bin_p)
        {
            bin_t& bin = **bin_p;
            for (uint32_t i = 0; i < bin.size(); ++i)
            {
                if (bin[i].visited)
                {
                    continue;
                }

                Neighbour currentPoint(*bin_p, &bin[i]);

                currentPoint.markVisited();

                uint32_t neighboursCount = regionQuery(mainNeighbourhood, currentPoint);
                if (neighboursCount >= minPoints)
                {
                    clusters.emplace_back();
                    clusters.back().reserve(5);
                    expandCluster(clusters.back(), mainNeighbourhood, currentPoint, minPoints);
                }
            }
        }

        cleanUpDiscretizationVolume();
    }

    void Dbscan::expandCluster(cluster_t& cluster, std::vector<Neighbour>& neighbourhood, Neighbour& currentPoint, uint16_t minPoints)
    {
        cluster.push_back(currentPoint.pointIndex());
        currentPoint.markIsMemberOfCluster();

        for (vector<Neighbour>::iterator neighbour_i = neighbourhood.begin(); neighbour_i != neighbourhood.end(); ++neighbour_i)
        { // cannot be made a for-each loop, since neighbourhood size changes
            Neighbour& neighbour = *neighbour_i;
            if (!neighbour.isVisited())
            {
                neighbour.markVisited();
                uint32_t neighboursCount = regionQuery(neighbourNeighbourhood, neighbour);
                if (neighboursCount >= minPoints)
                {
                    neighbourhood.insert(neighbourhood.end(), neighbourNeighbourhood.begin(), neighbourNeighbourhood.end());
                }
            }
            if (!neighbour.isMemberOfCluster())
            {
                cluster.push_back(neighbour.pointIndex());
                neighbour.markIsMemberOfCluster();
            }
        }
    }

    uint32_t Dbscan::regionQuery(std::vector<Neighbour>& nieghbourhood, const Neighbour& currentPoint)
    {
        uint32_t validNeighboursCount = 0;
        nieghbourhood.clear();

        const Vector4f& currentPointPos = currentPoint.point();
        for (int j = 0; j < 27; j++)
        {
            int neighbourBinIndexOffset = neighbourBinIndexOffsets[j];

            bin_t& neighbourBin = *(currentPoint.bin + neighbourBinIndexOffset);
            for (uint32_t i = 0; i < neighbourBin.size(); ++i)
            {
                const Vector4f& neighbourPos = neighbourBin[i].point;
                if ((currentPointPos - neighbourPos).squaredNorm() <= squaredEpsilon)
                {
                    validNeighboursCount++;
                    if (!neighbourBin[i].visitedAndMemberOfCluster)
                    {
                        nieghbourhood.emplace_back(&neighbourBin, &neighbourBin[i]);
                    }
                }
            }
        }

        return validNeighboursCount;
    }

    void Dbscan::fillDiscretizationVolume(const Eigen::Matrix3Xf& points)
    {
        uint32_t pointsCount = points.cols();

        for (uint32_t i = 0; i < pointsCount; ++i)
        {
            const Vector3f& point = points.col(i);
            const uint32_t index = getIndex(point);
            bin_t& bin = discretizationVolume[index];

            bin.emplace_back();
            auto& entry = bin.back();
            entry.point = Vector4f(point.x(), point.y(), point.z(), 0);
            entry.pointIndex = i;
            entry.visited = false;
            entry.isMemberOfCluster = false;
            entry.visitedAndMemberOfCluster = false;

            usedBins.insert(&bin);
        }
    }

    void Dbscan::cleanUpDiscretizationVolume()
    {
        const auto end = usedBins.end();
        for (auto bin_p = usedBins.begin(); bin_p != end; bin_p++)
        {
            (*bin_p)->clear();
        }
    }
} // namespace xgandalf
