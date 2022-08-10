/*
 * SparsePeakFinder.cpp
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

#include <xgandalf/BadInputException.h>
#include <xgandalf/SparsePeakFinder.h>

#include <math.h>
#include <sstream>

using namespace std;
using namespace Eigen;

namespace xgandalf
{
    SparsePeakFinder::SparsePeakFinder()
        : minDistanceBetweenRealPeaks(0)
        , minDistanceBetweenRealPeaks_squared(0)
        , binWidth(0)
        , binWidth_reciprocal(0)
        , binsPerDimension(0)
        , binCountMinus1(0)
    {
        precomputed = false;
    }

    SparsePeakFinder::SparsePeakFinder(float minDistanceBetweenRealPeaks, float maxPossiblePointNorm)
    {
        precompute(minDistanceBetweenRealPeaks, maxPossiblePointNorm);
    }

    void SparsePeakFinder::precompute(float minDistanceBetweenRealPeaks, float maxPossiblePointNorm)
    {
        this->minDistanceBetweenRealPeaks = minDistanceBetweenRealPeaks;

        this->minDistanceBetweenRealPeaks_squared = minDistanceBetweenRealPeaks * minDistanceBetweenRealPeaks;

        binWidth = sqrt(minDistanceBetweenRealPeaks * minDistanceBetweenRealPeaks / 3); // minDistanceBetweenRealPeaks is diagonal of the cube
        binWidth_reciprocal = 1 / binWidth;

        binsPerDimension = 2 * ceil(maxPossiblePointNorm / binWidth) + 2 + 1; //+2 for one extra border bin, where nothing should be inside.
        binCountMinus1 = binsPerDimension * binsPerDimension * binsPerDimension - 1;
        bin1Position.setConstant(-1.0f * binsPerDimension / 2 * binWidth);
        strides << 1, binsPerDimension, binsPerDimension * binsPerDimension;

        discretizationVolume.resize(binCountMinus1 + 1);

        int neighboursIndex = 0;
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

        if (binsPerDimension > 200)
        { // pow(2^23,1/3) == 200  (float 23 bit matissa)
            stringstream errStream;
            errStream
                << "Created discretizationVolume is too big. Function getIndex() will not work (because of the late casting of float to uint32 in getIndex).";
            throw BadInputException(errStream.str());
        }

        precomputed = true;
    }

    void SparsePeakFinder::findPeaks_fast(Matrix3Xf& pointPositions, RowVectorXf& pointValues)
    {
        if (!precomputed)
        {
            stringstream errStream;
            errStream << "The sparsePeakFinder is used without precomputation!" << endl;
            throw WrongUsageException(errStream.str());
        }

        fill(discretizationVolume.begin(), discretizationVolume.end(), bin_t({-1, numeric_limits<float>::lowest()})); // possibly with 0 faster

        for (int pointIndex = 0; pointIndex < pointPositions.cols(); pointIndex++)
        {
            uint32_t index = getIndex(pointPositions.col(pointIndex));
            bin_t& currentBin = discretizationVolume[index];
            if (pointValues[pointIndex] > currentBin.value)
            {
                currentBin.value = pointValues[pointIndex];
                currentBin.pointIndex = pointIndex;
            }
        }

        int peakCount = 0;
        Matrix3Xf peakPositions(3, pointPositions.cols());
        RowVectorXf peakValues(pointValues.size());

        // For not checking whether index is out of borders, an extra border bin is added. Only inner borders are checked for peaks
        uint32_t binsPerDimensionMinus1 = binsPerDimension - 1;
        for (uint32_t y = 1; y < binsPerDimensionMinus1; y++)
        {
            for (uint32_t z = 1; z < binsPerDimensionMinus1; z++)
            {
                int lineStartIndex = y * strides.y() + z * strides.z();
                for (uint32_t x = 1; x < binsPerDimensionMinus1; x++)
                {
                    int binIndex = lineStartIndex + x;
                    if (discretizationVolume[binIndex].value > numeric_limits<float>::lowest()) // something inside bin
                    {
                        bin_t* currentBin_p = &discretizationVolume[binIndex];
                        bin_t& currentBin = *currentBin_p;
                        bool isPeak = true;
                        for (int i = 0; i < 26; ++i)
                        {
                            bin_t& neighbourBin = *(currentBin_p + neighbourBinIndexOffsets[i]);
                            if (neighbourBin.value > currentBin.value &&
                                (pointPositions.col(neighbourBin.pointIndex) - pointPositions.col(currentBin.pointIndex)).squaredNorm() <
                                    minDistanceBetweenRealPeaks_squared)
                            {
                                isPeak = false;
                                break;
                            }
                        }

                        if (isPeak)
                        {
                            peakPositions.col(peakCount) = pointPositions.col(currentBin.pointIndex); // auf überschriebenes wird zugegriffen!!!
                            peakValues[peakCount] = pointValues[currentBin.pointIndex];
                            peakCount++;
                        }
                    }
                }
            }
        }

        peakPositions.conservativeResize(NoChange, peakCount);
        peakValues.conservativeResize(peakCount);

        pointPositions.swap(peakPositions);
        pointValues.swap(peakValues);
    }
} // namespace xgandalf
