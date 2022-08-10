/*
 * samplePointsGenerator.cpp
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
#include <xgandalf/SamplePointsGenerator.h>
#include <xgandalf/eigenDiskImport.h>
#include <xgandalf/eigenSTLContainers.h>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <vector>
#include <fstream>
#include <iostream>

namespace xgandalf
{
    using namespace std;
    using namespace Eigen;

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#ifndef PRECOMPUTED_DATA_DIR
#define PRECOMPUTED_DATA_DIR "define the path to the precomputed sample points in the cmake script"
#endif

    SamplePointsGenerator::SamplePointsGenerator()
    {
        precomputedSamplePointsPath = STRINGIFY(PRECOMPUTED_DATA_DIR);
    }

    inline static float getClosestArrayElement(ArrayXf arr, float value)
    {
        int minIndex;
        (arr - value).abs().minCoeff(&minIndex);
        return arr[minIndex];
    }

    // clang-format off
void SamplePointsGenerator::loadPrecomputedSamplePoints(Matrix3Xf& samplePoints, float unitPitch, float tolerance)
{
    Array< float, 1, Eigen::Dynamic > pitches, tolerances;

    stringstream fullPath;
    fullPath << precomputedSamplePointsPath << "/pitches";
    loadEigenMatrixFromDisk(pitches, fullPath.str());

    fullPath.str(string());
    fullPath << precomputedSamplePointsPath << "/tolerances";
    loadEigenMatrixFromDisk(tolerances, fullPath.str());

    fullPath.str(string());
    fullPath << precomputedSamplePointsPath <<
            "/pitch" << getClosestArrayElement(pitches, unitPitch) <<
            "_tolerance" << getClosestArrayElement(tolerances, tolerance);
    MatrixX3f samplePoints_T;
    loadEigenMatrixFromDisk(samplePoints_T, fullPath.str());
    samplePoints = samplePoints_T.transpose();  //workaround, because sample points are stored in a better human readable way
}
    // clang-format on

    void SamplePointsGenerator::getDenseGrid(Matrix3Xf& samplePoints, float unitPitch, float minRadius, float maxRadius)
    {
        EigenSTL::vector_Vector3f tmpSamplePoints;

        VectorXf xSamples, ySamples, zSamples;

        int samplesPerRadius = round(1 / unitPitch);
        float minRadiusSquared = minRadius * minRadius;
        float maxRadiusSquared = maxRadius * maxRadius;

        xSamples.setLinSpaced(samplesPerRadius * 2, -maxRadius, maxRadius);
        ySamples.setLinSpaced(samplesPerRadius * 2, -maxRadius, maxRadius);
        zSamples.setLinSpaced(samplesPerRadius, 0, maxRadius);

        for (int xIndex = 0; xIndex < xSamples.size(); xIndex++)
        {
            for (int yIndex = 0; yIndex < ySamples.size(); yIndex++)
            {
                for (int zIndex = 0; zIndex < zSamples.size(); zIndex++)
                {
                    Vector3f samplePoint(xSamples[xIndex], ySamples[yIndex], zSamples[zIndex]);
                    float squaredNorm = samplePoint.squaredNorm();
                    if (squaredNorm >= minRadiusSquared && squaredNorm <= maxRadiusSquared)
                    {
                        tmpSamplePoints.push_back(samplePoint);
                    }
                }
            }
        }

        samplePoints = Map<Matrix3Xf>(tmpSamplePoints[0].data(), 3, tmpSamplePoints.size()); // copy
    }

    void SamplePointsGenerator::getTightGrid(Matrix3Xf& samplePoints, float unitPitch, float tolerance, const VectorXf radii)
    {
        samplePoints.resize(3, 0);

        for (int i = 0; i < radii.size(); i++)
        {
            float radius = radii[i];
            float adaptedPitch = unitPitch / radii[i] * radii.maxCoeff();

            Matrix3Xf newSamplingPoints;
            loadPrecomputedSamplePoints(newSamplingPoints, adaptedPitch, tolerance);
            newSamplingPoints *= radius;

            Matrix3Xf tmp(3, samplePoints.cols() + newSamplingPoints.cols());
            tmp << samplePoints, newSamplingPoints;
            samplePoints = tmp;
        }
    }
} // namespace xgandalf
