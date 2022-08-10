/*
 * experimentSettings.cpp
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

#include <xgandalf/ExperimentSettings.h>
#include <assert.h>
#include <sstream>

#include <cmath>

#define M_PI 3.14159265358979323846

using namespace std;
using namespace Eigen;

namespace xgandalf
{
    ExperimentSettings::ExperimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
                                           float pixelLength_m, float detectorRadius_pixel, float minRealLatticeVectorLength_A,
                                           float maxRealLatticeVectorLength_A, float reflectionRadius_1_per_A)
        : reflectionRadius_1_per_A(reflectionRadius_1_per_A)
        , latticeParametersKnown(false)
        , minRealLatticeVectorLength_A(minRealLatticeVectorLength_A)
        , maxRealLatticeVectorLength_A(maxRealLatticeVectorLength_A)
    {
        constructFromGeometryFileValues(coffset_m, clen_mm, beamEenergy_eV, divergenceAngle_deg, nonMonochromaticity, pixelLength_m, detectorRadius_pixel);

        minRealLatticeDeterminant_A3 = pow(minRealLatticeVectorLength_A, 3) * 0.1;
        maxRealLatticeDeterminant_A3 = pow(maxRealLatticeVectorLength_A, 3);

        minReciprocalLatticeVectorLength_1A = 1 / maxRealLatticeVectorLength_A;
        maxReciprocalLatticeVectorLength_1A = 1 / minRealLatticeVectorLength_A;
        minReciprocalLatticeDeterminant_1A3 = 1 / maxRealLatticeDeterminant_A3;
        maxReciprocalLatticeDeterminant_1A3 = 1 / minRealLatticeDeterminant_A3;

        differentRealLatticeVectorLengths_A.resize(2);
        differentRealLatticeVectorLengths_A[0] = minRealLatticeVectorLength_A;
        differentRealLatticeVectorLengths_A[1] = maxRealLatticeVectorLength_A;
    }

    ExperimentSettings::ExperimentSettings(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                           float nonMonochromaticity, float minRealLatticeVectorLength_A, float maxRealLatticeVectorLength_A,
                                           float reflectionRadius_1_per_A)
        : reflectionRadius_1_per_A(reflectionRadius_1_per_A)
        , latticeParametersKnown(false)
        , minRealLatticeVectorLength_A(minRealLatticeVectorLength_A)
        , maxRealLatticeVectorLength_A(maxRealLatticeVectorLength_A)
    {
        constructFromPrecomputedValues(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity);

        minRealLatticeDeterminant_A3 = pow(minRealLatticeVectorLength_A, 3);
        maxRealLatticeDeterminant_A3 = pow(maxRealLatticeVectorLength_A, 3);

        minReciprocalLatticeVectorLength_1A = 1 / maxRealLatticeVectorLength_A;
        maxReciprocalLatticeVectorLength_1A = 1 / minRealLatticeVectorLength_A;
        minReciprocalLatticeDeterminant_1A3 = 1 / maxRealLatticeDeterminant_A3;
        maxReciprocalLatticeDeterminant_1A3 = 1 / minRealLatticeDeterminant_A3;

        differentRealLatticeVectorLengths_A.resize(2);
        differentRealLatticeVectorLengths_A[0] = minRealLatticeVectorLength_A;
        differentRealLatticeVectorLengths_A[1] = maxRealLatticeVectorLength_A;
    }

    ExperimentSettings::ExperimentSettings(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg, float nonMonochromaticity,
                                           float pixelLength_m, float detectorRadius_pixel, const Lattice& sampleReciprocalLattice_1A, float tolerance,
                                           float reflectionRadius_1_per_A)
        : reflectionRadius_1_per_A(reflectionRadius_1_per_A)
        , latticeParametersKnown(true)
        , sampleReciprocalLattice_1A(sampleReciprocalLattice_1A)
        , latticeParametersTolerance(tolerance)
    {
        constructFromGeometryFileValues(coffset_m, clen_mm, beamEenergy_eV, divergenceAngle_deg, nonMonochromaticity, pixelLength_m, detectorRadius_pixel);
        deduceValuesFromSampleReciprocalLattice();
    }

    ExperimentSettings::ExperimentSettings(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                           float nonMonochromaticity, const Lattice& sampleReciprocalLattice_1A, float tolerance,
                                           float reflectionRadius_1_per_A)
        : reflectionRadius_1_per_A(reflectionRadius_1_per_A)
        , latticeParametersKnown(true)
        , sampleReciprocalLattice_1A(sampleReciprocalLattice_1A)
        , latticeParametersTolerance(tolerance)
    {
        constructFromPrecomputedValues(beamEenergy_eV, detectorDistance_m, detectorRadius_m, divergenceAngle_deg, nonMonochromaticity);
        deduceValuesFromSampleReciprocalLattice();
    }

    void ExperimentSettings::deduceValuesFromSampleReciprocalLattice()
    {
        Lattice sampleReciprocalLatticeMinimized_1A = sampleReciprocalLattice_1A;
        sampleReciprocalLatticeMinimized_1A.minimize();

        sampleRealLattice_A = sampleReciprocalLattice_1A.getReciprocalLattice();
        Lattice sampleRealLatticeMinimized_A = sampleRealLattice_A;
        sampleRealLatticeMinimized_A.minimize();

        realLatticeVectorLengths_A = sampleRealLatticeMinimized_A.getBasisVectorNorms();
        realLatticeVectorAngles_deg = sampleRealLatticeMinimized_A.getBasisVectorAngles_deg();
        realLatticeVectorAnglesNormalized_deg = sampleRealLatticeMinimized_A.getBasisVectorAnglesNormalized_deg();
        realLatticeDeterminant_A3 = abs(sampleRealLatticeMinimized_A.det());
        reciprocalLatticeVectorLengths_1A = sampleReciprocalLatticeMinimized_1A.getBasisVectorNorms();
        reciprocalLatticeVectorAngles_deg = sampleReciprocalLatticeMinimized_1A.getBasisVectorAngles_deg();
        reciprocalLatticeVectorAnglesNormalized_deg = sampleReciprocalLatticeMinimized_1A.getBasisVectorAnglesNormalized_deg();
        reciprocalLatticeDeterminant_1A3 = abs(sampleReciprocalLatticeMinimized_1A.det());

        // norms are ordered due to minimization of matrix
        assert(reciprocalLatticeVectorLengths_1A[0] <= reciprocalLatticeVectorLengths_1A[1] &&
               reciprocalLatticeVectorLengths_1A[1] <= reciprocalLatticeVectorLengths_1A[2]);
        assert(realLatticeVectorLengths_A[0] <= realLatticeVectorLengths_A[1] && realLatticeVectorLengths_A[1] <= realLatticeVectorLengths_A[2]);

        float minSimilarityFactor = 0.96;
        if (realLatticeVectorLengths_A[0] / realLatticeVectorLengths_A[1] > minSimilarityFactor)
        {
            if (realLatticeVectorLengths_A[1] / realLatticeVectorLengths_A[2] > minSimilarityFactor)
            {
                differentRealLatticeVectorLengths_A.resize(1);
                differentRealLatticeVectorLengths_A[0] = realLatticeVectorLengths_A.mean();
            }
            else
            {
                differentRealLatticeVectorLengths_A.resize(2);
                differentRealLatticeVectorLengths_A[0] = realLatticeVectorLengths_A.head(2).mean();
                differentRealLatticeVectorLengths_A[1] = realLatticeVectorLengths_A[2];
            }
        }
        else if (realLatticeVectorLengths_A[1] / realLatticeVectorLengths_A[2] > minSimilarityFactor)
        {
            differentRealLatticeVectorLengths_A.resize(2);
            differentRealLatticeVectorLengths_A[0] = realLatticeVectorLengths_A[0];
            differentRealLatticeVectorLengths_A[1] = realLatticeVectorLengths_A.tail(2).mean();
        }
        else
        {
            differentRealLatticeVectorLengths_A = realLatticeVectorLengths_A;
        }

        minRealLatticeVectorLength_A = realLatticeVectorLengths_A.minCoeff() * (1 - latticeParametersTolerance);
        maxRealLatticeVectorLength_A = realLatticeVectorLengths_A.maxCoeff() * (1 + latticeParametersTolerance);
        minRealLatticeDeterminant_A3 =
            realLatticeDeterminant_A3 * pow((1 - latticeParametersTolerance), 3); // TODO: check: possibly better (1 - latticeParametersTolerance)
        maxRealLatticeDeterminant_A3 = realLatticeDeterminant_A3 * pow((1 + latticeParametersTolerance), 3);
        minReciprocalLatticeVectorLength_1A = reciprocalLatticeVectorLengths_1A.minCoeff() * (1 - latticeParametersTolerance);
        maxReciprocalLatticeVectorLength_1A = reciprocalLatticeVectorLengths_1A.maxCoeff() * (1 + latticeParametersTolerance);
        minReciprocalLatticeDeterminant_1A3 = reciprocalLatticeDeterminant_1A3 * (1 - latticeParametersTolerance);
        maxReciprocalLatticeDeterminant_1A3 = reciprocalLatticeDeterminant_1A3 * (1 + latticeParametersTolerance);
    }

    void ExperimentSettings::constructFromGeometryFileValues(float coffset_m, float clen_mm, float beamEenergy_eV, float divergenceAngle_deg,
                                                             float nonMonochromaticity, float pixelLength_m, float detectorRadius_pixel)
    {
        this->nonMonochromaticity = nonMonochromaticity;

        detectorDistance_m = clen_mm * 1e-3 + coffset_m;
        detectorRadius_m = detectorRadius_pixel * pixelLength_m;
        divergenceAngle_rad = divergenceAngle_deg * M_PI / 180;
        maxResolutionAngle_rad = atan(detectorRadius_m / detectorDistance_m);

        float h_Plank = 4.135667662e-15; // Planck constant [eV*s]
        float c_light = 299792458;       // speed of light [m/s]
        lambda_A = h_Plank * c_light / beamEenergy_eV * 1e10;
        lambdaShort_A = lambda_A * (1 - nonMonochromaticity / 2);
        lambdaLong_A = lambda_A * (1 + nonMonochromaticity / 2);

        reciprocal_lambda_1A = 1 / lambda_A;
        reciprocal_lambdaShort_1A = 1 / lambdaShort_A;
        reciprocal_lambdaLong_1A = 1 / lambdaLong_A;
    }

    void ExperimentSettings::constructFromPrecomputedValues(float beamEenergy_eV, float detectorDistance_m, float detectorRadius_m, float divergenceAngle_deg,
                                                            float nonMonochromaticity)
    {
        this->detectorDistance_m = detectorDistance_m;
        this->detectorRadius_m = detectorRadius_m;
        this->nonMonochromaticity = nonMonochromaticity;

        divergenceAngle_rad = divergenceAngle_deg * M_PI / 180;
        maxResolutionAngle_rad = atan(detectorRadius_m / detectorDistance_m);

        float h_Plank = 4.135667662e-15; // Planck constant [eV*s]
        float c_light = 299792458;       // speed of light [m/s]
        lambda_A = h_Plank * c_light / beamEenergy_eV * 1e10;
        lambdaShort_A = lambda_A * (1 - nonMonochromaticity / 2);
        lambdaLong_A = lambda_A * (1 + nonMonochromaticity / 2);

        reciprocal_lambda_1A = 1 / lambda_A;
        reciprocal_lambdaShort_1A = 1 / lambdaShort_A;
        reciprocal_lambdaLong_1A = 1 / lambdaLong_A;
    }

    float ExperimentSettings::getDetectorDistance_m() const
    {
        return detectorDistance_m;
    }

    float ExperimentSettings::getDetectorRadius_m() const
    {
        return detectorRadius_m;
    }

    const Eigen::ArrayXf& ExperimentSettings::getDifferentRealLatticeVectorLengths_A() const
    {
        return differentRealLatticeVectorLengths_A;
    }

    float ExperimentSettings::getDivergenceAngle_rad() const
    {
        return divergenceAngle_rad;
    }

    float ExperimentSettings::getLambda_A() const
    {
        return lambda_A;
    }

    float ExperimentSettings::getLambdaLong_A() const
    {
        return lambdaLong_A;
    }

    float ExperimentSettings::getLambdaShort_A() const
    {
        return lambdaShort_A;
    }

    bool ExperimentSettings::isLatticeParametersKnown() const
    {
        return latticeParametersKnown;
    }

    float ExperimentSettings::getMaxRealLatticeDeterminant_A3() const
    {
        return maxRealLatticeDeterminant_A3;
    }

    float ExperimentSettings::getMaxRealLatticeVectorLength_A() const
    {
        return maxRealLatticeVectorLength_A;
    }

    float ExperimentSettings::getMaxReciprocalLatticeDeterminant_1A3() const
    {
        return maxReciprocalLatticeDeterminant_1A3;
    }

    float ExperimentSettings::getMaxReciprocalLatticeVectorLength_1A() const
    {
        return maxReciprocalLatticeVectorLength_1A;
    }

    float ExperimentSettings::getMaxResolutionAngle_rad() const
    {
        return maxResolutionAngle_rad;
    }

    float ExperimentSettings::getReflectionRadius() const
    {
        return reflectionRadius_1_per_A;
    }

    float ExperimentSettings::getMinRealLatticeDeterminant_A3() const
    {
        return minRealLatticeDeterminant_A3;
    }

    float ExperimentSettings::getMinRealLatticeVectorLength_A() const
    {
        return minRealLatticeVectorLength_A;
    }

    float ExperimentSettings::getMinReciprocalLatticeDeterminant_1A3() const
    {
        return minReciprocalLatticeDeterminant_1A3;
    }

    float ExperimentSettings::getMinReciprocalLatticeVectorLength_1A() const
    {
        return minReciprocalLatticeVectorLength_1A;
    }

    float ExperimentSettings::getNonMonochromaticity() const
    {
        return nonMonochromaticity;
    }

    float ExperimentSettings::getRealLatticeDeterminant_A3() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return realLatticeDeterminant_A3;
    }

    const Eigen::Array3f& ExperimentSettings::getRealLatticeVectorAngles_deg() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return realLatticeVectorAngles_deg;
    }

    const Eigen::Array3f& ExperimentSettings::getRealLatticeVectorAnglesNormalized_deg() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return realLatticeVectorAnglesNormalized_deg;
    }

    const Eigen::Array3f& ExperimentSettings::getRealLatticeVectorLengths_A() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return realLatticeVectorLengths_A;
    }

    float ExperimentSettings::getReciprocalLambda_1A() const
    {
        return reciprocal_lambda_1A;
    }

    float ExperimentSettings::getReciprocalLambdaLong_1A() const
    {
        return reciprocal_lambdaLong_1A;
    }

    float ExperimentSettings::getReciprocalLambdaShort_1A() const
    {
        return reciprocal_lambdaShort_1A;
    }

    float ExperimentSettings::getReciprocalLatticeDeterminant_1A3() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return reciprocalLatticeDeterminant_1A3;
    }

    const Eigen::Array3f& ExperimentSettings::getReciprocalLatticeVectorAngles_deg() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return reciprocalLatticeVectorAngles_deg;
    }

    const Eigen::Array3f& ExperimentSettings::getReciprocalLatticeVectorAnglesNormalized_deg() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return reciprocalLatticeVectorAnglesNormalized_deg;
    }

    const Eigen::Array3f& ExperimentSettings::getReciprocalLatticeVectorLengths_1A() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return reciprocalLatticeVectorLengths_1A;
    }

    const Lattice& ExperimentSettings::getSampleRealLattice_A() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return sampleRealLattice_A;
    }

    const Lattice& ExperimentSettings::getSampleReciprocalLattice_1A() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return sampleReciprocalLattice_1A;
    }

    float ExperimentSettings::getTolerance() const
    {
        if (!latticeParametersKnown)
        {
            stringstream errStream;
            errStream << "Function can only be used, if the lattice parameters are known" << endl;
            throw WrongUsageException(errStream.str());
        }

        return latticeParametersTolerance;
    }
} // namespace xgandalf
