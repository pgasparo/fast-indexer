/*
 * InverseSpaceTransform.cpp
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

#define _USE_MATH_DEFINES
#include <cmath>

#include <xgandalf/InverseSpaceTransform.h>
#include <assert.h>
#include <iostream>

using namespace Eigen;
using namespace std;

namespace xgandalf
{
    static inline void function1(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument);
    static inline void function1_periodic(ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument,
                                          Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& closeToPoint, float maxCloseToPointDeviation);
    static inline void function2(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
    static inline void function3(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
    static inline void function4(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
    static inline void function5(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
    static inline void function6(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
    static inline void function7(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument);
    static inline void function8(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope);
    static inline void function9(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument);

    InverseSpaceTransform::InverseSpaceTransform()
        : inverseTransformEvaluationScalingFactor(0)
        , resultsUpToDate(false)
    {
        accuracyConstants.functionSelection = 0;
        accuracyConstants.optionalFunctionArgument = 1;
        accuracyConstants.localTransform = false;
        accuracyConstants.radialWeighting = false;
        accuracyConstants.maxCloseToPointDeviation = 0.15;
    }

    InverseSpaceTransform::InverseSpaceTransform(float maxCloseToPointDeviation)
        : inverseTransformEvaluationScalingFactor(0)
        , resultsUpToDate(false)
    {
        accuracyConstants.functionSelection = 0;
        accuracyConstants.optionalFunctionArgument = 1;
        accuracyConstants.localTransform = false;
        accuracyConstants.radialWeighting = false;
        accuracyConstants.maxCloseToPointDeviation = maxCloseToPointDeviation;
    }

    void InverseSpaceTransform::performTransform(const Matrix3Xf& positionsToEvaluate)
    {
        float pointsToTransformCount_inverse = 1 / (float)pointsToTransform.cols();

        ArrayXXf x = pointsToTransform.transpose() * positionsToEvaluate;
        onePeriodicFunction(x);

        //    cout << slope << endl << endl << pointsToTransform << endl << endl << pointsToTransformWeights << endl << endl;
        Matrix3Xf fullGradient;
        if (accuracyConstants.localTransform)
        {
            fullGradient = (pointsToTransform.array().rowwise() * pointsToTransformWeights.array()).matrix() * slope.matrix();
            functionEvaluation = functionEvaluation * closeToPoint.matrix().cast<float>().array();
            slope = slope * closeToPoint.matrix().cast<float>().array();
        }
        //    cout << slope << endl << endl << functionEvaluation << endl << endl << fullGradient << endl << endl;

        gradient = (pointsToTransform.array().rowwise() * pointsToTransformWeights.array()).matrix() * slope.matrix();
        inverseTransformEvaluation = pointsToTransformWeights * functionEvaluation.matrix() * inverseTransformEvaluationScalingFactor;

        if (closeToPoint.rows() <= 255)
        {
            closeToPointsCount = closeToPoint.matrix().cast<uint8_t>().colwise().sum().cast<float>();
        }
        else
        {
            closeToPointsCount = closeToPoint.matrix().cast<uint16_t>().colwise().sum().cast<float>();
        }
        //    cout << gradient << endl << endl << inverseTransformEvaluation << endl << endl << closeToPoint << endl << endl;

        if (accuracyConstants.localTransform)
        {
            for (int i = 0; i < closeToPointsCount.size(); i++)
            {
                gradient.col(i) =
                    (closeToPointsCount(i) != 0) ? (gradient.col(i) * (1.0f / closeToPointsCount(i))) : (fullGradient.col(i) * pointsToTransformCount_inverse);
            }
        }
        else
        {
            gradient = gradient * pointsToTransformCount_inverse;
        }

        closeToPointsCount = closeToPointsCount * pointsToTransformCount_inverse;

        resultsUpToDate = true;
    }

    // cos(x * 2*pi).^optionalFunctionArgument
    static inline void function1(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument)
    {
        assert((optionalFunctionArgument - round(optionalFunctionArgument)) == 0);
        assert((int)optionalFunctionArgument % 2 != 0); // An even optionanFunctionArgument does not make sense with this function.

        if (optionalFunctionArgument == 1)
        {
            functionEvaluation = cos(x * (2 * M_PI));
            slope = -sin(x * (2 * M_PI));
        }
        else
        {
            functionEvaluation = pow(cos(x * (2 * M_PI)), (int)optionalFunctionArgument); // can be faster using manual pow for integer exponent!
            float n = optionalFunctionArgument;
            float scaling = pow(n, n / 2) / pow((n - 1), (int)(n - 1) / 2);
            slope = -scaling * sin(x * (2 * M_PI)) * pow(cos(x * (2 * M_PI)), (int)(n - 1));
        }
    }

    // cos(x * 2*pi).^optionalFunctionArgument
    static inline void function1_periodic(ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument,
                                          Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& closeToPoint, float maxCloseToPointDeviation)
    {
        assert((optionalFunctionArgument - round(optionalFunctionArgument)) == 0);
        assert((int)optionalFunctionArgument % 2 != 0); // An even optionanFunctionArgument does not make sense with this function.

        if (optionalFunctionArgument == 1)
        {
            functionEvaluation = cos(x * (2 * M_PI));
            slope = -sin(x * (2 * M_PI));

            float threshold = cos(maxCloseToPointDeviation * (2 * M_PI)); // can be precomputed
            closeToPoint = functionEvaluation > threshold;
        }
        else
        {
            functionEvaluation = pow(cos(x * (2 * M_PI)), (int)optionalFunctionArgument); // can be faster using manual pow for integer exponent!
            float n = optionalFunctionArgument;
            float scaling = pow(n, n / 2) / pow((n - 1), (int)(n - 1) / 2);
            slope = -scaling * sin(x * (2 * M_PI)) * pow(cos(x * (2 * M_PI)), (int)(n - 1));

            float threshold = pow(cos(maxCloseToPointDeviation * (2 * M_PI)), (int)optionalFunctionArgument);
            closeToPoint = functionEvaluation > threshold;
        }
    }

    // sawtooth
    static inline void function2(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
    {
        functionEvaluation = (-4 * abs(x) + 1);
        slope = -x.sign();
    }

    static inline void function3(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
    {
        auto x_p4 = x.square().square();

        functionEvaluation = (32 * x_p4 - 16 * x.square() + 1);
        slope = 128 * x.cube() - 32 * x;
    }

    // amplitude from 2, slope from 3
    static inline void function4(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
    {
        functionEvaluation = (-4 * abs(x) + 1);
        slope = 128 * x.cube() - 32 * x;
    }

    static inline void function5(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
    {
        functionEvaluation = -8 * x.square() + 1;
        slope = -2 * x;
    }

    static inline void function6(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
    {
        functionEvaluation = -32 * x.square().square() + 1;
        slope = -8 * x.cube();
    }

    static inline void function7(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument)
    {
        functionEvaluation = -(abs(x) - optionalFunctionArgument / 2).sign();
        slope = -x.sign() * (-functionEvaluation + 1) / 2;
    }

    static inline void function8(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope)
    {
        functionEvaluation = 8 * ((abs(x) - 0.5)).square() - 1;
        slope = 2 * (abs(x) - 0.5) * x.sign();
    }

    static inline void function9(const ArrayXXf& x, ArrayXXf& functionEvaluation, ArrayXXf& slope, float optionalFunctionArgument)
    {
        if (optionalFunctionArgument - round(optionalFunctionArgument) == 0)
        {
            //        functionEvaluation = pow(1 - 2*abs(x), (int) optionalFunctionArgument) * 2 - 1;
            //        slope = -x * pow(1 - 2*abs(x), (int) optionalFunctionArgument - 1) / (abs(x) + 0.0001);
            int exponent = (int)optionalFunctionArgument;

            ArrayXXf base = 1 - 2 * abs(x); // TODO: auto may be faster or slower... check!
            switch (exponent)
            { // just for performance, in case the compiler does not recognize the integer exponent
                case 1:
                    functionEvaluation = base * 2 - 1;
                    slope = -1 * x.sign();
                    break;
                case 2:
                    functionEvaluation = base.square() * 2 - 1;
                    slope = -x * base / (abs(x) + 0.0001);
                    break;
                case 3:
                    functionEvaluation = base.cube() * 2 - 1;
                    slope = -x * base.square() / (abs(x) + 0.0001);
                    break;
                case 4:
                    functionEvaluation = base.square().square() * 2 - 1;
                    slope = -x * base.cube() / (abs(x) + 0.0001);
                    break;
                case 5:
                    functionEvaluation = base.cube() * base.square() * 2 - 1;
                    slope = -x * base.square().square() / (abs(x) + 0.0001);
                    break;
                case 6:
                    functionEvaluation = base.cube().square() * 2 - 1;
                    slope = -x * base.cube() * base.square() / (abs(x) + 0.0001);
                    break;
                case 7:
                    functionEvaluation = base.cube() * base.square().square() * 2 - 1;
                    slope = -x * base.cube().square() / (abs(x) + 0.0001);
                    break;
                case 8:
                    functionEvaluation = base.square().square().square() * 2 - 1;
                    slope = -x * base.cube() * base.square().square() / (abs(x) + 0.0001);
                    break;
                case 9:
                    functionEvaluation = base.cube().cube() * 2 - 1;
                    slope = -x * base.square().square().square() / (abs(x) + 0.0001);
                    break;
                case 10:
                    functionEvaluation = base.cube().cube() * base * 2 - 1;
                    slope = -x * base.cube().cube() / (abs(x) + 0.0001);
                    break;
                case 11:
                    functionEvaluation = base.cube().cube() * base.square() * 2 - 1;
                    slope = -x * base.cube().cube() * base / (abs(x) + 0.0001);
                    break;
                case 12:
                    functionEvaluation = base.cube().square().square() * 2 - 1;
                    slope = -x * base.cube().cube() * base.square() / (abs(x) + 0.0001);
                    break;
                default:
                    functionEvaluation = pow(base, exponent) * 2 - 1;
                    slope = -x * pow(base, exponent - 1) / (abs(x) + 0.0001);
            }
        }
        else
        {
            functionEvaluation = pow(1 - 2 * abs(x), optionalFunctionArgument) * 2 - 1;
            slope = -x * pow(1 - 2 * abs(x), optionalFunctionArgument - 1) / (abs(x) + 0.0001);
        }
    }

    void InverseSpaceTransform::onePeriodicFunction(ArrayXXf& x)
    {
        switch (accuracyConstants.functionSelection)
        {
            case 1:
                function1_periodic(x, functionEvaluation, slope, accuracyConstants.optionalFunctionArgument, closeToPoint,
                                   accuracyConstants.maxCloseToPointDeviation);
                break;
            case 2:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function2(x, functionEvaluation, slope);
                break;
            case 3:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function3(x, functionEvaluation, slope);
                break;
            case 4:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function4(x, functionEvaluation, slope);
                break;
            case 5:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function5(x, functionEvaluation, slope);
                break;
            case 6:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function6(x, functionEvaluation, slope);
                break;
            case 7:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function7(x, functionEvaluation, slope, accuracyConstants.optionalFunctionArgument);
                break;
            case 8:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function8(x, functionEvaluation, slope);
                break;
            case 9:
                x = x - round(x);
                closeToPoint = abs(x) < accuracyConstants.maxCloseToPointDeviation;
                function9(x, functionEvaluation, slope, accuracyConstants.optionalFunctionArgument);
                break;
            default:
                stringstream errStream;
                errStream << "Selected function is not available.";
                throw BadInputException(errStream.str());
        }
    }

    void InverseSpaceTransform::setPointsToTransform(const Matrix3Xf& pointsToTransform)
    {
        resultsUpToDate = false;
        this->pointsToTransform = pointsToTransform;

        if (pointsToTransform.cols() != pointsToTransformWeights.cols())
        { // TODO: actually not a good choise...
            pointsToTransformWeights_userPreset = RowVectorXf::Ones(pointsToTransform.cols());
            update_pointsToTransformWeights();
        }
    }

    void InverseSpaceTransform::setPointsToTransformWeights(const RowVectorXf& pointsToTransformWeights)
    {
        resultsUpToDate = false;
        pointsToTransformWeights_userPreset = pointsToTransformWeights;
        update_pointsToTransformWeights();
    }

    void InverseSpaceTransform::update_pointsToTransformWeights()
    {
        if (accuracyConstants.radialWeighting)
        {
            Array<float, 1, Dynamic> radialWeight = pointsToTransform.colwise().squaredNorm().array().rsqrt();
            pointsToTransformWeights = (pointsToTransformWeights_userPreset.array() * radialWeight).matrix();
        }
        else
        {
            pointsToTransformWeights = pointsToTransformWeights_userPreset;
        }
        inverseTransformEvaluationScalingFactor = 1 / pointsToTransformWeights.sum();
    }

    void InverseSpaceTransform::setFunctionSelection(int functionSelection)
    {
        if (accuracyConstants.functionSelection != functionSelection)
        {
            resultsUpToDate = false;
            accuracyConstants.functionSelection = functionSelection;
        }
    }
    void InverseSpaceTransform::setOptionalFunctionArgument(float optionalFunctionArgument)
    {
        if (accuracyConstants.optionalFunctionArgument != optionalFunctionArgument)
        {
            resultsUpToDate = false;
            accuracyConstants.optionalFunctionArgument = optionalFunctionArgument;
        }
    }
    void InverseSpaceTransform::setLocalTransformFlag()
    {
        if (accuracyConstants.localTransform == false)
        {
            resultsUpToDate = false;
            accuracyConstants.localTransform = true;
        }
    }
    void InverseSpaceTransform::clearLocalTransformFlag()
    {
        if (accuracyConstants.localTransform == true)
        {
            resultsUpToDate = false;
            accuracyConstants.localTransform = false;
        }
    }
    void InverseSpaceTransform::setRadialWeightingFlag()
    {
        if (accuracyConstants.radialWeighting == false)
        {
            resultsUpToDate = false;
            accuracyConstants.radialWeighting = true;
        }

        update_pointsToTransformWeights();
    }
    void InverseSpaceTransform::clearRadialWeightingFlag()
    {
        if (accuracyConstants.radialWeighting == true)
        {
            resultsUpToDate = false;
            accuracyConstants.radialWeighting = false;
        }

        update_pointsToTransformWeights();
    }

    void InverseSpaceTransform::setMaxCloseToPointDeviation(float maxCloseToPointDeviation)
    {
        assert(maxCloseToPointDeviation < 0.5);
        if (accuracyConstants.maxCloseToPointDeviation != maxCloseToPointDeviation)
        {
            resultsUpToDate = false;
            accuracyConstants.maxCloseToPointDeviation = maxCloseToPointDeviation;
        }
    }

    Eigen::Matrix3Xf& InverseSpaceTransform::getGradient()
    {
        if (resultsUpToDate)
        {
            return gradient;
        }
        else
        {
            stringstream errStream;
            errStream << "Gradient not up to date, call performTransform() first.";
            throw BadInputException(errStream.str());
        }
    }
    Eigen::RowVectorXf& InverseSpaceTransform::getInverseTransformEvaluation()
    {
        if (resultsUpToDate)
        {
            return inverseTransformEvaluation;
        }
        else
        {
            stringstream errStream;
            errStream << "InverseTransformEvaluation not up to date, call performTransform() first.";
            throw BadInputException(errStream.str());
        }
    }
    Eigen::RowVectorXf& InverseSpaceTransform::getCloseToPointsCount()
    {
        if (resultsUpToDate)
        {
            return closeToPointsCount;
        }
        else
        {
            stringstream errStream;
            errStream << "CloseToPointsCount not up to date, call performTransform() first.";
            throw BadInputException(errStream.str());
        }
    }

    vector<vector<uint16_t>>& InverseSpaceTransform::getPointsCloseToEvaluationPositions_indices()
    {
        const int typicalMaxClosePointsCount = 100;
        pointsCloseToEvaluationPositions_indices.resize(closeToPoint.cols(), vector<uint16_t>(typicalMaxClosePointsCount));

        const auto end = pointsCloseToEvaluationPositions_indices.end();
        for (auto it = pointsCloseToEvaluationPositions_indices.begin(); it < end; ++it)
        {
            it->clear();
        }

        if (resultsUpToDate)
        {
            for (int pointIndex = 0; pointIndex < closeToPoint.rows(); pointIndex++)
            {
                for (int evaluationPositionIndex = 0; evaluationPositionIndex < closeToPoint.cols(); evaluationPositionIndex++)
                {
                    if (closeToPoint(pointIndex, evaluationPositionIndex))
                    {
                        pointsCloseToEvaluationPositions_indices[evaluationPositionIndex].push_back(pointIndex);
                    }
                }
            }

            return pointsCloseToEvaluationPositions_indices;
        }
        else
        {
            stringstream errStream;
            errStream << "closeToPoint not up to date, call performTransform() first.";
            throw BadInputException(errStream.str());
        }
    }
} // namespace xgandalf
