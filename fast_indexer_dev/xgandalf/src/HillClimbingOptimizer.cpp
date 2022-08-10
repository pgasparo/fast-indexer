/*
 * HillClimbingOptimizer.cpp
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

#include <xgandalf/HillClimbingOptimizer.h>
#include <cstddef>
#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

namespace xgandalf
{
    HillClimbingOptimizer::HillClimbingOptimizer()
        : transform()
        , hillClimbingAccuracyConstants()
    {
    }

    void HillClimbingOptimizer::performOptimization(const Matrix3Xf& pointsToTransform, Matrix3Xf& positionsToOptimize)
    {
        //    std::ofstream ofs("workfolder/tmp", std::ofstream::out);
        //    ofs << positionsToOptimize.transpose().eval() << endl;

        float& gamma = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.gamma;
        float& maxStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.maxStep;
        float& minStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.minStep;

        const float gamma_initial = gamma;
        const float maxStep_initial = maxStep;
        const float minStep_initial = minStep;

        const int initialIterationCount = hillClimbingAccuracyConstants.initialIterationCount;
        const int calmDownIterationCount = hillClimbingAccuracyConstants.calmDownIterationCount;
        const float calmDownFactor = hillClimbingAccuracyConstants.calmDownFactor;
        const int localFitIterationCount = hillClimbingAccuracyConstants.localFitIterationCount;
        const int localCalmDownIterationCount = hillClimbingAccuracyConstants.localCalmDownIterationCount;
        const float localCalmDownFactor = hillClimbingAccuracyConstants.localCalmDownFactor;

        transform.setPointsToTransform(pointsToTransform);
        const uint32_t maxPositionsPerIteration = 100; // TODO: find sweet spot. Maybe choose dependent on pointsToTransform.cols()
        Matrix3Xf positionsToOptimize_local;
        for (int64_t positionsProcessedCount = 0; positionsProcessedCount < positionsToOptimize.cols(); positionsProcessedCount += maxPositionsPerIteration)
        {
            uint32_t remainingPositionsCount = positionsToOptimize.cols() - positionsProcessedCount;
            uint32_t positionsCount_local = min(maxPositionsPerIteration, remainingPositionsCount);
            positionsToOptimize_local = positionsToOptimize.block(0, positionsProcessedCount, 3, positionsCount_local);

            previousStepDirection = Matrix3Xf::Zero(3, positionsToOptimize_local.cols());
            previousStepLength = Array<float, 1, Eigen::Dynamic>::Constant(1, positionsToOptimize_local.cols(), minStep + (maxStep - minStep) / 4);

            for (int i = 0; i < initialIterationCount; i++)
            {
                transform.clearLocalTransformFlag();
                transform.setRadialWeightingFlag();
                bool useStepOrthogonalization = true;

                performOptimizationStep(positionsToOptimize_local, useStepOrthogonalization);
                //        ofs << positionsToOptimize.transpose().eval() << endl;
            }

            for (int i = 0; i < calmDownIterationCount; i++)
            {
                transform.clearLocalTransformFlag();
                transform.setRadialWeightingFlag();
                bool useStepOrthogonalization = true;

                maxStep = maxStep * calmDownFactor;
                minStep = minStep * calmDownFactor;
                gamma = gamma * calmDownFactor;

                performOptimizationStep(positionsToOptimize_local, useStepOrthogonalization);
                //        ofs << positionsToOptimize.transpose().eval() << endl;
            }

            for (int i = 0; i < localFitIterationCount; i++)
            {
                transform.setLocalTransformFlag();
                transform.clearRadialWeightingFlag();
                bool useStepOrthogonalization = true;

                performOptimizationStep(positionsToOptimize_local, useStepOrthogonalization);
                //        ofs << positionsToOptimize.transpose().eval() << endl;
            }

            for (int i = 0; i < localCalmDownIterationCount; i++)
            {
                transform.setLocalTransformFlag();
                transform.clearRadialWeightingFlag();
                bool useStepOrthogonalization = false;

                maxStep = maxStep * localCalmDownFactor;
                minStep = minStep * localCalmDownFactor;
                gamma = gamma * localCalmDownFactor;

                performOptimizationStep(positionsToOptimize_local, useStepOrthogonalization);
                //        ofs << positionsToOptimize.transpose().eval() << endl;
            }

            positionsToOptimize.block(0, positionsProcessedCount, 3, positionsCount_local) = positionsToOptimize_local;

            gamma = gamma_initial;
            maxStep = maxStep_initial;
            minStep = minStep_initial;
        }

        // can be optimized! Does not always need to compute slope, closeToPoints and gradient
        transform.performTransform(positionsToOptimize);
        lastInverseTransformEvaluation = transform.getInverseTransformEvaluation();
    }

    void HillClimbingOptimizer::performOptimizationStep(Matrix3Xf& positionsToOptimize, bool useStepOrthogonalization)
    {
        transform.performTransform(positionsToOptimize);
        Matrix3Xf& gradient = transform.getGradient();
        RowVectorXf& closeToPointsCount = transform.getCloseToPointsCount();
        RowVectorXf& inverseTransformEvaluation = transform.getInverseTransformEvaluation();
        computeStep(gradient, closeToPointsCount, inverseTransformEvaluation, useStepOrthogonalization);
        //    cout << "step: " << endl << step << endl << endl;
        positionsToOptimize += step;
        //    cout << "positionsToOptimize: " << endl << positionsToOptimize << endl << endl;
    }

    void HillClimbingOptimizer::setPointsToTransformWeights(const RowVectorXf& pointsToTransformWeights)
    {
        transform.setPointsToTransformWeights(pointsToTransformWeights);
    }

    RowVectorXf& HillClimbingOptimizer::getLastInverseTransformEvaluation()
    {
        return lastInverseTransformEvaluation;
    }

    RowVectorXf& HillClimbingOptimizer::getCloseToPointsCount()
    {
        return transform.getCloseToPointsCount();
    }

    vector<vector<uint16_t>>& HillClimbingOptimizer::getPointsCloseToEvaluationPositions_indices()
    {
        return transform.getPointsCloseToEvaluationPositions_indices();
    }

    void HillClimbingOptimizer::computeStep(Matrix3Xf& gradient, RowVectorXf& closeToPointsCount, RowVectorXf& inverseTransformEvaluation,
                                            bool useStepOrthogonalization)
    {
        // reuse memory for processing for performance reasons
        Matrix3Xf& stepDirection = gradient;
        RowVectorXf& closeToPointsFactor = closeToPointsCount;
        RowVectorXf& functionEvaluationFactor = inverseTransformEvaluation;

        stepDirection.colwise().normalize();
        //    cout << "stepDirection " << endl << stepDirection << endl << endl<< previousStepDirection << endl << endl;

        Array<float, 1, Dynamic> directionChange = (stepDirection.array() * previousStepDirection.array()).matrix().colwise().sum();
        //    cout << "dirChange " << endl << directionChange << endl << endl;

        Array<float, 1, Dynamic> stepDirectionFactor = ((directionChange + 1) / 2).square().square() * 1.5 + 0.5; // directionChange in [-1 1]
        closeToPointsFactor = (closeToPointsCount.array() * (-1) + 0.8).square() * 6 + 0.5;                       // closeToPoints in [0 1]
        functionEvaluationFactor = ((inverseTransformEvaluation.array() * (-1) + 0.8) / 2).cube() * 4 + 0.3;      // functionEvaluation in [-1 1]
        //    cout << stepDirectionFactor << endl << endl << closeToPointsFactor << endl << endl << functionEvaluationFactor << endl;

        //    cout << "stepDirection " << endl << stepDirection << endl << endl;
        if (useStepOrthogonalization)
        {
            for (int i = 0; i < closeToPointsCount.size(); i++)
            {
                if (directionChange(i) < -0.4)
                {
                    stepDirection.col(i) = (stepDirection.col(i) + previousStepDirection.col(i)).normalized();
                    stepDirectionFactor(i) = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.directionChangeFactor;
                }
            }
        }
        //    cout << "stepDirection " << endl << stepDirection << endl << endl;

        previousStepDirection = stepDirection;

        const float minStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.minStep;
        const float maxStep = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.maxStep;
        const float gamma = hillClimbingAccuracyConstants.stepComputationAccuracyConstants.gamma;

        Array<float, 1, Eigen::Dynamic>& stepLength = previousStepLength; // reuse memory
        stepLength = ((0.5 * (minStep + (maxStep - minStep) * gamma) + 0.5 * previousStepLength) * stepDirectionFactor * closeToPointsFactor.array() *
                      functionEvaluationFactor.array())
                         .min(maxStep)
                         .max(minStep);

        step = stepDirection.array().rowwise() * stepLength;
    }

    void HillClimbingOptimizer::setStepComputationAccuracyConstants(stepComputationAccuracyConstants_t stepComputationAccuracyConstants)
    {
        hillClimbingAccuracyConstants.stepComputationAccuracyConstants = stepComputationAccuracyConstants;
    }

    void HillClimbingOptimizer::setHillClimbingAccuracyConstants(hillClimbingAccuracyConstants_t accuracyConstants)
    {
        hillClimbingAccuracyConstants = accuracyConstants;

        transform.setFunctionSelection(accuracyConstants.functionSelection);
        transform.setOptionalFunctionArgument(accuracyConstants.optionalFunctionArgument);
        transform.setMaxCloseToPointDeviation(accuracyConstants.maxCloseToPointDeviation);
    }
} // namespace xgandalf
