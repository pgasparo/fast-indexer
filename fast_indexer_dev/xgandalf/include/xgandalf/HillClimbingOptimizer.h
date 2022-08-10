/*
 * HillClimbingOptimizer.h
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

#ifndef HILLCLIMBINGOPTIMIZER_H_
#define HILLCLIMBINGOPTIMIZER_H_

#include <Eigen/Dense>
#include <xgandalf/InverseSpaceTransform.h>

namespace xgandalf
{
    class HillClimbingOptimizer
    {
      public:
        typedef struct
        {
            float gamma;
            float maxStep;
            float minStep;
            float directionChangeFactor;
        } stepComputationAccuracyConstants_t;

        typedef struct
        {
            int initialIterationCount;
            int calmDownIterationCount;
            float calmDownFactor;
            int localFitIterationCount;
            int localCalmDownIterationCount;
            float localCalmDownFactor;

            int functionSelection;
            float optionalFunctionArgument;
            float maxCloseToPointDeviation;

            stepComputationAccuracyConstants_t stepComputationAccuracyConstants;
        } hillClimbingAccuracyConstants_t;

        HillClimbingOptimizer();

        void performOptimization(const Eigen::Matrix3Xf& pointsToTransform, Eigen::Matrix3Xf& positionsToOptimize);
        Eigen::RowVectorXf& getLastInverseTransformEvaluation();
        Eigen::RowVectorXf& getCloseToPointsCount();
        std::vector<std::vector<uint16_t>>& getPointsCloseToEvaluationPositions_indices();

        void setHillClimbingAccuracyConstants(hillClimbingAccuracyConstants_t accuracyConstants);

        // optional
        void setPointsToTransformWeights(const Eigen::RowVectorXf& pointsToTransformWeights);

      public:
        void setStepComputationAccuracyConstants(stepComputationAccuracyConstants_t stepComputationAccuracyConstants);

        // watch out! gradient, closeToPointsCount and inverseTransformEvaluation are changed in this function (for performance reasons)!
        void computeStep(Eigen::Matrix3Xf& gradient, Eigen::RowVectorXf& closeToPointsCount, Eigen::RowVectorXf& inverseTransformEvaluation,
                         bool useStepOrthogonalization);

        void performOptimizationStep(Eigen::Matrix3Xf& positionsToOptimize, bool useStepOrthogonalization);

        InverseSpaceTransform transform;
        hillClimbingAccuracyConstants_t hillClimbingAccuracyConstants;

        // interna
        Eigen::Matrix3Xf step;

        Eigen::Matrix3Xf previousStepDirection;
        Eigen::Array<float, 1, Eigen::Dynamic> previousStepLength;

        Eigen::RowVectorXf lastInverseTransformEvaluation;
    };
} // namespace xgandalf

#endif /* HILLCLIMBINGOPTIMIZER_H_ */
