/*
 * InverseSpaceTransform.h
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

#ifndef INVERSESPACETRANSFORM_H_
#define INVERSESPACETRANSFORM_H_

#include "BadInputException.h"
#include <Eigen/Dense>
#include <ctype.h>
#include <vector>

namespace xgandalf
{

    class InverseSpaceTransform
    {
      public:
        typedef struct
        {
            int functionSelection;
            float optionalFunctionArgument;
            bool localTransform;
            bool radialWeighting;

            float maxCloseToPointDeviation;
        } accuracyConstants_t;

        InverseSpaceTransform();
        InverseSpaceTransform(float maxCloseToPointDeviation);

        void performTransform(const Eigen::Matrix3Xf& positionsToEvaluate);

        void setPointsToTransform(const Eigen::Matrix3Xf& pointsToTransform);
        void setPointsToTransformWeights(const Eigen::RowVectorXf& pointsToTransformWeights);

        void setMaxCloseToPointDeviation(float maxCloseToPointDeviation);
        void setFunctionSelection(int functionSelection);
        void setOptionalFunctionArgument(float optionalFunctionArgument);
        void setLocalTransformFlag();
        void clearLocalTransformFlag();
        void setRadialWeightingFlag();
        void clearRadialWeightingFlag();

        Eigen::Matrix3Xf& getGradient();
        Eigen::RowVectorXf& getInverseTransformEvaluation();
        Eigen::RowVectorXf& getCloseToPointsCount();

        std::vector<std::vector<uint16_t>>& getPointsCloseToEvaluationPositions_indices();

      private:
        void onePeriodicFunction(Eigen::ArrayXXf& x);

        void update_pointsToTransformWeights();

        Eigen::Matrix3Xf pointsToTransform;
        Eigen::RowVectorXf pointsToTransformWeights_userPreset;
        Eigen::RowVectorXf pointsToTransformWeights;

        accuracyConstants_t accuracyConstants;

        // output
        Eigen::Matrix3Xf gradient;
        Eigen::RowVectorXf inverseTransformEvaluation;
        Eigen::RowVectorXf closeToPointsCount;

        // interna
        Eigen::ArrayXXf functionEvaluation;
        Eigen::ArrayXXf slope;
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> closeToPoint;

        std::vector<std::vector<uint16_t>> pointsCloseToEvaluationPositions_indices;

        float inverseTransformEvaluationScalingFactor;

        bool resultsUpToDate;
    };
} // namespace xgandalf
#endif /* INVERSESPACETRANSFORM_H_ */
