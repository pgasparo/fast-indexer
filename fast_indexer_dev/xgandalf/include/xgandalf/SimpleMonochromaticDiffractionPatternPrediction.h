/*
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

#pragma once


#include "SimpleMonochromaticProjection.h"
#include <Eigen/Dense>

namespace xgandalf
{

    class SimpleMonochromaticDiffractionPatternPrediction
    {
      public:
        SimpleMonochromaticDiffractionPatternPrediction(const ExperimentSettings& experimentSettings);

        void getPeaksOnEwaldSphere(Eigen::Matrix3Xf& peaksOnEwaldSphere, Eigen::Matrix3Xi& millerIndices, const Lattice& lattice);
        void predictPattern(Eigen::Matrix2Xf& predictedPeaks, Eigen::Matrix3Xi& millerIndices, Eigen::Matrix3Xf& projectionDirections, const Lattice& lattice);

      private:
        SimpleMonochromaticProjection reciprocalToRealProjection;

        float maxResolutionAngle;
        float reflectionRadius;
        float reciprocalLambdaShort, reciprocalLambdaLong;
        float reciprocalLambdaShort_extended_squared, reciprocalLambdaLong_extended_squared;
        float detectorDistance;
    };

} // namespace xgandalf