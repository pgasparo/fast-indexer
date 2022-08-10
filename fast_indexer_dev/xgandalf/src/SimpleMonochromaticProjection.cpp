/*
 * SimpleMonochromaticProjection.cpp
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

#include <xgandalf/SimpleMonochromaticProjection.h>

namespace xgandalf
{
    using namespace std;
    using namespace Eigen;

    SimpleMonochromaticProjection::SimpleMonochromaticProjection(const ExperimentSettings& experimentSettings)
        : ReciprocalToRealProjection(experimentSettings)
    {
        reciprocalLambda_1A = experimentSettings.getReciprocalLambda_1A();
    }


    void SimpleMonochromaticProjection::project(Matrix2Xf& projectedPeaks, const Matrix3Xf& reciprocalPeaks)
    {
        RowVectorXf yzSquaredNorms = reciprocalPeaks.bottomRows(2).colwise().squaredNorm();
        // RowVectorXf rayOriginsX = (reciprocalPeaks.row(0) + (yzSquaredNorms.array() / reciprocalPeaks.row(0).array()).matrix()) / 2;
        RowVectorXf rayOriginsX;
        rayOriginsX.setConstant(yzSquaredNorms.size(), -1 * reciprocalLambda_1A);

        projectedPeaks =
            reciprocalPeaks.bottomRows(2).array().rowwise() / (reciprocalPeaks.row(0) - rayOriginsX).array() * experimentSettings.getDetectorDistance_m();
    }
} // namespace xgandalf
