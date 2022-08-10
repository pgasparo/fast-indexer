/*
 * ReciprocalSpaceTranslation.cpp
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

#include <xgandalf/DetectorToReciprocalSpaceTransform.h>

using namespace Eigen;
using namespace std;

namespace xgandalf
{
    DetectorToReciprocalSpaceTransform::DetectorToReciprocalSpaceTransform(const ExperimentSettings& experimentSettings)
    {
        reciprocal_lambda_1A = experimentSettings.getReciprocalLambda_1A();
        detectorDistance_m = experimentSettings.getDetectorDistance_m();
    }

    void DetectorToReciprocalSpaceTransform::computeReciprocalPeaksFromDetectorPeaks(Matrix3Xf& reciprocalPeaks_A, const Matrix2Xf& detectorPeaks_m)
    {
        Matrix3Xf backprojectionDirectionVectors(3, detectorPeaks_m.cols());
        backprojectionDirectionVectors.row(0).setConstant(detectorDistance_m);
        backprojectionDirectionVectors.row(1) = -1 * detectorPeaks_m.row(0); // detector x-coordinate is -y coordinate in reciprocal space
        backprojectionDirectionVectors.row(2) = detectorPeaks_m.row(1);

        reciprocalPeaks_A = backprojectionDirectionVectors.colwise().normalized() * reciprocal_lambda_1A;
        reciprocalPeaks_A.row(0) -= RowVectorXf::Constant(reciprocalPeaks_A.cols(), reciprocal_lambda_1A);
    }
} // namespace xgandalf
