/*
 * DetectorToReciprocalSpaceTransform.h
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

#ifndef DETECTORTORECIPROCALSPACETRANSFORM_H_
#define DETECTORTORECIPROCALSPACETRANSFORM_H_

#include <Eigen/Dense>
#include <xgandalf/ExperimentSettings.h>

namespace xgandalf
{
    class DetectorToReciprocalSpaceTransform
    {
      public:
        DetectorToReciprocalSpaceTransform(const ExperimentSettings& experimentSettings);

        // coordinate system same as reciprocal x-z
        void computeReciprocalPeaksFromDetectorPeaks(Eigen::Matrix3Xf& reciprocalPeaks_A, const Eigen::Matrix2Xf& detectorPeaks_m);

      private:
        float reciprocal_lambda_1A;
        float detectorDistance_m;
    };
} // namespace xgandalf
#endif /* DETECTORTORECIPROCALSPACETRANSFORM_H_ */
