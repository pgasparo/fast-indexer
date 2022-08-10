/*
 * samplePointsFiltering.h
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

#ifndef SAMPLEPOINTSFILTERING_H_
#define SAMPLEPOINTSFILTERING_H_

#include "BadInputException.h"
#include <Eigen/Dense>

namespace xgandalf
{
    void filterSamplePointsForNorm(Eigen::Matrix3Xf& samplePoints, const Eigen::ArrayXf& allowedNorms, float allowedTolerance);
}

#endif /* SAMPLEPOINTSFILTERING_H_ */
