/*
 * pointAutocorrelation.h
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

#ifndef POINTAUTOCORRELATION_H_
#define POINTAUTOCORRELATION_H_

#include <Eigen/Dense>
#include <limits>

namespace xgandalf
{

    // all autocorrelation results will have only half of the possible points, since symmetric points (at z < 0) will be removed

    void getPointAutocorrelation(Eigen::Matrix3Xf& autocorrelationPoints, const Eigen::Matrix3Xf& points, float minNormInAutocorrelation,
                                 float maxNormInAutocorrelation);

    void getPointAutocorrelation(Eigen::Matrix3Xf& autocorrelationPoints, Eigen::VectorXi& centerPointIndices, Eigen::VectorXi& shiftedPointIndices,
                                 const Eigen::Matrix3Xf& points, float minNormInAutocorrelation, float maxNormInAutocorrelation);

} // namespace xgandalf
#endif /* POINTAUTOCORRELATION_H_ */
