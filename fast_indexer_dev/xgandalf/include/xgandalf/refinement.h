/*
 * refinement.h
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


#ifndef REFINEMENT_H_
#define REFINEMENT_H_

#include <Eigen/Dense>

namespace xgandalf
{

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void getGradient_reciprocalPeakMatch_meanDist(Eigen::Matrix3f& gradient, const Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N);

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void getGradient_detectorAngleMatch(Eigen::Matrix3f& gradient, const Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N);

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void getGradient_reciprocalPeakMatch_meanSquaredDist(Eigen::Matrix3f& gradient, const Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M,
                                                         const Eigen::Matrix3Xf& N);

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanSquaredDist(Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N);

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanDist_peaksAndAngle(Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N);

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    // TODO: Function can be made faster easily... nevertheless double is really necessary...
    void getGradient_detectorAngleMatchFixedLattice(Eigen::Vector3f& rotationAnglesGradient_detectorAngle,
                                                    Eigen::Vector3f& rotationAnglesGradient_distFromPoints, const Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M,
                                                    const Eigen::Matrix3Xf& N);
    // B: reciprocal basis
    // B_sample: sample reciprocal basis (with correct lattice parameters)
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanSquaredDist_fixedBasisParameters_kabsch(Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N,
                                                                           const Eigen::Matrix3f& B_sample);

    // B: reciprocal basis
    // B_sample: sample reciprocal basis (with correct lattice parameters)
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanSquaredDist_fixedBasisParameters(Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N,
                                                                    const Eigen::Matrix3f& B_sample);

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanDist_detectorAngleMatchFixedParameters(Eigen::Matrix3f& B, const Eigen::Matrix3Xf& M, const Eigen::Matrix3Xf& N);

} // namespace xgandalf

#endif /* REFINEMENT_H_ */