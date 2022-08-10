/*
 * eigenSTLContainers.h
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

#ifndef EIGENSTLCONTAINERS_H_
#define EIGENSTLCONTAINERS_H_

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <vector>

namespace xgandalf
{
    namespace EigenSTL
    {
        typedef std::vector<Eigen::Vector2f, Eigen::aligned_allocator<Eigen::Vector2f>> vector_Vector2f;
        typedef std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> vector_Vector2d;
        typedef std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> vector_Vector3f;
        typedef std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> vector_Vector3i;
        typedef std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vector_Vector3d;
        typedef std::vector<Eigen::Vector4f, Eigen::aligned_allocator<Eigen::Vector4f>> vector_Vector4f;

        typedef std::vector<Eigen::ArrayXXd, Eigen::aligned_allocator<Eigen::ArrayXXd>> vector_ArrayXXd;

        typedef std::vector<Eigen::Matrix3Xd, Eigen::aligned_allocator<Eigen::Matrix3Xd>> vector_Matrix3Xd;

        typedef std::vector<Eigen::RowVectorXd, Eigen::aligned_allocator<Eigen::RowVectorXd>> vector_RowVectorXd;
    } // namespace EigenSTL
} // namespace xgandalf
#endif /* EIGENSTLCONTAINERS_H_ */
