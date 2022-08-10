/*
 * Lattice.h
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

#ifndef LATTICE_H_
#define LATTICE_H_

#include <Eigen/Dense>
#include <iostream>

namespace xgandalf
{
    class Lattice
    {
      public:
        Lattice();
        Lattice(const Eigen::Matrix3f& basis);
        Lattice(const Eigen::Vector3f& a, const Eigen::Vector3f& b, const Eigen::Vector3f& c);

        Lattice& minimize();

        inline float det() const
        {
            return basis.determinant();
        }

        inline const Eigen::Matrix3f& getBasis() const
        {
            return basis;
        }

        inline Eigen::Vector3f getBasisVectorNorms() const
        {
            return basis.colwise().norm();
        }

        Eigen::Vector3f getBasisVectorAngles_deg() const;
        Eigen::Vector3f getBasisVectorAnglesNormalized_deg() const;

        inline Lattice getReciprocalLattice() const
        {
            return Lattice(basis.transpose().inverse().eval());
        }

        friend std::ostream& operator<<(std::ostream& os, const Lattice& lattice);

        void reorder(const Eigen::Vector3f prototypeNorms, const Eigen::Vector3f prototypeAngles_deg);
        void reorder(const Lattice prototypeLattice);
        void normalizeAngles();

      private:
        Eigen::Matrix3f basis;

      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
} // namespace xgandalf
#endif /* LATTICE_H_ */
