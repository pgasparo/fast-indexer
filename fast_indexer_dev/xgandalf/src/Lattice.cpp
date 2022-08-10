/*
 * Lattice.cpp
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

#include <xgandalf/Lattice.h>
#include <algorithm>
#include <assert.h>
#include <inttypes.h>
#include <limits>
#include <math.h>

#ifndef M_1_PI
#define M_1_PI 0.318309886183790671538
#endif // !M_1_PI

using namespace std;
using namespace Eigen;

namespace xgandalf
{
    static inline void minimize2DLattice(Matrix3f& basis);
    static inline void sortTwoColumns_ascending(Matrix3f& basis, Vector3f& squaredNorms, float colNumber1, float colNumber2);
    static inline void sortColumnsByNorm_ascending(Matrix3f& basis);
    static inline void getMinA(Matrix3f& basis, Vector3f& minA, float& minALengthSquared);

    Lattice::Lattice() {}

    Lattice::Lattice(const Matrix3f& basis)
        : basis(basis)
    {
    }
    Lattice::Lattice(const Vector3f& a, const Vector3f& b, const Vector3f& c)
    {
        basis << a, b, c;
    }

    static inline void minimize2DLattice(Matrix3f& basis)
    {
        if (basis.col(0).squaredNorm() < basis.col(1).squaredNorm())
        { // will always be true in the current arrangement!
            basis.col(0).swap(basis.col(1));
        }

        do
        {
            float r = round(basis.col(0).dot(basis.col(1)) / basis.col(1).squaredNorm());
            Vector3f tmp = basis.col(0) - r * basis.col(1);

            basis.col(0) = basis.col(1);
            basis.col(1) = tmp;
        } while (basis.col(1).squaredNorm() < basis.col(0).squaredNorm());
    }

    static inline void sortTwoColumns_ascending(Matrix3f& basis, Vector3f& squaredNorms, float colNumber1, float colNumber2)
    {
        if (squaredNorms[colNumber1] > squaredNorms[colNumber2])
        {
            basis.col(colNumber1).swap(basis.col(colNumber2));

            float tmp = squaredNorms[colNumber1];
            squaredNorms[colNumber1] = squaredNorms[colNumber2];
            squaredNorms[colNumber2] = tmp;
        }
    }

    static inline void sortColumnsByNorm_ascending(Matrix3f& basis)
    {
        // simple bubble sort implementation
        Vector3f squaredNorms = basis.colwise().squaredNorm();
        sortTwoColumns_ascending(basis, squaredNorms, 0, 1);
        sortTwoColumns_ascending(basis, squaredNorms, 1, 2);
        sortTwoColumns_ascending(basis, squaredNorms, 0, 1);
    }

    static inline void getMinA(Matrix3f& basis, Vector3f& minA, float& minALengthSquared)
    {
        const Vector3f& b1 = basis.col(0);
        const Vector3f& b2 = basis.col(1);
        const Vector3f& b3 = basis.col(2);

        float b232 = b2.dot(b3) / b2.squaredNorm();
        float b122 = b1.dot(b2) / b2.squaredNorm();
        float b131 = b1.dot(b3) / b1.squaredNorm();
        float b121 = b1.dot(b2) / b1.squaredNorm();

        float y2 = -(b232 - b122 * b131) / (1 - b121 * b122);
        float y1 = -(b131 - b121 * b232) / (1 - b121 * b122);

        minALengthSquared = numeric_limits<float>::max();
        for (float i_1 = -1; i_1 <= 1; i_1++)
        {
            for (float i_2 = -1; i_2 <= 1; i_2++)
            {
                float x1 = roundf(y1 + i_1);
                float x2 = roundf(y2 + i_2);
                if ((abs(x1 - y1) <= 1) & (abs(x2 - y2) <= 1))
                {
                    Vector3f a = b3 + x2 * b2 + x1 * b1;
                    const float aLengthSquared = a.squaredNorm();

                    if (aLengthSquared < minALengthSquared)
                    {
                        minA = a;
                        minALengthSquared = aLengthSquared;
                    }
                }
            }
        }
    }

    // algorithm implemented after http://www.csie.nuk.edu.tw/~cychen/Lattices/A%203-Dimensional%20Lattice%20Reduction%20Algorithm.pdf
    Lattice& Lattice::minimize()
    {
        Vector3f minA(0, 0, 0); // initialization not required, but if not done, compiler issues warning

        bool terminationFlag = false;
        while (terminationFlag == false)
        {
            sortColumnsByNorm_ascending(basis);

            minimize2DLattice(basis);

            float minALengthSquared;
            getMinA(basis, minA, minALengthSquared);

            if (minALengthSquared >= basis.col(2).squaredNorm())
            {
                terminationFlag = true;
            }
            else
            {
                basis.col(2) = minA;
            }
        }

        sortColumnsByNorm_ascending(basis);

        return *this;
    }

    Vector3f Lattice::getBasisVectorAngles_deg() const
    {
        const Vector3f& a = basis.col(0);
        const Vector3f& b = basis.col(1);
        const Vector3f& c = basis.col(2);

        Vector3f angles;
        angles[0] = atan2(b.cross(c).norm(), b.dot(c)) * (M_1_PI * 180);
        angles[1] = atan2(a.cross(c).norm(), a.dot(c)) * (M_1_PI * 180);
        angles[2] = atan2(a.cross(b).norm(), a.dot(b)) * (M_1_PI * 180);

        return angles;
    }

    Vector3f Lattice::getBasisVectorAnglesNormalized_deg() const
    {
        const Vector3f& a = basis.col(0);
        const Vector3f& b = basis.col(1);
        const Vector3f& c = basis.col(2);

        Vector3f angles;
        angles[0] = atan2(b.cross(c).norm(), abs(b.dot(c))) * (M_1_PI * 180);
        angles[1] = atan2(a.cross(c).norm(), abs(a.dot(c))) * (M_1_PI * 180);
        angles[2] = atan2(a.cross(b).norm(), abs(a.dot(b))) * (M_1_PI * 180);

        return angles;
    }

    std::ostream& operator<<(std::ostream& os, const Lattice& lattice)
    {
        os << lattice.basis;
        return os;
    }

    // reorder such, that reordered basis differs from the prototype lattice by a rotation
    void Lattice::reorder(const Lattice prototypeLattice)
    {
        Matrix3f prototypeBasis_inv = prototypeLattice.getBasis().inverse();

        // clang-format off
	Array<int, 3, 6> indexPermutations;
	indexPermutations <<
		0, 0, 1, 1, 2, 2,
		1, 2, 0, 2, 0, 1,
		2, 1, 2, 0, 1, 0;

	Array<float, 3, 8> signPermutations;
	signPermutations <<
		1,  1,  1,  1, -1, -1, -1, -1,
		1,  1, -1, -1,  1,  1, -1, -1,
		1, -1,  1, -1,  1, -1,  1, -1;
    

    Matrix3f permutedBasis;
    Matrix3f bestPermutedBasis;
	Matrix3f resultingRotationMatrix;
	float smallestDefect = numeric_limits<float>::max();
    for (int indexPermutation = 0; indexPermutation < 6; indexPermutation++)
    {
        for (int signPermutation = 0; signPermutation < 8; signPermutation++)
        {
            permutedBasis << basis.col(indexPermutations(0, indexPermutation))*signPermutations(0, signPermutation),
					         basis.col(indexPermutations(1, indexPermutation))*signPermutations(1, signPermutation),
						     basis.col(indexPermutations(2, indexPermutation))*signPermutations(2, signPermutation);

			resultingRotationMatrix = permutedBasis*prototypeBasis_inv;

			//cout << "resultingRotationMatrix \n" << resultingRotationMatrix << endl;

			float defect = (resultingRotationMatrix.transpose() * resultingRotationMatrix - Matrix3f::Identity()).squaredNorm();
			//cout << "defect: " << defect<<endl;
			if (defect < smallestDefect) {
				bestPermutedBasis = permutedBasis;
				smallestDefect = defect;
				//cout << "taken\n" << permutedBasis 
				//	<< endl << endl<< prototypeLattice.getBasis() 
				//	<< endl << endl << prototypeBasis_inv 
				//	<< endl << endl << resultingRotationMatrix << endl << endl  ;
			}
        }
    }

	basis = bestPermutedBasis;

        // clang-format on
    }

    void Lattice::reorder(const Eigen::Vector3f prototypeNorms, const Eigen::Vector3f prototypeAngles_deg)
    {
        // cout << "prototype angles" << endl << prototypeAngles_deg << endl;

        Vector3f prototypeAnglesNormalized_deg = prototypeAngles_deg;
        for (int i = 0; i < 3; i++)
        {
            if (prototypeAnglesNormalized_deg(i) > 90)
            {
                prototypeAnglesNormalized_deg(i) = 180 - prototypeAnglesNormalized_deg(i);
            }
        }

        // cout << "prototype angles normalized" << endl << prototypeAnglesNormalized_deg << endl;

        Array<float, 6, 1> prototypeLatticeParameters;
        prototypeLatticeParameters << prototypeNorms, prototypeAnglesNormalized_deg;
        Array<float, 6, 1> prototypeLatticeParametersInverse = 1.0f / prototypeLatticeParameters;

        Vector3f n = getBasisVectorNorms();
        Vector3f a = getBasisVectorAnglesNormalized_deg();

        // clang-format off
	Array<int, 3, 6> indexPermutations;
	indexPermutations <<
		0, 0, 1, 1, 2, 2,
		1, 2, 0, 2, 0, 1,
		2, 1, 2, 0, 1, 0;

	Array<float, 6, 6> allPermutations;
	auto& i = indexPermutations;
	allPermutations <<
		n[i(0, 0)], n[i(0, 1)], n[i(0, 2)], n[i(0, 3)], n[i(0, 4)], n[i(0, 5)],
		n[i(1, 0)], n[i(1, 1)], n[i(1, 2)], n[i(1, 3)], n[i(1, 4)], n[i(1, 5)],
		n[i(2, 0)], n[i(2, 1)], n[i(2, 2)], n[i(2, 3)], n[i(2, 4)], n[i(2, 5)],
		a[i(0, 0)], a[i(0, 1)], a[i(0, 2)], a[i(0, 3)], a[i(0, 4)], a[i(0, 5)],
		a[i(1, 0)], a[i(1, 1)], a[i(1, 2)], a[i(1, 3)], a[i(1, 4)], a[i(1, 5)],
		a[i(2, 0)], a[i(2, 1)], a[i(2, 2)], a[i(2, 3)], a[i(2, 4)], a[i(2, 5)];
        // clang-format on

        // find best permutation
        auto rasiduals = ((allPermutations.colwise() - prototypeLatticeParameters).colwise() * prototypeLatticeParametersInverse).abs();
        auto maxResiduals = rasiduals.colwise().maxCoeff();

        int bestPermutationIndex;
        maxResiduals.minCoeff(&bestPermutationIndex);

        Matrix3f reorderedBasis;
        reorderedBasis << basis.col(indexPermutations(0, bestPermutationIndex)), basis.col(indexPermutations(1, bestPermutationIndex)),
            basis.col(indexPermutations(2, bestPermutationIndex));
        basis = reorderedBasis;

        // change sign of vectors
        Matrix3f bestNegatedBasis;
        float bestNegatedBasisResidual = numeric_limits<float>::max();
        for (int l = -1; l <= 1; l += 2)
        {
            for (int m = -1; m <= 1; m += 2)
            {
                Matrix3f negatedBasis;
                negatedBasis << basis.col(0), basis.col(1) * l, basis.col(2) * m;
                Lattice negatedLattice(negatedBasis);

                float maxResidual = (negatedLattice.getBasisVectorAngles_deg() - prototypeAngles_deg).cwiseAbs().maxCoeff();

                if (maxResidual < bestNegatedBasisResidual)
                {
                    bestNegatedBasisResidual = maxResidual;
                    bestNegatedBasis = negatedBasis;
                }
            }
        }
        basis = bestNegatedBasis;

        // cout << "best residual: " << bestNegatedBasisResidual << endl;

        // cout << "angles" << endl << getBasisVectorAngles_deg() << endl;
    }

    void Lattice::normalizeAngles()
    {
        for (int l = -1; l <= 1; l += 2)
        {
            for (int m = -1; m <= 1; m += 2)
            {
                Matrix3f negatedBasis;
                negatedBasis << basis.col(0), basis.col(1) * l, basis.col(2) * m;
                Lattice negatedLattice(negatedBasis);

                if ((negatedLattice.getBasisVectorAngles_deg().array() <= 90).all())
                {
                    basis = negatedBasis;
                    return;
                }
            }
        }

        // not possible to have all < 90°, get the combination with smallest sum
        float smallestSum = 3 * 180;
        for (int l = -1; l <= 1; l += 2)
        {
            for (int m = -1; m <= 1; m += 2)
            {
                Matrix3f negatedBasis;
                negatedBasis << basis.col(0), basis.col(1) * l, basis.col(2) * m;
                Lattice negatedLattice(negatedBasis);

                float anglesSum = negatedLattice.getBasisVectorAngles_deg().sum();
                if (anglesSum <= smallestSum)
                {
                    smallestSum = anglesSum;
                    basis = negatedBasis;
                }
            }
        }
    }
} // namespace xgandalf
