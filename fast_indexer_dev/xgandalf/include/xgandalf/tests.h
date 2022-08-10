/*
 * tests.h
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

#ifndef TESTS_H_
#define TESTS_H_

namespace xgandalf
{

    void test();
    void testPatternPrediction();
    void test_fixedBasisRefinementKabsch();
    void test_fixedBasisRefinement();
    void test_mixedGradientDescentRefinement();
    void test_gradientDescentRefinement();
    void test_getGradient();
    void test_latticeReorder();
    void test_crystfelAdaption2();
    void test_crystfelAdaption();
    void test_filterSamplePointsForNorm();
    void test_indexerAutocorrPrefit();
    void test_indexerPlain();
    void test_dbscan();
    void test_pointAutocorrelation();
    void test_latticeAssembler();
    void test_sparsePeakFinder();
    void test_hillClimbing();
    void test_computeStep();
    void test_InverseSpaceTransform();

} // namespace xgandalf

#endif /* TESTS_H_ */
