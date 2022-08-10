/*
 * main.cpp
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

#include <cmath>

#include <xgandalf/HillClimbingOptimizer.h>
#include <xgandalf/InverseSpaceTransform.h>
#include <xgandalf/Lattice.h>
#include <xgandalf/LatticeAssembler.h>
#include <xgandalf/SamplePointsGenerator.h>
#include <xgandalf/SparsePeakFinder.h>
#include <xgandalf/eigenDiskImport.h>
#include <xgandalf/pointAutocorrelation.h>
#include <xgandalf/tests.h>

#include <Eigen/Dense>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace xgandalf;
using namespace std;
using namespace Eigen;

int main()
{
    try
    {
        // testPatternPrediction();
        // test_filterSamplePointsForNorm();
        // test_indexerAutocorrPrefit();
        // test_indexerPlain();
        // test_crystfelAdaption();
        // test_crystfelAdaption2();
        // test_latticeReorder();
        // test_gradientDescentRefinement();
        // test_mixedGradientDescentRefinement();
        // test_fixedBasisRefinement();
        // test_fixedBasisRefinementKabsch();
        test_hillClimbing();
    }
    catch (exception& e)
    {
        cout << e.what();
    }

    cout << endl << "done";
    getchar();

    return 0;
}
