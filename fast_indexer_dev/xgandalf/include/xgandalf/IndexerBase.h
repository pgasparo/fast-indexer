/*
 * IndexerBase.h
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

#ifndef INDEXERBASE_H_
#define INDEXERBASE_H_

#include "Dbscan.h"
#include "DetectorToReciprocalSpaceTransform.h"
#include "ExperimentSettings.h"
#include "Lattice.h"
#include "LatticeAssembler.h"
#include "SamplePointsGenerator.h"
#include "SparsePeakFinder.h"
#include <Eigen/Dense>

namespace xgandalf
{
    class IndexerBase
    {
      public:
        IndexerBase(const ExperimentSettings& experimentSettings);

        virtual ~IndexerBase() = default;

        virtual void index(std::vector<Lattice>& assembledLattices, const Eigen::Matrix3Xf& reciprocalPeaks_1_per_A) = 0;

      protected:
        void keepSamplePointsWithHighEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation, float minEvaluation);
        void keepSamplePointsWithHighestEvaluation(Eigen::Matrix3Xf& samplePoints, Eigen::RowVectorXf& samplePointsEvaluation,
                                                   uint32_t maxToTakeCount); // output is sorted

        ExperimentSettings experimentSettings;
        SamplePointsGenerator samplePointsGenerator;

        LatticeAssembler latticeAssembler;

      private:
        std::vector<uint32_t> sortIndices; // to avoid frequent reallocation
    };
} // namespace xgandalf
#endif /* INDEXERBASE_H_ */
