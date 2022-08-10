/*
 * samplePointsGenerator.cpp
 *
 * refinement.cpp
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

#include <xgandalf/refinement.h>
#include <iostream>

using namespace Eigen;
using namespace std;

namespace xgandalf
{
    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void getGradient_reciprocalPeakMatch_meanDist(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        Matrix3Xf d = B * M - N;
        RowVectorXf d_norms = d.colwise().norm();

        Array<float, 1, Dynamic> denominator_inv = 1 / (d_norms.array());

        for (int i = 0; i < denominator_inv.size(); ++i)
        {
            if (!isfinite(denominator_inv[i]))
            {
                denominator_inv = 1e-30;
            }
        }

        Array<float, 1, Dynamic> numerator;
        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                numerator = d.row(row).array() * M.row(col).array();

                gradient(row, col) = (numerator * denominator_inv).mean();
            }
        }
    }


    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    // TODO: Function can be made faster easily... nevertheless double is really necessary...
    void getGradient_detectorAngleMatch(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        Matrix2Xd detectorPeakDirections = N.bottomRows(2).colwise().normalized().cast<double>();

        Matrix3Xd Md = M.cast<double>();

        Matrix3Xd predictedPoints = B.cast<double>() * Md;
        RowVectorXd projectionNorms = (predictedPoints.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
        Matrix2Xd predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionNorms.array();
        RowVectorXd defect = (predictedPoints.bottomRows(2) - predictedPointsProjected).colwise().norm();
        double meanDefect = defect.mean();

        // now do the numeric differentiatiation
        double differentiatiationShift = 1e-8; // should be small enough, but not too small
        double differentiatiationShift_inv = 1 / differentiatiationShift;
        Matrix3Xd predictedPoints_shifted;
        for (int row = 1; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                predictedPoints_shifted = predictedPoints;
                predictedPoints_shifted.row(row) += differentiatiationShift * Md.row(col);


                projectionNorms = (predictedPoints_shifted.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
                predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionNorms.array();
                defect = (predictedPoints_shifted.bottomRows(2) - predictedPointsProjected).colwise().norm();
                double meanDefect_shifted = defect.mean();

                gradient(row, col) = (meanDefect_shifted - meanDefect) * differentiatiationShift_inv;
            }
        }

        // shifts in x direction do not change the angle
        gradient(0, 0) = 0;
        gradient(0, 1) = 0;
        gradient(0, 2) = 0;
    }

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void getGradient_reciprocalPeakMatch_meanSquaredDist(Matrix3f& gradient, const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        Matrix3Xf d = B * M - N;
        RowVectorXf d_norms = d.colwise().norm();

        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                gradient(row, col) = (d.row(row).array() * M.row(col).array()).sum();
            }
        }
    }

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanSquaredDist(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        B = M.transpose().colPivHouseholderQr().solve(N.transpose()).transpose();
    }

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanDist_peaksAndAngle(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        Matrix3f gradient;
        Matrix3f summedGradient;

        float reciprocalPeakDistWeight = 1;
        float reciprocalPeakAngleWeight = 0.5;

        float stepLength = B.maxCoeff() * 0.003;
        for (int i = 0; i < 150; i++)
        {
            getGradient_reciprocalPeakMatch_meanDist(gradient, B, M, N);
            summedGradient = gradient * reciprocalPeakDistWeight;
            getGradient_detectorAngleMatch(gradient, B, M, N);
            summedGradient += gradient * reciprocalPeakAngleWeight;
            float maxCoeff = summedGradient.cwiseAbs().maxCoeff();
            if (maxCoeff < 1e-10)
            {
                break;
            }
            summedGradient /= maxCoeff;
            B = B - stepLength * summedGradient;

            if (i >= 75)
            {
                stepLength *= 0.93;
            }
        }
    }

    // B: reciprocal basis
    // B_sample: sample reciprocal basis (with correct lattice parameters)
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanSquaredDist_fixedBasisParameters(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N, const Matrix3f& B_sample)
    {
        JacobiSVD<MatrixXf> svd(M.rows(), M.cols(), ComputeThinU | ComputeThinV);

        svd.compute(B_sample * M);
        MatrixXf S1 = svd.singularValues();
        Matrix3f U1 = svd.matrixU();
        MatrixXf V1 = svd.matrixV();

        svd.compute(N);
        MatrixXf S2 = svd.singularValues();
        Matrix3f U2 = svd.matrixU();
        MatrixXf V2 = svd.matrixV();

        Matrix3f permutationFlipMatrix = (V1.transpose() * V2).array().round(); // maybe not necessary
        U2 = U2 * permutationFlipMatrix;

        if (abs(permutationFlipMatrix.determinant()) != 1)
        {
            B.setIdentity();
            cout << "Not able to fit found basis to sample basis! \n" << (V1.transpose() * V2).eval() << endl;
            return;
        }

        Matrix3f estimatedR = U2 * U1.transpose();

        B = estimatedR * B_sample;
    }

    // B: reciprocal basis
    // B_sample: sample reciprocal basis (with correct lattice parameters)
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanSquaredDist_fixedBasisParameters_kabsch(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N, const Matrix3f& B_sample)
    {
        const Matrix3Xf& Q = N;
        Matrix3Xf P = B_sample * M;

        Matrix3f H = P * Q.transpose();

        JacobiSVD<MatrixXf> svd(H.rows(), H.cols(), ComputeFullU | ComputeFullV);
        svd.compute(H);
        Matrix3f U = svd.matrixU();
        Matrix3f V = svd.matrixV();

        float d = (V * U.transpose()).determinant();

        Matrix3f rightHandInsurance = Matrix3f::Identity();
        rightHandInsurance(2, 2) = d;

        Matrix3f estimatedR = V * rightHandInsurance * U.transpose();

        B = estimatedR * B_sample;
    }


    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    // TODO: Function can be made faster easily... nevertheless double is really necessary...
    void getGradient_detectorAngleMatchFixedLattice(Vector3f& rotationAnglesGradient_detectorAngle, Vector3f& rotationAnglesGradient_distFromPoints,
                                                    const Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        Matrix3Xd Md = M.cast<double>();
        Matrix3Xd Nd = N.cast<double>();
        Matrix3Xd Bd = B.cast<double>();

        double meanDefectInitial_distFromPoint = (Bd * Md - Nd).colwise().norm().mean();

        Matrix2Xd detectorPeakDirections = N.bottomRows(2).colwise().normalized().cast<double>();

        Matrix3Xd predictedPoints = Bd.cast<double>() * Md;
        RowVectorXd projectionNorms = (predictedPoints.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
        Matrix2Xd predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionNorms.array();
        RowVectorXd defect_detectorAngle = (predictedPoints.bottomRows(2) - predictedPointsProjected).colwise().norm();
        double meanDefectInitial_detectorAngle = defect_detectorAngle.mean();

        // now do the numeric differentiatiation
        double differentiatiationShift = 1e-8; // should be small enough, but not too small
        double differentiatiationShift_inv = 1 / differentiatiationShift;
        Matrix3Xd predictedPoints_shifted;
        for (int i = 0; i < 3; i++)
        {
            Vector3d axis(0, 0, 0);
            axis[i] = 1;
            predictedPoints_shifted = AngleAxisd(differentiatiationShift, axis).toRotationMatrix() * predictedPoints;

            projectionNorms = (predictedPoints_shifted.bottomRows(2).cwiseProduct(detectorPeakDirections)).colwise().sum(); // colwise dot product
            predictedPointsProjected = detectorPeakDirections.array().rowwise() * projectionNorms.array();
            defect_detectorAngle = (predictedPoints_shifted.bottomRows(2) - predictedPointsProjected).colwise().norm();

            rotationAnglesGradient_detectorAngle(i) = (defect_detectorAngle.mean() - meanDefectInitial_detectorAngle) * differentiatiationShift_inv;
            rotationAnglesGradient_distFromPoints(i) =
                ((predictedPoints_shifted - Nd).colwise().norm().mean() - meanDefectInitial_distFromPoint) * differentiatiationShift_inv;
        }
    }

    // B: reciprocal basis
    // M: miller indices (reciprocal)
    // N: reciprocal peaks
    void refineReciprocalBasis_meanDist_detectorAngleMatchFixedParameters(Matrix3f& B, const Matrix3Xf& M, const Matrix3Xf& N)
    {
        Vector3f rotationAnglesGradient_detectorAngle;
        Vector3f rotationAnglesGradient_distFromPoints;
        Vector3f summedGradient;

        float reciprocalPeakDistWeight = 1;
        float reciprocalPeakAngleWeight = 0.3;

        float stepLength = 0.05 / 180 * 3.14;
        for (int i = 0; i < 150; i++)
        {
            getGradient_detectorAngleMatchFixedLattice(rotationAnglesGradient_detectorAngle, rotationAnglesGradient_distFromPoints, B, M, N);
            summedGradient =
                rotationAnglesGradient_detectorAngle * reciprocalPeakAngleWeight + rotationAnglesGradient_distFromPoints * reciprocalPeakDistWeight;

            float maxCoeff = summedGradient.cwiseAbs().maxCoeff();
            if (maxCoeff < 1e-10)
            {
                break;
            }
            summedGradient /= maxCoeff;
            summedGradient *= -stepLength;
            Matrix3f rotation;
            rotation = AngleAxisf(summedGradient[0], Vector3f::UnitX()) * AngleAxisf(summedGradient[1], Vector3f::UnitY()) *
                       AngleAxisf(summedGradient[2], Vector3f::UnitZ());

            B = rotation * B;

            if (i >= 75)
            {
                stepLength *= 0.93;
            }
        }
    }
} // namespace xgandalf
