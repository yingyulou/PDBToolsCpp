/*
    Math.hpp
    ========
        Math functions implementation.
*/

#pragma once

#include <cmath>
#include <vector>
#include <Eigen/Dense>

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::vector;
using Eigen::RowVector3d;
using Eigen::Matrix3d;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::JacobiSVD;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;


////////////////////////////////////////////////////////////////////////////////
// Convert Radians To Degrees
////////////////////////////////////////////////////////////////////////////////

double Degrees(double radiansAngle)
{
    return radiansAngle * 180. / M_PI;
}


////////////////////////////////////////////////////////////////////////////////
// Convert Degrees To Radians
////////////////////////////////////////////////////////////////////////////////

double Radians(double degreesAngle)
{
    return degreesAngle * M_PI / 180.;
}


////////////////////////////////////////////////////////////////////////////////
// Calc Vector Angle
////////////////////////////////////////////////////////////////////////////////

double CalcVectorAngle(const RowVector3d &coordA, const RowVector3d &coordB)
{
    return acos(coordA.dot(coordB) / (coordA.norm() * coordB.norm()));
}


////////////////////////////////////////////////////////////////////////////////
// Calc Rotation Matrix (Right Multiply Matrix)
////////////////////////////////////////////////////////////////////////////////

Matrix3d CalcRotationMatrix(const RowVector3d &rotationAxis, double rotationAngle)
{
    auto normRotationAxis = rotationAxis.normalized();

    double x = normRotationAxis[0], y = normRotationAxis[1], z = normRotationAxis[2],
        s = sin(rotationAngle), c = cos(rotationAngle), one_c = 1 - c;

    return (Matrix3d() <<
        c + x * x * one_c, x * y * one_c + z * s, x * z * one_c - y * s,
        x * y * one_c - z * s, c + y * y * one_c, y * z * one_c + x * s,
        x * z * one_c + y * s, y * z * one_c - x * s, c + z * z * one_c).finished();
}


////////////////////////////////////////////////////////////////////////////////
// Calc Rotation Matrix By Two Vector (Right Multiply Matrix)
////////////////////////////////////////////////////////////////////////////////

Matrix3d CalcRotationMatrixByTwoVector(const RowVector3d &refCoord,
    const RowVector3d &tarCoord)
{
    return CalcRotationMatrix(tarCoord.cross(refCoord),
        CalcVectorAngle(tarCoord, refCoord));
}


////////////////////////////////////////////////////////////////////////////////
// Calc Dihedral Angle
////////////////////////////////////////////////////////////////////////////////

double CalcDihedralAngle(const RowVector3d &coordA, const RowVector3d &coordB,
    const RowVector3d &coordC, const RowVector3d &coordD)
{
    RowVector3d AB = coordB - coordA;
    RowVector3d AC = coordC - coordA;
    RowVector3d DB = coordB - coordD;
    RowVector3d DC = coordC - coordD;

    RowVector3d ABAC = AB.cross(AC);
    RowVector3d DBDC = DB.cross(DC);

    double dihedralAngle = CalcVectorAngle(ABAC, DBDC);

    // Calc Sign
    RowVector3d OA = coordA - coordB;
    RowVector3d OC = coordC - coordB;
    RowVector3d OD = coordD - coordB;

    double rotationAngle = CalcVectorAngle(OC, RowVector3d(1., 0., 0.));

    Matrix3d rotationMatrix = CalcRotationMatrix(
        OC.cross(RowVector3d(1., 0., 0.)), rotationAngle);

    OA *= rotationMatrix;
    OD *= rotationMatrix;

    OA[0] = 0.;
    OD[0] = 0.;

    rotationAngle  = CalcVectorAngle(OA, RowVector3d(0., 0., 1.));
    rotationMatrix = CalcRotationMatrix(OA.cross(RowVector3d(0., 0., 1.)), rotationAngle);

    OD *= rotationMatrix;

    if (OD[1] > 0.)
    {
        dihedralAngle = -dihedralAngle;
    }

    return dihedralAngle;
}


////////////////////////////////////////////////////////////////////////////////
// Calc RMSD (Root-Mean-Square Deviation)
////////////////////////////////////////////////////////////////////////////////

double CalcRMSD(const Matrix<double, Dynamic, 3> &coordArrayA,
    const Matrix<double, Dynamic, 3> &coordArrayB)
{
    return sqrt((coordArrayA - coordArrayB).array().square().sum() / coordArrayA.rows());
}


////////////////////////////////////////////////////////////////////////////////
// Calc Superimpose Rotation Matrix (Kabsch Algorithm)
////////////////////////////////////////////////////////////////////////////////

void CalcSuperimposeRotationMatrix(RowVector3d &refCenterCoord,
    Matrix3d &rotationMatrix, RowVector3d &tarCenterCoord,
    const Matrix<double, Dynamic, 3> &refCoordArray,
    const Matrix<double, Dynamic, 3> &tarCoordArray)
{
    refCenterCoord = refCoordArray.colwise().mean();
    tarCenterCoord = tarCoordArray.colwise().mean();

    JacobiSVD<Matrix3d> svd(
        (tarCoordArray.rowwise() - tarCenterCoord).transpose() *
        (refCoordArray.rowwise() - refCenterCoord),
        ComputeFullU | ComputeFullV);

    Matrix3d U = svd.matrixU(), V = svd.matrixV().transpose();

    if (U.determinant() * V.determinant() < 0.)
    {
        U.col(2) = -U.col(2);
    }

    rotationMatrix = U * V;
}


////////////////////////////////////////////////////////////////////////////////
// Calc RMSD After Superimpose A <= B
////////////////////////////////////////////////////////////////////////////////

double CalcRMSDAfterSuperimpose(const Matrix<double, Dynamic, 3> &refCoordArray,
    const Matrix<double, Dynamic, 3> &tarCoordArray)
{
    RowVector3d refCenterCoord, tarCenterCoord;
    Matrix3d rotationMatrix;

    CalcSuperimposeRotationMatrix(refCenterCoord, rotationMatrix, tarCenterCoord,
        refCoordArray, tarCoordArray);

    return CalcRMSD(refCoordArray, (((tarCoordArray.rowwise() - tarCenterCoord) *
        rotationMatrix).rowwise() + refCenterCoord).eval());
}


}  // End namespace PDBTools
