/*
    Math.hpp
    ========
        Math functions implementation.
*/

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
#include <Eigen/Dense>

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using std::vector;
using std::min;
using std::max;
using std::tuple;
using Eigen::RowVector3d;
using Eigen::Matrix3d;
using Eigen::MatrixX3d;
using Eigen::JacobiSVD;
using Eigen::ComputeFullU;
using Eigen::ComputeFullV;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert Radians To Degrees
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double degrees(double radiansAngle)
{
    return radiansAngle * 180. / M_PI;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert Degrees To Radians
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double radians(double degreesAngle)
{
    return degreesAngle * M_PI / 180.;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc Vector Angle
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calcVectorAngle(const RowVector3d &coordA, const RowVector3d &coordB)
{
    return acos(min(max(coordA.dot(coordB) / (coordA.norm() * coordB.norm()), -1.), 1.));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc Rotation Matrix (Right Multiply Matrix)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix3d calcRotationMatrix(const RowVector3d &rotationAxis, double rotationAngle)
{
    auto normRotationAxis = rotationAxis.normalized();

    double x = normRotationAxis[0], y = normRotationAxis[1], z = normRotationAxis[2],
        s = sin(rotationAngle), c = cos(rotationAngle), _1c = 1. - c;

    return (Matrix3d() <<
        c + x * x * _1c, x * y * _1c + z * s, x * z * _1c - y * s,
        x * y * _1c - z * s, c + y * y * _1c, y * z * _1c + x * s,
        x * z * _1c + y * s, y * z * _1c - x * s, c + z * z * _1c).finished();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc Rotation Matrix By Two Vector (Right Multiply Matrix)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix3d calcRotationMatrixByTwoVector(const RowVector3d &tarCoord, const RowVector3d &srcCoord)
{
    return calcRotationMatrix(srcCoord.cross(tarCoord), calcVectorAngle(srcCoord, tarCoord));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc Dihedral Angle
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calcDihedralAngle(const RowVector3d &coordA, const RowVector3d &coordB, const RowVector3d &coordC, const RowVector3d &coordD)
{
    RowVector3d AB = coordB - coordA;
    RowVector3d AC = coordC - coordA;
    RowVector3d DB = coordB - coordD;
    RowVector3d DC = coordC - coordD;

    RowVector3d ABAC = AB.cross(AC);
    RowVector3d DBDC = DB.cross(DC);

    double dihedralAngle = calcVectorAngle(ABAC, DBDC);

    // Calc Sign
    RowVector3d OA = coordA - coordB;
    RowVector3d OC = coordC - coordB;
    RowVector3d OD = coordD - coordB;

    double rotationAngle = calcVectorAngle(OC, RowVector3d(1., 0., 0.));

    Matrix3d rotationMatrix = calcRotationMatrix(OC.cross(RowVector3d(1., 0., 0.)), rotationAngle);

    OA *= rotationMatrix;
    OD *= rotationMatrix;

    OA[0] = 0.;
    OD[0] = 0.;

    rotationAngle  = calcVectorAngle(OA, RowVector3d(0., 0., 1.));
    rotationMatrix = calcRotationMatrix(OA.cross(RowVector3d(0., 0., 1.)), rotationAngle);

    OD *= rotationMatrix;

    if (OD[1] > 0.)
    {
        dihedralAngle = -dihedralAngle;
    }

    return dihedralAngle;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc RMSD (Root-Mean-Square Deviation)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calcRMSD(const MatrixX3d &coordMatrixA, const MatrixX3d &coordMatrixB)
{
    return sqrt((coordMatrixA - coordMatrixB).array().square().sum() / coordMatrixA.rows());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc Superimpose Rotation Matrix (Kabsch Algorithm)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tuple<RowVector3d, Matrix3d, RowVector3d> calcSuperimposeRotationMatrix(
    const MatrixX3d &tarCoordMatrix, const MatrixX3d &srcCoordMatrix)
{
    RowVector3d srcCenterCoord = srcCoordMatrix.colwise().mean();
    RowVector3d tarCenterCoord = tarCoordMatrix.colwise().mean();

    JacobiSVD<Matrix3d> svd((srcCoordMatrix.rowwise() - srcCenterCoord).transpose() *
        (tarCoordMatrix.rowwise() - tarCenterCoord), ComputeFullU | ComputeFullV);

    Matrix3d U = svd.matrixU(), V = svd.matrixV().transpose();

    if (U.determinant() * V.determinant() < 0.)
    {
        U.col(2) = -U.col(2);
    }

    auto rotationMatrix = U * V;

    return {srcCenterCoord, rotationMatrix, tarCenterCoord};
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calc RMSD After Superimpose A <= B
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calcRMSDAfterSuperimpose(const MatrixX3d &tarCoordMatrix, const MatrixX3d &srcCoordMatrix)
{
    auto [srcCenterCoord, rotationMatrix, tarCenterCoord] = calcSuperimposeRotationMatrix(tarCoordMatrix, srcCoordMatrix);

    return calcRMSD(tarCoordMatrix, (((srcCoordMatrix.rowwise() - srcCenterCoord) * rotationMatrix).rowwise() + tarCenterCoord).eval());
}


}  // End namespace PDBTools
