/*
    MathUtil.hpp
    ============
        Math utils implementation.
*/

#ifndef __PDBTOOLS_MATH_UTIL_HPP
#define __PDBTOOLS_MATH_UTIL_HPP

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

inline double Degrees(double radiansAngle)
{
    return radiansAngle * 180. / M_PI;
}


////////////////////////////////////////////////////////////////////////////////
// Convert Degrees To Radians
////////////////////////////////////////////////////////////////////////////////

inline double Radians(double degreesAngle)
{
    return degreesAngle * M_PI / 180.;
}


////////////////////////////////////////////////////////////////////////////////
// Calc Vector Angle
////////////////////////////////////////////////////////////////////////////////

inline double CalcVectorAngle(const RowVector3d &coordA, const RowVector3d &coordB)
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
// Calc Rotation Matrix By Two Vector (From A To B, Right Multiply Matrix)
////////////////////////////////////////////////////////////////////////////////

inline Matrix3d CalcRotationMatrixByTwoVector(const RowVector3d &coordA,
    const RowVector3d &coordB)
{
    return CalcRotationMatrix(coordA.cross(coordB), CalcVectorAngle(coordA, coordB));
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

    if (rotationAngle != 0.)
    {
        Matrix3d rotationMatrix = CalcRotationMatrix(
            OC.cross(RowVector3d(1., 0., 0.)), rotationAngle);

        OA *= rotationMatrix;
        OD *= rotationMatrix;
    }

    OA[0] = 0.;
    OD[0] = 0.;

    rotationAngle = CalcVectorAngle(OA, RowVector3d(0., 0., 1.));

    if (rotationAngle != 0.)
    {
        OD *= CalcRotationMatrix(OA.cross(RowVector3d(0., 0., 1.)), rotationAngle);
    }

    if (OD[1] > 0.)
    {
        dihedralAngle = -dihedralAngle;
    }

    return dihedralAngle;
}


////////////////////////////////////////////////////////////////////////////////
// Calc RMSD (Root-Mean-Square Deviation)
////////////////////////////////////////////////////////////////////////////////

inline double CalcRMSD(const Matrix<double, Dynamic, 3> &coordArrayA,
    const Matrix<double, Dynamic, 3> &coordArrayB)
{
    return sqrt((coordArrayA - coordArrayB).array().square().sum() / coordArrayA.rows());
}


////////////////////////////////////////////////////////////////////////////////
// Calc Superimpose Rotation Matrix (Kabsch Algorithm)
////////////////////////////////////////////////////////////////////////////////

void CalcSuperimposeRotationMatrix(const Matrix<double, Dynamic, 3> &sourceCoordArray,
    const Matrix<double, Dynamic, 3> &targetCoordArray, RowVector3d &sourceCenterCoord,
    Matrix3d &rotationMatrix, RowVector3d &targetCenterCoord)
{
    sourceCenterCoord = sourceCoordArray.colwise().mean();
    targetCenterCoord = targetCoordArray.colwise().mean();

    JacobiSVD<Matrix3d> svd(
        (sourceCoordArray.rowwise() - sourceCenterCoord).transpose() *
        (targetCoordArray.rowwise() - targetCenterCoord),
        ComputeFullU | ComputeFullV);

    Matrix3d U = svd.matrixU(), V = svd.matrixV().transpose();

    if (U.determinant() * V.determinant() < 0.)
    {
        U.col(2) = -U.col(2);
    }

    rotationMatrix = U * V;
}


////////////////////////////////////////////////////////////////////////////////
// Calc RMSD After Superimpose A => B
////////////////////////////////////////////////////////////////////////////////

double CalcRMSDAfterSuperimpose(const Matrix<double, Dynamic, 3> &coordArrayA,
    const Matrix<double, Dynamic, 3> &coordArrayB)
{
    RowVector3d sourceCenterCoord;
    Matrix3d rotationMatrix;
    RowVector3d targetCenterCoord;

    CalcSuperimposeRotationMatrix(coordArrayA, coordArrayB, sourceCenterCoord,
        rotationMatrix, targetCenterCoord);

    return CalcRMSD((((coordArrayA.rowwise() - sourceCenterCoord) * rotationMatrix).rowwise() +
        targetCenterCoord).eval(), coordArrayB);
}


////////////////////////////////////////////////////////////////////////////////
// Transform vector<RowVector3d *> To Matrix<double, Dynamic, 3>
////////////////////////////////////////////////////////////////////////////////

Matrix<double, Dynamic, 3> Coord2Matrix(const vector<RowVector3d *> &coordPtrList)
{
    Matrix<double, Dynamic, 3> coordMatrix(coordPtrList.size(), 3);

    for (int idx = 0; idx < coordPtrList.size(); idx++)
    {
        coordMatrix.row(idx) = *coordPtrList[idx];
    }

    return coordMatrix;
}


Matrix<double, Dynamic, 3> Coord2Matrix(const vector<const RowVector3d *> &coordPtrList)
{
    Matrix<double, Dynamic, 3> coordMatrix(coordPtrList.size(), 3);

    for (int idx = 0; idx < coordPtrList.size(); idx++)
    {
        coordMatrix.row(idx) = *coordPtrList[idx];
    }

    return coordMatrix;
}


}  // End namespace PDBTools


#endif  // __PDBTOOLS_MATH_UTIL_HPP
