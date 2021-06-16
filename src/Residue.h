/*
    Residue.h
    =========
        Class Residue header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include "StructBase.h"
#include "NotAtom.h"
#include "NotProtein.h"
#include "Chain.h"
#include "Atom.h"
#include "Constants.hpp"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::vector;
using std::unordered_map;
using Eigen::RowVector3d;
using Eigen::Matrix3d;


////////////////////////////////////////////////////////////////////////////////
// Class Residue
////////////////////////////////////////////////////////////////////////////////

class Residue: public __StructBase<Residue>, public __NotAtom<Residue, Atom>,
    public __NotProtein<Residue, Chain>
{
public:

    // Attribute
    string name;
    int num;
    string ins;
    Chain *owner;
    vector<Atom *> sub;


    // Constructor
    explicit Residue(const string &resName = "", int resNum = 0,
        const string &resIns = "", Chain *resOwner = nullptr);


    // str
    string str() const;


    // Copy
    Residue *Copy();


    // GetResidues
    vector<Residue *> GetResidues();


    // GetAtoms
    vector<Atom *> GetAtoms();


    // compNum
    string compNum();
    Residue *compNum(int resNum, const string &resIns = "");


    // subMap
    unordered_map<string, Atom *> subMap();


    // coordMap
    unordered_map<string, RowVector3d> coordMap();


    // Calc Backbone Dihedral Angle
    double CalcBBDihedralAngle(DIH dihedralEnum);


    // Calc Backbone Rotation Matrix By Delta Angle
    Residue *CalcBBRotationMatrixByDeltaAngle(
        RowVector3d &moveCoord, Matrix3d &rotationMatrix,
        DIH dihedralEnum, SIDE sideEnum, double deltaAngle);


    // Calc Backbone Rotation Matrix By Target Angle
    Residue *CalcBBRotationMatrixByTargetAngle(
        RowVector3d &moveCoord, Matrix3d &rotationMatrix,
        DIH dihedralEnum, SIDE sideEnum, double targetAngle);


    // Get Backbone Rotation Atom Pointer
    vector<Atom *> GetBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum);


    // Rotate Backbone Dihedral Angle By Delta Angle
    Residue *RotateBBDihedralAngleByDeltaAngle(DIH dihedralEnum,
        SIDE sideEnum, double deltaAngle);


    // Rotate Backbone Dihedral Angle By Target Angle
    Residue *RotateBBDihedralAngleByTargetAngle(DIH dihedralEnum,
        SIDE sideEnum, double targetAngle);


    // Calc Side Chain Dihedral Angle
    double CalcSCDihedralAngle(int dihedralIdx);


    // Calc Side Chain Rotation Matrix By Delta Angle
    Residue *CalcSCRotationMatrixByDeltaAngle(
        RowVector3d &moveCoord, Matrix3d &rotationMatrix,
        int dihedralIdx, double deltaAngle);


    // Calc Side Chain Rotation Matrix By Target Angle
    Residue *CalcSCRotationMatrixByTargetAngle(
        RowVector3d &moveCoord, Matrix3d &rotationMatrix,
        int dihedralIdx, double targetAngle);


    // Get Side Chain Rotation Atom Pointer
    vector<Atom *> GetSCRotationAtomPtr(int dihedralIdx);


    // Rotate Side Chain Dihedral Angle By Delta Angle
    Residue *RotateSCDihedralAngleByDeltaAngle(int dihedralIdx, double deltaAngle);


    // Rotate Side Chain Dihedral Angle By Target Angle
    Residue *RotateSCDihedralAngleByTargetAngle(int dihedralIdx, double targetAngle);


    // Destructor
    ~Residue();
};


}  // End namespace PDBTools
