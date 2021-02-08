/*
    Residue.h
    =========
        Class Residue header.
*/

#ifndef __PDBTOOLS_RESIDUE_H
#define __PDBTOOLS_RESIDUE_H

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
    Residue *Copy() const;


    // GetResidues
    vector<Residue *> GetResidues();
    vector<const Residue *> GetResidues() const;


    // GetAtoms
    vector<Atom *> GetAtoms();
    vector<const Atom *> GetAtoms() const;


    // compNum
    string compNum() const;
    Residue *compNum(int resNum, const string &resIns = "");


    // subMap
    unordered_map<string, Atom *> subMap();
    unordered_map<string, const Atom *> subMap() const;


    // coordMap
    unordered_map<string, RowVector3d *> coordMap();
    unordered_map<string, const RowVector3d *> coordMap() const;


    // Calc Backbone Dihedral Angle
    double CalcBBDihedralAngle(DIH dihedralEnum) const;


    // Calc Backbone Rotation Matrix By Delta Angle
    Residue *CalcBBRotationMatrixByDeltaAngle(DIH dihedralEnum,
        SIDE sideEnum, double deltaAngle, RowVector3d &moveCoord,
        Matrix3d &rotationMatrix);

    const Residue *CalcBBRotationMatrixByDeltaAngle(DIH dihedralEnum,
        SIDE sideEnum, double deltaAngle, RowVector3d &moveCoord,
        Matrix3d &rotationMatrix) const;


    // Calc Backbone Rotation Matrix By Target Angle
    Residue *CalcBBRotationMatrixByTargetAngle(DIH dihedralEnum,
        SIDE sideEnum, double targetAngle, RowVector3d &moveCoord,
        Matrix3d &rotationMatrix);

    const Residue *CalcBBRotationMatrixByTargetAngle(DIH dihedralEnum,
        SIDE sideEnum, double targetAngle, RowVector3d &moveCoord,
        Matrix3d &rotationMatrix) const;


    // Get Backbone Rotation Atom Pointer
    vector<Atom *> GetBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum);
    vector<const Atom *> GetBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum) const;


    // Rotate Backbone Dihedral Angle By Delta Angle
    Residue *RotateBBDihedralAngleByDeltaAngle(DIH dihedralEnum,
        SIDE sideEnum, double deltaAngle);


    // Rotate Backbone Dihedral Angle By Target Angle
    Residue *RotateBBDihedralAngleByTargetAngle(DIH dihedralEnum,
        SIDE sideEnum, double targetAngle);


    // Calc Side Chain Dihedral Angle
    double CalcSCDihedralAngle(int dihedralIdx) const;


    // Calc Side Chain Rotation Matrix By Delta Angle
    Residue *CalcSCRotationMatrixByDeltaAngle(int dihedralIdx, double deltaAngle,
        RowVector3d &moveCoord, Matrix3d &rotationMatrix);

    const Residue *CalcSCRotationMatrixByDeltaAngle(int dihedralIdx, double deltaAngle,
        RowVector3d &moveCoord, Matrix3d &rotationMatrix) const;


    // Calc Side Chain Rotation Matrix By Target Angle
    Residue *CalcSCRotationMatrixByTargetAngle(int dihedralIdx, double targetAngle,
        RowVector3d &moveCoord, Matrix3d &rotationMatrix);

    const Residue *CalcSCRotationMatrixByTargetAngle(int dihedralIdx, double targetAngle,
        RowVector3d &moveCoord, Matrix3d &rotationMatrix) const;


    // Get Side Chain Rotation Atom Pointer
    vector<Atom *> GetSCRotationAtomPtr(int dihedralIdx);
    vector<const Atom *> GetSCRotationAtomPtr(int dihedralIdx) const;


    // Rotate Side Chain Dihedral Angle By Delta Angle
    Residue *RotateSCDihedralAngleByDeltaAngle(int dihedralIdx, double deltaAngle);


    // Rotate Side Chain Dihedral Angle By Target Angle
    Residue *RotateSCDihedralAngleByTargetAngle(int dihedralIdx, double targetAngle);


    // Destructor
    ~Residue();
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_RESIDUE_H
