/*
    Residue.h
    =========
        Class Residue header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <Eigen/Dense>
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
using std::pair;
using Eigen::RowVector3d;
using Eigen::Matrix3d;


////////////////////////////////////////////////////////////////////////////////
// Class Residue
////////////////////////////////////////////////////////////////////////////////

class Residue: public __NotAtom<Residue, Atom>,
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
    Residue *copy();


    // GetResidues
    vector<Residue *> getResidues();


    // GetAtoms
    vector<Atom *> getAtoms();


    // compNum
    string compNum();
    Residue *compNum(int resNum, const string &resIns = "");
    Residue *compNum(const pair<int, string> &compNumPair);


    // subMap
    unordered_map<string, Atom *> subMap();


    // coordMap
    unordered_map<string, RowVector3d> coordMap();


    // Calc Backbone Dihedral Angle
    double calcBBDihedralAngle(DIH dihedralEnum);


    // Calc Backbone Rotation Matrix By Delta Angle
    pair<RowVector3d, Matrix3d> calcBBRotationMatrixByDeltaAngle(
        DIH dihedralEnum, SIDE sideEnum, double deltaAngle);


    // Calc Backbone Rotation Matrix By Target Angle
    pair<RowVector3d, Matrix3d> calcBBRotationMatrixByTargetAngle(
        DIH dihedralEnum, SIDE sideEnum, double targetAngle);


    // Get Backbone Rotation Atom Pointer
    vector<Atom *> getBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum);


    // Rotate Backbone Dihedral Angle By Delta Angle
    Residue *rotateBBDihedralAngleByDeltaAngle(DIH dihedralEnum,
        SIDE sideEnum, double deltaAngle);


    // Rotate Backbone Dihedral Angle By Target Angle
    Residue *rotateBBDihedralAngleByTargetAngle(DIH dihedralEnum,
        SIDE sideEnum, double targetAngle);


    // Calc Side Chain Dihedral Angle
    double calcSCDihedralAngle(int dihedralIdx);


    // Calc Side Chain Rotation Matrix By Delta Angle
    pair<RowVector3d, Matrix3d> calcSCRotationMatrixByDeltaAngle(
        int dihedralIdx, double deltaAngle);


    // Calc Side Chain Rotation Matrix By Target Angle
    pair<RowVector3d, Matrix3d> calcSCRotationMatrixByTargetAngle(
        int dihedralIdx, double targetAngle);


    // Get Side Chain Rotation Atom Pointer
    vector<Atom *> getSCRotationAtomPtr(int dihedralIdx);


    // Rotate Side Chain Dihedral Angle By Delta Angle
    Residue *rotateSCDihedralAngleByDeltaAngle(int dihedralIdx, double deltaAngle);


    // Rotate Side Chain Dihedral Angle By Target Angle
    Residue *rotateSCDihedralAngleByTargetAngle(int dihedralIdx, double targetAngle);


    // Dump
    Residue *dump(const string &dumpFilePath, const string &fileMode = "w");


    // Destructor
    ~Residue();
};


}  // End namespace PDBTools
