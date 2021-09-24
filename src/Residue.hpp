/*
    Residue.hpp
    ===========
        Class Residue implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cctype>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "Residue.h"
#include "Chain.h"
#include "Atom.h"
#include "Math.hpp"
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
using boost::format;
using Eigen::RowVector3d;
using Eigen::Matrix3d;


////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////

Residue::Residue(const string &name, int num, const string &ins, Chain *owner):
    __name (name),
    __num  (num),
    __ins  (ins),
    __owner(owner)
{
    if (owner)
    {
        owner->sub().push_back(this);
    }
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __name
////////////////////////////////////////////////////////////////////////////////

string &Residue::name()
{
    return __name;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __num
////////////////////////////////////////////////////////////////////////////////

int Residue::num()
{
    return __num;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __ins
////////////////////////////////////////////////////////////////////////////////

string &Residue::ins()
{
    return __ins;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __owner
////////////////////////////////////////////////////////////////////////////////

Chain *Residue::owner()
{
    return __owner;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __sub
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> &Residue::sub()
{
    return __sub;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __name
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::name(const string &val)
{
    __name = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __num
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::num(int val)
{
    __num = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __ins
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::ins(const string &val)
{
    __ins = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __owner
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::owner(Chain *val)
{
    __owner = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __sub
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::sub(const vector<Atom *> &val)
{
    __sub = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: compNum
////////////////////////////////////////////////////////////////////////////////

string Residue::compNum()
{
    return (format("%d%s") %
        __num              %
        __ins
    ).str();
}


////////////////////////////////////////////////////////////////////////////////
// Setter: compNum
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::compNum(int num, const string &ins)
{
    __num = num;
    __ins = ins;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: compNum (by compNumPair)
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::compNum(const pair<int, string> &compNumPair)
{
    __num = compNumPair.first;
    __ins = compNumPair.second;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::copy()
{
    auto copyResPtr = new Residue(__name, __num, __ins);

    for (auto atomPtr: __sub)
    {
        auto copyAtomPtr = atomPtr->copy();

        copyAtomPtr->owner(copyResPtr);

        copyResPtr->__sub.push_back(copyAtomPtr);
    }

    return copyResPtr;
}


////////////////////////////////////////////////////////////////////////////////
// GetResidues
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> Residue::getResidues()
{
    return {this};
}


////////////////////////////////////////////////////////////////////////////////
// GetAtoms
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Residue::getAtoms()
{
    return __sub;
}


////////////////////////////////////////////////////////////////////////////////
// subMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, Atom *> Residue::subMap()
{
    unordered_map<string, Atom *> atomPtrMap;

    for (auto atomPtr: __sub)
    {
        atomPtrMap.emplace(atomPtr->name(), atomPtr);
    }

    return atomPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// coordMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, RowVector3d> Residue::coordMap()
{
    unordered_map<string, RowVector3d> coordPtrMap;

    for (auto atomPtr: __sub)
    {
        coordPtrMap.emplace(atomPtr->name(), atomPtr->coord());
    }

    return coordPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// Calc Backbone Dihedral Angle
////////////////////////////////////////////////////////////////////////////////

double Residue::calcBBDihedralAngle(DIH dihedralEnum)
{
    auto atomCoordMap = coordMap();

    if (dihedralEnum == DIH::L)
    {
        return calcDihedralAngle(
            prev()->coordMap().at("C"),
            atomCoordMap.at("N"),
            atomCoordMap.at("CA"),
            atomCoordMap.at("C")
        );
    }
    else
    {
        return calcDihedralAngle(
            atomCoordMap.at("N"),
            atomCoordMap.at("CA"),
            atomCoordMap.at("C"),
            next()->coordMap().at("N")
        );
    }
}


////////////////////////////////////////////////////////////////////////////////
// Calc Backbone Rotation Matrix By Delta Angle
////////////////////////////////////////////////////////////////////////////////

pair<RowVector3d, Matrix3d> Residue::calcBBRotationMatrixByDeltaAngle(
    DIH dihedralEnum, SIDE sideEnum, double deltaAngle)
{
    RowVector3d moveCoord;
    Matrix3d rotationMatrix;

    auto atomCoordMap = coordMap();

    if (sideEnum == SIDE::L)
    {
        deltaAngle = -deltaAngle;
    }

    if (dihedralEnum == DIH::L)
    {
        moveCoord = atomCoordMap.at("N");
        rotationMatrix = calcRotationMatrix(atomCoordMap.at("CA") - moveCoord, deltaAngle);
    }
    else
    {
        moveCoord = atomCoordMap.at("CA");
        rotationMatrix = calcRotationMatrix(atomCoordMap.at("C") - moveCoord, deltaAngle);
    }

    return {moveCoord, rotationMatrix};
}


////////////////////////////////////////////////////////////////////////////////
// Calc Backbone Rotation Matrix By Target Angle
////////////////////////////////////////////////////////////////////////////////

pair<RowVector3d, Matrix3d> Residue::calcBBRotationMatrixByTargetAngle(
    DIH dihedralEnum, SIDE sideEnum, double targetAngle)
{
    return calcBBRotationMatrixByDeltaAngle(dihedralEnum, sideEnum,
        targetAngle - calcBBDihedralAngle(dihedralEnum));
}


////////////////////////////////////////////////////////////////////////////////
// Get Backbone Rotation Atom Pointer
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Residue::getBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum)
{
    vector<Atom *> rotationAtomObjList;
    auto iterInOwner = iter();

    if (sideEnum == SIDE::L)
    {
        for (auto resIter = __owner->sub().begin(); resIter != iterInOwner; resIter++)
        {
            for (auto atomPtr: **resIter)
            {
                rotationAtomObjList.push_back(atomPtr);
            }
        }

        if (dihedralEnum == DIH::R)
        {
            for (auto atomPtr: __sub)
            {
                if (atomPtr->name() != "CA" && atomPtr->name() != "C" &&
                    atomPtr->name() != "O" && atomPtr->name() != "OXT")
                {
                    rotationAtomObjList.push_back(atomPtr);
                }
            }
        }
    }
    else
    {
        if (dihedralEnum == DIH::L)
        {
            for (auto atomPtr: __sub)
            {
                if (atomPtr->name() != "N" && atomPtr->name() != "CA")
                {
                    rotationAtomObjList.push_back(atomPtr);
                }
            }
        }
        else
        {
            for (auto atomPtr: __sub)
            {
                if (atomPtr->name() == "O" || atomPtr->name() == "OXT")
                {
                    rotationAtomObjList.push_back(atomPtr);
                }
            }
        }

        for (auto resIter = iterInOwner + 1; resIter != __owner->sub().end(); resIter++)
        {
            for (auto atomPtr: **resIter)
            {
                rotationAtomObjList.push_back(atomPtr);
            }
        }
    }

    return rotationAtomObjList;
}


////////////////////////////////////////////////////////////////////////////////
// Rotate Backbone Dihedral Angle By Delta Angle
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::rotateBBDihedralAngleByDeltaAngle(DIH dihedralEnum,
    SIDE sideEnum, double deltaAngle)
{
    auto [moveCoord, rotationMatrix] = calcBBRotationMatrixByDeltaAngle(
        dihedralEnum, sideEnum, deltaAngle);

    for (auto atomPtr: getBBRotationAtomPtr(dihedralEnum, sideEnum))
    {
        atomPtr->coord((atomPtr->coord() - moveCoord) * rotationMatrix + moveCoord);
    }

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Rotate Backbone Dihedral Angle By Target Angle
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::rotateBBDihedralAngleByTargetAngle(DIH dihedralEnum,
    SIDE sideEnum, double targetAngle)
{
    return rotateBBDihedralAngleByDeltaAngle(dihedralEnum, sideEnum,
        targetAngle - calcBBDihedralAngle(dihedralEnum));
}


////////////////////////////////////////////////////////////////////////////////
// Calc Side Chain Dihedral Angle
////////////////////////////////////////////////////////////////////////////////

double Residue::calcSCDihedralAngle(int dihedralIdx)
{
    auto atomCoordMap = coordMap();

    return calcDihedralAngle(
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx)[0]],
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx)[1]],
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx)[2]],
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx)[3]]);
}


////////////////////////////////////////////////////////////////////////////////
// Calc Side Chain Rotation Matrix By Delta Angle
////////////////////////////////////////////////////////////////////////////////

pair<RowVector3d, Matrix3d> Residue::calcSCRotationMatrixByDeltaAngle(
    int dihedralIdx, double deltaAngle)
{
    auto atomCoordMap = coordMap();

    auto moveCoord = atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(
        __name).at(dihedralIdx)[1]];

    auto rotationMatrix = calcRotationMatrix(atomCoordMap[
        __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx)[2]] -
        moveCoord, deltaAngle);

    return {moveCoord, rotationMatrix};
}


////////////////////////////////////////////////////////////////////////////////
// Calc Side Chain Rotation Matrix By Target Angle
////////////////////////////////////////////////////////////////////////////////

pair<RowVector3d, Matrix3d> Residue::calcSCRotationMatrixByTargetAngle(
    int dihedralIdx, double targetAngle)
{
    return calcSCRotationMatrixByDeltaAngle(dihedralIdx,
        targetAngle - calcSCDihedralAngle(dihedralIdx));
}


////////////////////////////////////////////////////////////////////////////////
// Get Side Chain Rotation Atom Pointer
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Residue::getSCRotationAtomPtr(int dihedralIdx)
{
    vector<Atom *> rotationAtomObjList;

    unordered_set<string> rotationAtomNameSet(
        __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx).cbegin() + 3,
        __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(__name).at(dihedralIdx).cend());

    for (auto atomPtr: __sub)
    {
        if (rotationAtomNameSet.count(atomPtr->name()))
        {
            rotationAtomObjList.push_back(atomPtr);
        }
    }

    return rotationAtomObjList;
}


////////////////////////////////////////////////////////////////////////////////
// Rotate Side Chain Dihedral Angle By Delta Angle
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::rotateSCDihedralAngleByDeltaAngle(int dihedralIdx, double deltaAngle)
{
    auto [moveCoord, rotationMatrix] = calcSCRotationMatrixByDeltaAngle(
        dihedralIdx, deltaAngle);

    for (auto atomPtr: getSCRotationAtomPtr(dihedralIdx))
    {
        atomPtr->coord((atomPtr->coord() - moveCoord) * rotationMatrix + moveCoord);
    }

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Rotate Side Chain Dihedral Angle By Target Angle
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::rotateSCDihedralAngleByTargetAngle(int dihedralIdx, double targetAngle)
{
    return rotateSCDihedralAngleByDeltaAngle(dihedralIdx, targetAngle -
        calcSCDihedralAngle(dihedralIdx));
}


////////////////////////////////////////////////////////////////////////////////
// Dump
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::dump(const string &dumpFilePath, const string &fileMode)
{
    string chainName;

    if (__owner)
    {
        chainName = __owner->name();
    }

    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto atomPtr: __sub)
    {
        if (isdigit(atomPtr->name()[0]) || atomPtr->name().size() == 4)
        {
            fprintf(fo, "ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n",
                atomPtr->num(),
                atomPtr->name().c_str(),
                atomPtr->alt().c_str(),
                __name.c_str(),
                chainName.c_str(),
                __num,
                __ins.c_str(),
                atomPtr->coord()[0],
                atomPtr->coord()[1],
                atomPtr->coord()[2],
                atomPtr->occ().c_str(),
                atomPtr->tempF().c_str(),
                atomPtr->ele().c_str(),
                atomPtr->chg().c_str()
            );
        }
        else
        {
            fprintf(fo, "ATOM  %5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n",
                atomPtr->num(),
                atomPtr->name().c_str(),
                atomPtr->alt().c_str(),
                __name.c_str(),
                chainName.c_str(),
                __num,
                __ins.c_str(),
                atomPtr->coord()[0],
                atomPtr->coord()[1],
                atomPtr->coord()[2],
                atomPtr->occ().c_str(),
                atomPtr->tempF().c_str(),
                atomPtr->ele().c_str(),
                atomPtr->chg().c_str()
            );
        }
    }

    fclose(fo);

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////

Residue::~Residue()
{
    for (auto subPtr: __sub)
    {
        delete subPtr;
    }
}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Residue::__str() const
{
    return (format("<Residue object: %d%s %s, at 0x%p>") %
        __num                                            %
        __ins                                            %
        __name                                           %
        this
    ).str();
}


}  // End namespace PDBTools
