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

Residue::Residue(const string &resName, int resNum, const string &resIns,
    Chain *resOwner): name(resName), num(resNum), ins(resIns), owner(resOwner)
{
    if (resOwner)
    {
        resOwner->sub.push_back(this);
    }
}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Residue::str() const
{
    return (format("<Residue object: %d%s %s, at 0x%p>") % num % ins % name % this).str();
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Residue *Residue::copy()
{
    auto copyResPtr = new Residue(name, num, ins);

    for (auto atomPtr: sub)
    {
        auto copyAtomPtr = atomPtr->copy();
        copyAtomPtr->owner = copyResPtr;
        copyResPtr->sub.push_back(copyAtomPtr);
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
    return sub;
}


////////////////////////////////////////////////////////////////////////////////
// compNum
////////////////////////////////////////////////////////////////////////////////

string Residue::compNum()
{
    return (format("%d%s") % num % ins).str();
}


Residue *Residue::compNum(int resNum, const string &resIns)
{
    num = resNum;
    ins = resIns;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// subMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, Atom *> Residue::subMap()
{
    unordered_map<string, Atom *> atomPtrMap;

    for (auto atomPtr: sub)
    {
        atomPtrMap.emplace(atomPtr->name, atomPtr);
    }

    return atomPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// coordMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, RowVector3d> Residue::coordMap()
{
    unordered_map<string, RowVector3d> coordPtrMap;

    for (auto atomPtr: sub)
    {
        coordPtrMap.emplace(atomPtr->name, atomPtr->coord);
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
        return calcDihedralAngle(pre()->coordMap().at("C"), atomCoordMap.at("N"),
            atomCoordMap.at("CA"), atomCoordMap.at("C"));
    }
    else
    {
        return calcDihedralAngle(atomCoordMap.at("N"), atomCoordMap.at("CA"),
            atomCoordMap.at("C"), next()->coordMap().at("N"));
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
        for (auto i = owner->sub.begin(); i != iterInOwner; i++)
        {
            for (auto atomPtr: (*i)->sub)
            {
                rotationAtomObjList.push_back(atomPtr);
            }
        }

        if (dihedralEnum == DIH::R)
        {
            for (auto atomPtr: sub)
            {
                if (atomPtr->name != "CA" && atomPtr->name != "C" &&
                    atomPtr->name != "O" && atomPtr->name != "OXT")
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
            for (auto atomPtr: sub)
            {
                if (atomPtr->name != "N" && atomPtr->name != "CA")
                {
                    rotationAtomObjList.push_back(atomPtr);
                }
            }
        }
        else
        {
            for (auto atomPtr: sub)
            {
                if (atomPtr->name == "O" || atomPtr->name == "OXT")
                {
                    rotationAtomObjList.push_back(atomPtr);
                }
            }
        }

        for (auto i = iterInOwner + 1; i != owner->sub.end(); i++)
        {
            for (auto atomPtr: (*i)->sub)
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
        atomPtr->coord = (atomPtr->coord - moveCoord) * rotationMatrix + moveCoord;
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
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx)[0]],
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx)[1]],
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx)[2]],
        atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx)[3]]);
}


////////////////////////////////////////////////////////////////////////////////
// Calc Side Chain Rotation Matrix By Delta Angle
////////////////////////////////////////////////////////////////////////////////

pair<RowVector3d, Matrix3d> Residue::calcSCRotationMatrixByDeltaAngle(
    int dihedralIdx, double deltaAngle)
{
    auto atomCoordMap = coordMap();

    auto moveCoord = atomCoordMap[__RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(
        name).at(dihedralIdx)[1]];

    auto rotationMatrix = calcRotationMatrix(atomCoordMap[
        __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx)[2]] -
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
        __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx).cbegin() + 3,
        __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP.at(name).at(dihedralIdx).cend());

    for (auto atomPtr: sub)
    {
        if (rotationAtomNameSet.count(atomPtr->name))
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
        atomPtr->coord = (atomPtr->coord - moveCoord) * rotationMatrix + moveCoord;
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
// Destructor
////////////////////////////////////////////////////////////////////////////////

Residue::~Residue()
{
    for (auto subPtr: sub)
    {
        delete subPtr;
    }
}


}  // End namespace PDBTools
