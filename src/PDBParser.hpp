/*
    PDBParser.hpp
    =============
        PDB parser functions implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <climits>
#include <fstream>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <Eigen/Dense>
#include "Protein.h"
#include "Chain.h"
#include "Residue.h"
#include "Atom.h"
#include "Util.hpp"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::vector;
using std::ifstream;
using std::getline;
using std::stoi;
using std::stod;
using std::filesystem::path;
using boost::algorithm::trim_copy;
using Eigen::RowVector3d;


////////////////////////////////////////////////////////////////////////////////
// Load PDB File
////////////////////////////////////////////////////////////////////////////////

Protein *Load(const string &pdbFilePath, bool parseHBool = false)
{
    Protein *proPtr = new Protein(path(pdbFilePath).stem().string());
    Chain *chainPtr;
    Residue *resPtr;

    string lastChainName = " ";
    string lastResName   = " ";
    int lastResNum       = INT_MAX;
    string lastResIns    = " ";

    string atomName;
    int atomNum;
    string atomAltLoc;
    string resName;
    string chainName;
    int resNum;
    string resIns;
    RowVector3d atomCoord;
    string atomOccupancy;
    string atomTempFactor;
    string atomElement;
    string atomCharge;

    ifstream f(pdbFilePath);
    string line;

    while (getline(f, line))
    {
        if (line.compare(0, 4, "ATOM") != 0)
        {
            continue;
        }

        atomName = trim_copy(line.substr(12, 4));

        if (IsH(atomName) && !parseHBool)
        {
            continue;
        }

        atomNum    = stoi(line.substr(6, 5));
        atomAltLoc = trim_copy(line.substr(16, 1));
        resName    = trim_copy(line.substr(17, 3));
        chainName  = trim_copy(line.substr(21, 1));
        resNum     = stoi(line.substr(22, 4));
        resIns     = trim_copy(line.substr(26, 1));

        atomCoord << stod(line.substr(30, 8)),
                     stod(line.substr(38, 8)),
                     stod(line.substr(46, 8));

        atomOccupancy  = trim_copy(line.substr(54, 6));
        atomTempFactor = trim_copy(line.substr(60, 6));
        atomElement    = trim_copy(line.substr(76, 2));
        atomCharge     = trim_copy(line.substr(78, 2));

        if (chainName != lastChainName)
        {
            lastChainName = chainName;
            lastResNum    = resNum;
            lastResName   = resName;
            lastResIns    = resIns;
            chainPtr      = new Chain(chainName, proPtr);
            resPtr        = new Residue(resName, resNum, resIns, chainPtr);
        }
        else if (lastResNum != resNum || lastResName != resName || lastResIns != resIns)
        {
            lastResNum  = resNum;
            lastResName = resName;
            lastResIns  = resIns;
            resPtr      = new Residue(resName, resNum, resIns, chainPtr);
        }

        new Atom(atomName, atomNum, atomCoord, atomAltLoc, atomOccupancy,
            atomTempFactor, atomElement, atomCharge, resPtr);
    }

    f.close();

    return proPtr;
}


////////////////////////////////////////////////////////////////////////////////
// Load PDB File With Model
////////////////////////////////////////////////////////////////////////////////

vector<Protein *> LoadModel(const string &pdbFilePath, bool parseHBool = false)
{
    string proName = path(pdbFilePath).stem().string();
    Protein *proPtr = new Protein(proName);
    vector<Protein *> proPtrList {proPtr};
    Chain *chainPtr;
    Residue *resPtr;

    string lastChainName = " ";
    string lastResName   = " ";
    int lastResNum       = INT_MAX;
    string lastResIns    = " ";

    string atomName;
    int atomNum;
    string atomAltLoc;
    string resName;
    string chainName;
    int resNum;
    string resIns;
    RowVector3d atomCoord;
    string atomOccupancy;
    string atomTempFactor;
    string atomElement;
    string atomCharge;

    ifstream f(pdbFilePath);
    string line;

    while (getline(f, line))
    {
        if (line.compare(0, 5, "MODEL") == 0)
        {
            proPtr = new Protein(proName, stoi(line.substr(10, 4)));
            proPtrList.push_back(proPtr);

            lastChainName = " ";
            lastResName   = " ";
            lastResNum    = INT_MAX;
            lastResIns    = " ";

            continue;
        }
        else if (line.compare(0, 4, "ATOM") != 0)
        {
            continue;
        }

        atomName = trim_copy(line.substr(12, 4));

        if (IsH(atomName) && !parseHBool)
        {
            continue;
        }

        atomNum    = stoi(line.substr(6, 5));
        atomAltLoc = trim_copy(line.substr(16, 1));
        resName    = trim_copy(line.substr(17, 3));
        chainName  = trim_copy(line.substr(21, 1));
        resNum     = stoi(line.substr(22, 4));
        resIns     = trim_copy(line.substr(26, 1));

        atomCoord << stod(line.substr(30, 8)),
                     stod(line.substr(38, 8)),
                     stod(line.substr(46, 8));

        atomOccupancy  = trim_copy(line.substr(54, 6));
        atomTempFactor = trim_copy(line.substr(60, 6));
        atomElement    = trim_copy(line.substr(76, 2));
        atomCharge     = trim_copy(line.substr(78, 2));

        if (chainName != lastChainName)
        {
            lastChainName = chainName;
            lastResNum    = resNum;
            lastResName   = resName;
            lastResIns    = resIns;
            chainPtr      = new Chain(chainName, proPtr);
            resPtr        = new Residue(resName, resNum, resIns, chainPtr);
        }
        else if (lastResNum != resNum || lastResName != resName || lastResIns != resIns)
        {
            lastResNum  = resNum;
            lastResName = resName;
            lastResIns  = resIns;
            resPtr      = new Residue(resName, resNum, resIns, chainPtr);
        }

        new Atom(atomName, atomNum, atomCoord, atomAltLoc, atomOccupancy,
            atomTempFactor, atomElement, atomCharge, resPtr);
    }

    f.close();

    if (proPtrList[0]->sub.size() == 0)
    {
        delete proPtrList[0];
        proPtrList.erase(proPtrList.begin());
    }

    return proPtrList;
}


}  // End namespace PDBTools
