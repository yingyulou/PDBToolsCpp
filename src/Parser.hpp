/*
    Parser.hpp
    ==========
        Parser functions implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <climits>
#include <fstream>
#include <filesystem>
#include <stdexcept>
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
using std::runtime_error;
using boost::algorithm::trim;
using Eigen::RowVector3d;


////////////////////////////////////////////////////////////////////////////////
// Load PDB File
////////////////////////////////////////////////////////////////////////////////

Protein *load(const string &pdbFilePath, bool parseHBool = false)
{
    ifstream f(pdbFilePath);
    string line;

    if (!f)
    {
        throw runtime_error(pdbFilePath + " not exists");
    }

    Protein *proPtr = new Protein(path(pdbFilePath).stem().string());
    Chain *chainPtr;
    Residue *resPtr;

    string lastChainName = " ";
    string lastResName   = " ";
    int lastResNum       = INT_MAX;
    string lastResIns    = " ";

    while (getline(f, line))
    {
        if (line.substr(0, 4) != "ATOM")
        {
            continue;
        }

        string atomName = line.substr(12, 4);

        trim(atomName);

        if (isH(atomName) && !parseHBool)
        {
            continue;
        }

        int atomNum       = stoi(line.substr(6, 5));
        string atomAltLoc = line.substr(16, 1);
        string resName    = line.substr(17, 3);
        string chainName  = line.substr(21, 1);
        int resNum        = stoi(line.substr(22, 4));
        string resIns     = line.substr(26, 1);

        trim(atomAltLoc);
        trim(resName);
        trim(chainName);
        trim(resIns);

        RowVector3d atomCoord(
            stod(line.substr(30, 8)),
            stod(line.substr(38, 8)),
            stod(line.substr(46, 8)));

        string atomOccupancy  = line.size() > 54 ? line.substr(54, 6) : "";
        string atomTempFactor = line.size() > 60 ? line.substr(60, 6) : "";
        string atomElement    = line.size() > 76 ? line.substr(76, 2) : "";
        string atomCharge     = line.size() > 78 ? line.substr(78, 2) : "";

        trim(atomOccupancy);
        trim(atomTempFactor);
        trim(atomElement);
        trim(atomCharge);

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

vector<Protein *> loadModel(const string &pdbFilePath, bool parseHBool = false)
{
    ifstream f(pdbFilePath);
    string line;

    if (!f)
    {
        throw runtime_error(pdbFilePath + " not exists");
    }

    string proName = path(pdbFilePath).stem().string();
    Protein *proPtr = new Protein(proName);
    vector<Protein *> proPtrList {proPtr};
    Chain *chainPtr;
    Residue *resPtr;

    string lastChainName = " ";
    string lastResName   = " ";
    int lastResNum       = INT_MAX;
    string lastResIns    = " ";

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

        string atomName = line.substr(12, 4);

        trim(atomName);

        if (isH(atomName) && !parseHBool)
        {
            continue;
        }

        int atomNum       = stoi(line.substr(6, 5));
        string atomAltLoc = line.substr(16, 1);
        string resName    = line.substr(17, 3);
        string chainName  = line.substr(21, 1);
        int resNum        = stoi(line.substr(22, 4));
        string resIns     = line.substr(26, 1);

        trim(atomAltLoc);
        trim(resName);
        trim(chainName);
        trim(resIns);

        RowVector3d atomCoord(
            stod(line.substr(30, 8)),
            stod(line.substr(38, 8)),
            stod(line.substr(46, 8)));

        string atomOccupancy  = line.size() > 54 ? line.substr(54, 6) : "";
        string atomTempFactor = line.size() > 60 ? line.substr(60, 6) : "";
        string atomElement    = line.size() > 76 ? line.substr(76, 2) : "";
        string atomCharge     = line.size() > 78 ? line.substr(78, 2) : "";

        trim(atomOccupancy);
        trim(atomTempFactor);
        trim(atomElement);
        trim(atomCharge);

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

    if (proPtrList[0]->sub().size() == 0)
    {
        delete proPtrList[0];
        proPtrList.erase(proPtrList.begin());
    }

    return proPtrList;
}


}  // End namespace PDBTools
