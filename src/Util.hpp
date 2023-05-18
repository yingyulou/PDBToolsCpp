/*
    Util.hpp
    ========
        Utility functions implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
#include <regex>
#include <utility>
#include "Protein.h"
#include "Chain.h"
#include "Residue.h"
#include "Atom.h"
#include "Constants.hpp"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::stoi;
using std::vector;
using std::ostream;
using std::regex_search;
using std::regex_match;
using std::smatch;
using std::pair;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// operator<< (Protein)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream &operator<<(ostream &os, const Protein &proObj)
{
    return os << proObj.__str();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// operator<< (Chain)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream &operator<<(ostream &os, const Chain &chainObj)
{
    return os << chainObj.__str();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// operator<< (Residue)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream &operator<<(ostream &os, const Residue &resObj)
{
    return os << resObj.__str();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// operator<< (Atom)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream &operator<<(ostream &os, const Atom &atomObj)
{
    return os << atomObj.__str();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Is H
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isH(const string &atomName)
{
    return regex_search(atomName, __H_RE);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Split CompNum
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

pair<int, string> splitCompNum(const string &compNumStr)
{
    smatch smatchObj;

    regex_match(compNumStr, smatchObj, __COMP_NUM_RE);

    return {stoi(smatchObj[1]), smatchObj[2]};
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dump Str
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
string dumpStr(const T &structPtrList)
{
    string pdbStr;

    for (auto structPtr: structPtrList)
    {
        pdbStr += structPtr->dumpStr();
    }

    return pdbStr;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dump
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void dump(const T &structPtrList, const string &dumpFilePath, const string &fileMode = "w")
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto structPtr: structPtrList)
    {
        fprintf(fo, "%s", structPtr->dumpStr().c_str());
    }

    fclose(fo);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fasta Str
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
string fastaStr(const T &structPtrList)
{
    string fastaStr;

    for (auto structPtr: structPtrList)
    {
        fastaStr += structPtr->fastaStr();
    }

    return fastaStr;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dump Fasta
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void dumpFasta(const T &structPtrList, const string &dumpFilePath, const string &fileMode = "w")
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto structPtr: structPtrList)
    {
        fprintf(fo, "%s", structPtr->fastaStr().c_str());
    }

    fclose(fo);
}


}  // End namespace PDBTools
