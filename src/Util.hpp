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
#include "StructBase.h"
#include "Constants.hpp"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::stoi;
using std::vector;
using std::ostream;
using std::regex_search;
using std::regex_match;
using std::smatch;
using std::pair;


////////////////////////////////////////////////////////////////////////////////
// operator<<
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
ostream &operator<<(ostream &os, const __StructBase<SelfType> &structObj)
{
    return os << static_cast<const SelfType &>(structObj).str();
}


////////////////////////////////////////////////////////////////////////////////
// If An Atom Name is H
////////////////////////////////////////////////////////////////////////////////

bool IsH(const string &atomName)
{
    return regex_search(atomName, __H_RE);
}


////////////////////////////////////////////////////////////////////////////////
// Split CompNum To ResNum + ResIns
////////////////////////////////////////////////////////////////////////////////

pair<int, string> SplitCompNum(const string &compNumStr)
{
    smatch smatchObj;

    regex_match(compNumStr, smatchObj, __COMP_NUM_RE);

    return {stoi(smatchObj[1]), smatchObj[2]};
}


////////////////////////////////////////////////////////////////////////////////
// Get Dump String Of Struct Object List
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
string Dumpls(const vector<SelfType *> &structPtrList)
{
    string dumpStr;

    for (auto structPtr: structPtrList)
    {
        dumpStr += structPtr->Dumps();
    }

    return dumpStr;
}


////////////////////////////////////////////////////////////////////////////////
// Dump Struct Object List To PDB File
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
void Dumpl(const vector<SelfType *> &structPtrList, const string &dumpFilePath,
    const string &fileMode = "w")
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto structPtr: structPtrList)
    {
        fprintf(fo, "%s", structPtr->Dumps().c_str());
    }

    fclose(fo);
}


////////////////////////////////////////////////////////////////////////////////
// Get Fasta String Of Struct Object List
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
string DumpFastals(const vector<SelfType *> &structPtrList)
{
    string dumpStr;

    for (auto structPtr: structPtrList)
    {
        dumpStr += structPtr->fasta();
    }

    return dumpStr;
}


////////////////////////////////////////////////////////////////////////////////
// Dump Struct Object List To Fasta File
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
void DumpFastal(const vector<SelfType *> &structPtrList, const string &dumpFilePath,
    const string &fileMode = "w")
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto structPtr: structPtrList)
    {
        fprintf(fo, "%s", structPtr->fasta().c_str());
    }

    fclose(fo);
}


}  // End namespace PDBTools
