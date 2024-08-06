/*
    Util.hpp
    ========
        Utility functions implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <charconv>
#include <cstdio>
#include <cstdint>
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
using std::from_chars;
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
    static const uint64_t __H_MASK_ARRAY[] = {0xe000000000000ull, 0x100ull, 0ull, 0ull};

    return (__H_MASK_ARRAY[(uint8_t)atomName[0] >> 6] >> ((uint8_t)atomName[0] & 0x3f)) & 0x1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Split CompNum
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

pair<int, string> splitCompNum(const string &compNumStr)
{
    int resNum;

    if (isalpha(compNumStr.back()))
    {
        from_chars(compNumStr.data(), compNumStr.data() + compNumStr.size() - 1, resNum);

        return {resNum, string(1, compNumStr.back())};
    }
    else
    {
        from_chars(compNumStr.data(), compNumStr.data() + compNumStr.size(), resNum);

        return {resNum, ""};
    }
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
