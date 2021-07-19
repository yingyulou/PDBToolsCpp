/*
    StructBase.hpp
    ==============
        Class __StructBase implementation.
*/

#pragma once

#include <string>
#include <cstdio>
#include "StructBase.h"
#include "Predeclaration.h"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;


////////////////////////////////////////////////////////////////////////////////
// Dump
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
SelfType *__StructBase<SelfType>::dump(const string &dumpFilePath,
    const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto atomPtr: static_cast<SelfType *>(this)->getAtoms())
    {
        fprintf(fo, atomPtr->dumps().c_str());
    }

    fclose(fo);

    return static_cast<SelfType *>(this);
}


template <>
Atom *__StructBase<Atom>::dump(const string &dumpFilePath,
    const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    fprintf(fo, static_cast<Atom *>(this)->dumps().c_str());

    fclose(fo);

    return static_cast<Atom *>(this);
}


}  // End namespace PDBTools
