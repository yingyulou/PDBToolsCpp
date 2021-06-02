/*
    BaseStruct.hpp
    ==============
        Class __BaseStruct implementation.
*/

#pragma once

#include <string>
#include <cstdio>
#include "BaseStruct.h"
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
SelfType *__BaseStruct<SelfType>::Dump(const string &dumpFilePath,
    const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        fprintf(fo, atomPtr->Dumps().c_str());
    }

    fclose(fo);

    return static_cast<SelfType *>(this);
}


template <>
Atom *__BaseStruct<Atom>::Dump(const string &dumpFilePath,
    const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    fprintf(fo, static_cast<Atom *>(this)->Dumps().c_str());

    fclose(fo);

    return static_cast<Atom *>(this);
}


}  // End namespace PDBTools
