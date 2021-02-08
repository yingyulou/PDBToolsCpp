/*
    StructBase.hpp
    ==============
        Class __StructBase implementation.
*/

#ifndef __PDBTOOLS_STRUCT_BASE_HPP
#define __PDBTOOLS_STRUCT_BASE_HPP

#include <string>
#include <cstdio>
#include "StructBase.h"

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
SelfType *__StructBase<SelfType>::Dump(const string &dumpFilePath,
    const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());
    fprintf(fo, static_cast<SelfType *>(this)->Dumps().c_str());
    fclose(fo);

    return static_cast<SelfType *>(this);
}


}  // End namespace PDBTools


#endif  // __PDBTOOLS_STRUCT_BASE_HPP
