/*
    StructBase.h
    ============
        Class __StructBase header.
*/

#ifndef __PDBTOOLS_STRUCT_BASE_H
#define __PDBTOOLS_STRUCT_BASE_H

#include <string>

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;


////////////////////////////////////////////////////////////////////////////////
// Class __StructBase
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType>
class __StructBase
{
public:

    // Dump
    SelfType *Dump(const string &dumpFilePath, const string &fileMode = "w");
    const SelfType *Dump(const string &dumpFilePath, const string &fileMode = "w") const;
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_STRUCT_BASE_H
