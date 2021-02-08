/*
    NotProtein.h
    ============
        Class __NotProtein header.
*/

#ifndef __PDBTOOLS_NOT_PROTEIN_H
#define __PDBTOOLS_NOT_PROTEIN_H

#include <vector>

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::vector;


////////////////////////////////////////////////////////////////////////////////
// Class __NotProtein
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
class __NotProtein
{
public:

    // iter
    typename vector<SelfType *>::iterator iter();
    typename vector<SelfType *>::const_iterator iter() const;


    // pre
    SelfType *pre(int shiftLen = 1);
    const SelfType *pre(int shiftLen = 1) const;


    // next
    SelfType *next(int shiftLen = 1);
    const SelfType *next(int shiftLen = 1) const;


    // Remove
    typename vector<SelfType *>::iterator Remove(bool deteleBool = true);
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_NOT_PROTEIN_H
