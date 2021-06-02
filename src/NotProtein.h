/*
    NotProtein.h
    ============
        Class __NotProtein header.
*/

#pragma once

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


    // pre
    SelfType *pre(int shiftLen = 1);


    // next
    SelfType *next(int shiftLen = 1);


    // Remove
    typename vector<SelfType *>::iterator Remove(bool deteleBool = true);
};


}  // End namespace PDBTools
