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

    // Iter
    typename vector<SelfType *>::iterator iter();


    // Prev
    SelfType *prev(int shiftLen = 1);


    // Next
    SelfType *next(int shiftLen = 1);


    // Remove
    typename vector<SelfType *>::iterator remove(bool deteleBool = true);
};


}  // End namespace PDBTools
