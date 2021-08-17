/*
    NotProtein.hpp
    ==============
        Class __NotProtein implementation.
*/

#pragma once

#include <vector>
#include <algorithm>
#include "NotProtein.h"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::vector;
using std::find;


////////////////////////////////////////////////////////////////////////////////
// Iter
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
typename vector<SelfType *>::iterator __NotProtein<SelfType, OwnerType>::iter()
{
    return find(static_cast<SelfType *>(this)->owner()->sub().begin(),
        static_cast<SelfType *>(this)->owner()->sub().end(),
        static_cast<SelfType *>(this));
}


////////////////////////////////////////////////////////////////////////////////
// Prev
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
SelfType *__NotProtein<SelfType, OwnerType>::prev(int shiftLen)
{
    return *(iter() - shiftLen);
}


////////////////////////////////////////////////////////////////////////////////
// Next
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
SelfType *__NotProtein<SelfType, OwnerType>::next(int shiftLen)
{
    return *(iter() + shiftLen);
}


////////////////////////////////////////////////////////////////////////////////
// Remove
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
typename vector<SelfType *>::iterator __NotProtein<SelfType, OwnerType>::remove(
    bool deteleBool)
{
    auto eraseIter = static_cast<SelfType *>(this)->owner()->sub().erase(iter());

    if (deteleBool)
    {
        delete static_cast<SelfType *>(this);
    }

    return eraseIter;
}


}  // End namespace PDBTools
