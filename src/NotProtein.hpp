/*
    NotProtein.hpp
    ==============
        Class __NotProtein implementation.
*/

#ifndef __PDBTOOLS_NOT_PROTEIN_HPP
#define __PDBTOOLS_NOT_PROTEIN_HPP

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
// iter
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
typename vector<SelfType *>::iterator __NotProtein<SelfType, OwnerType>::iter()
{
    return find(static_cast<SelfType *>(this)->owner->sub.begin(),
        static_cast<SelfType *>(this)->owner->sub.end(),
        static_cast<SelfType *>(this));
}


template <typename SelfType, typename OwnerType>
typename vector<SelfType *>::const_iterator __NotProtein<SelfType, OwnerType>::iter() const
{
    return find(static_cast<const SelfType *>(this)->owner->sub.begin(),
        static_cast<const SelfType *>(this)->owner->sub.end(),
        static_cast<const SelfType *>(this));
}


////////////////////////////////////////////////////////////////////////////////
// pre
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
SelfType *__NotProtein<SelfType, OwnerType>::pre(int shiftLen)
{
    return *(iter() - shiftLen);
}


template <typename SelfType, typename OwnerType>
const SelfType *__NotProtein<SelfType, OwnerType>::pre(int shiftLen) const
{
    return *(iter() - shiftLen);
}


////////////////////////////////////////////////////////////////////////////////
// next
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
SelfType *__NotProtein<SelfType, OwnerType>::next(int shiftLen)
{
    return *(iter() + shiftLen);
}


template <typename SelfType, typename OwnerType>
const SelfType *__NotProtein<SelfType, OwnerType>::next(int shiftLen) const
{
    return *(iter() + shiftLen);
}


////////////////////////////////////////////////////////////////////////////////
// Remove
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
typename vector<SelfType *>::iterator __NotProtein<SelfType, OwnerType>::Remove(
    bool deteleBool)
{
    auto eraseIter = static_cast<SelfType *>(this)->owner->sub.erase(iter());

    if (deteleBool)
    {
        delete static_cast<SelfType *>(this);
    }

    return eraseIter;
}


}  // End namespace PDBTools


#endif  // __PDBTOOLS_NOT_PROTEIN_HPP
