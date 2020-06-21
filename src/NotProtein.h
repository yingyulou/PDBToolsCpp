/*
    NotProtein.h
    ============
        Class __NotProtein header.
*/

#ifndef __PDBTOOLS_NOT_PROTEIN_H
#define __PDBTOOLS_NOT_PROTEIN_H

#include <algorithm>

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::find;


////////////////////////////////////////////////////////////////////////////////
// Class __NotProtein
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename OwnerType>
class __NotProtein
{
public:

    // iter
    typename vector<SelfType *>::iterator iter()
    {
        return find(static_cast<SelfType *>(this)->owner->sub.begin(),
            static_cast<SelfType *>(this)->owner->sub.end(),
            static_cast<SelfType *>(this));
    }


    typename vector<SelfType *>::const_iterator iter() const
    {
        return find(static_cast<const SelfType *>(this)->owner->sub.begin(),
            static_cast<const SelfType *>(this)->owner->sub.end(),
            static_cast<const SelfType *>(this));
    }


    // pre
    SelfType *pre(int shiftLen = 1) { return *(iter() - shiftLen); }
    const SelfType *pre(int shiftLen = 1) const { return *(iter() - shiftLen); }


    // next
    SelfType *next(int shiftLen = 1) { return *(iter() + shiftLen); }
    const SelfType *next(int shiftLen = 1) const { return *(iter() + shiftLen); }


    // Remove
    typename vector<SelfType *>::iterator Remove(bool deteleBool = true)
    {
        auto eraseIter = static_cast<SelfType *>(this)->owner->sub.erase(iter());

        if (deteleBool)
        {
            delete static_cast<SelfType *>(this);
        }

        return eraseIter;
    }
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_NOT_PROTEIN_H
