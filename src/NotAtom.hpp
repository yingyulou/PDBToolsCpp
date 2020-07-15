/*
    NotAtom.hpp
    ===========
        Class __NotAtom implementation.
*/

#ifndef __PDBTOOLS_NOT_ATOM_HPP
#define __PDBTOOLS_NOT_ATOM_HPP

#include <string>
#include <vector>
#include <unordered_set>
#include <cstdio>
#include <iterator>
#include <initializer_list>
#include <Eigen/Dense>
#include "NotAtom.h"
#include "Struct.h"
#include "Constants.hpp"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::vector;
using std::unordered_set;
using std::distance;
using std::initializer_list;
using Eigen::RowVector3d;


////////////////////////////////////////////////////////////////////////////////
// FilterAtoms
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
vector<Atom *> __NotAtom<SelfType, SubType>::FilterAtoms(
    const unordered_set<string> &atomNameSet)
{
    vector<Atom *> atomPtrList;

    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        if (atomNameSet.count(atomPtr->name))
        {
            atomPtrList.push_back(atomPtr);
        }
    }

    return atomPtrList;
}


template <typename SelfType, typename SubType>
vector<const Atom *> __NotAtom<SelfType, SubType>::FilterAtoms(
    const unordered_set<string> &atomNameSet) const
{
    vector<const Atom *> atomPtrList;

    for (auto atomPtr: static_cast<const SelfType *>(this)->GetAtoms())
    {
        if (atomNameSet.count(atomPtr->name))
        {
            atomPtrList.push_back(atomPtr);
        }
    }

    return atomPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// GetAtomsCoord
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
vector<RowVector3d *> __NotAtom<SelfType, SubType>::GetAtomsCoord()
{
    vector<RowVector3d *> coordPtrList;

    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        coordPtrList.push_back(&atomPtr->coord);
    }

    return coordPtrList;
}


template <typename SelfType, typename SubType>
vector<const RowVector3d *> __NotAtom<SelfType, SubType>::GetAtomsCoord() const
{
    vector<const RowVector3d *> coordPtrList;

    for (auto atomPtr: static_cast<const SelfType *>(this)->GetAtoms())
    {
        coordPtrList.push_back(&atomPtr->coord);
    }

    return coordPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// FilterAtomsCoord
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
vector<RowVector3d *> __NotAtom<SelfType, SubType>::FilterAtomsCoord(
    const unordered_set<string> &atomNameSet)
{
    vector<RowVector3d *> coordPtrList;

    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        if (atomNameSet.count(atomPtr->name))
        {
            coordPtrList.push_back(&atomPtr->coord);
        }
    }

    return coordPtrList;
}


template <typename SelfType, typename SubType>
vector<const RowVector3d *> __NotAtom<SelfType, SubType>::FilterAtomsCoord(
    const unordered_set<string> &atomNameSet) const
{
    vector<const RowVector3d *> coordPtrList;

    for (auto atomPtr: static_cast<const SelfType *>(this)->GetAtoms())
    {
        if (atomNameSet.count(atomPtr->name))
        {
            coordPtrList.push_back(&atomPtr->coord);
        }
    }

    return coordPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// Dumps
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::Dumps() const
{
    string dumpStr;

    for (auto atomPtr: static_cast<const SelfType *>(this)->GetAtoms())
    {
        dumpStr += atomPtr->Dumps();
    }

    return dumpStr;
}


////////////////////////////////////////////////////////////////////////////////
// center
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
RowVector3d __NotAtom<SelfType, SubType>::center() const
{
    auto coordPtrList = GetAtomsCoord();
    double sumX = 0., sumY = 0., sumZ = 0., coordLen = coordPtrList.size();

    for (auto coordPtr: coordPtrList)
    {
        sumX += (*coordPtr)[0];
        sumY += (*coordPtr)[1];
        sumZ += (*coordPtr)[2];
    }

    return RowVector3d(sumX / coordLen, sumY / coordLen, sumZ / coordLen);
}


////////////////////////////////////////////////////////////////////////////////
// MoveCenter
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::MoveCenter()
{
    auto centerCoord = center();

    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        atomPtr->coord -= centerCoord;
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// seq
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::seq() const
{
    string seqStr;

    for (auto resPtr: static_cast<const SelfType *>(this)->GetResidues())
    {
        seqStr += RESIDUE_NAME_THREE_TO_ONE_MAP.at(resPtr->name);
    }

    return seqStr;
}


////////////////////////////////////////////////////////////////////////////////
// DumpFasta
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::DumpFasta(const string &dumpFilePath,
    const string &titleStr, const string &fileMode)
{
    return const_cast<__NotAtom *>(const_cast<const __NotAtom *>(this)->
        DumpFasta(dumpFilePath, fileMode, titleStr));
}


template <typename SelfType, typename SubType>
const SelfType *__NotAtom<SelfType, SubType>::DumpFasta(const string &dumpFilePath,
    const string &titleStr, const string &fileMode) const
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());
    fprintf(fo, "%s", fasta(titleStr).c_str());
    fclose(fo);

    return static_cast<const SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Renum Residues
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::RenumResidues(int startNum)
{
    for (auto resPtr: static_cast<SelfType *>(this)->GetResidues())
    {
        resPtr->compNum(startNum++, "");
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Renum Atoms
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::RenumAtoms(int startNum)
{
    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        atomPtr->num = startNum++;
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// PushBack
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::PushBack(const SubType *subPtr)
{
    SubType *copySubPtr = subPtr->Copy();
    copySubPtr->owner = static_cast<SelfType *>(this);
    static_cast<SelfType *>(this)->sub.push_back(copySubPtr);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Insert
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::Insert(const_iterator insertIter,
    SubType *subPtr)
{
    SubType *copySubPtr = subPtr->Copy();
    copySubPtr->owner = static_cast<SelfType *>(this);
    static_cast<SelfType *>(this)->sub.insert(insertIter, copySubPtr);

    return static_cast<SelfType *>(this);
}


template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::Insert(const_iterator insertIter,
    int n, SubType *subPtr)
{
    auto insertFirstIter = static_cast<SelfType *>(this)->sub.insert(insertIter,
        n, subPtr);

    for (int _ = 0; _ < n; _++)
    {
        *insertFirstIter = (*insertFirstIter)->Copy();
        (*insertFirstIter)->owner = static_cast<SelfType *>(this);
        insertFirstIter++;
    }

    return static_cast<SelfType *>(this);
}


template <typename SelfType, typename SubType>
template <typename ForwardIterator>
SelfType *__NotAtom<SelfType, SubType>::Insert(const_iterator insertIter,
    ForwardIterator firstIter, ForwardIterator lastIter)
{
    auto insertFirstIter = static_cast<SelfType *>(this)->sub.insert(insertIter,
        firstIter, lastIter);

    for (int _ = 0; _ < distance(firstIter, lastIter); _++)
    {
        *insertFirstIter = (*insertFirstIter)->Copy();
        (*insertFirstIter)->owner = static_cast<SelfType *>(this);
        insertFirstIter++;
    }

    return static_cast<SelfType *>(this);
}


template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::Insert(const_iterator insertIter,
    initializer_list<SubType *> initializerList)
{
    auto insertFirstIter = static_cast<SelfType *>(this)->sub.insert(insertIter,
        initializerList);

    for (int _ = 0; _ < initializerList.size(); _++)
    {
        *insertFirstIter = (*insertFirstIter)->Copy();
        (*insertFirstIter)->owner = static_cast<SelfType *>(this);
        insertFirstIter++;
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// MoveBack
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::MoveBack(SubType *subPtr)
{
    static_cast<SelfType *>(this)->sub.push_back(subPtr);

    if (subPtr->owner)
    {
        subPtr->Remove(false);
    }

    subPtr->owner = static_cast<SelfType *>(this);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// MoveInsert
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::MoveInsert(const_iterator insertIter,
    SubType *subPtr)
{
    static_cast<SelfType *>(this)->sub.insert(insertIter, subPtr);

    if (subPtr->owner)
    {
        subPtr->Remove(false);
    }

    subPtr->owner = static_cast<SelfType *>(this);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// RemoveAlt
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::RemoveAlt()
{
    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        if (atomPtr->alt == "A")
        {
            atomPtr->alt = "";
        }
        else if (atomPtr->alt != "")
        {
            atomPtr->Remove();
        }
    }

    return static_cast<SelfType *>(this);
}


}  // End namespace PDBTools


#endif  // __PDBTOOLS_NOT_ATOM_HPP
