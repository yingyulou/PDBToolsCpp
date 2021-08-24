/*
    NotAtom.hpp
    ===========
        Class __NotAtom implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_set>
#include <cstdio>
#include <iterator>
#include <initializer_list>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "NotAtom.h"
#include "Predecl.h"
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
using boost::format;
using Eigen::RowVector3d;
using Eigen::Matrix;
using Eigen::Dynamic;


////////////////////////////////////////////////////////////////////////////////
// Begin
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
typename vector<SubType *>::iterator __NotAtom<SelfType, SubType>::begin()
{
    return static_cast<SelfType *>(this)->sub().begin();
}


////////////////////////////////////////////////////////////////////////////////
// End
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
typename vector<SubType *>::iterator __NotAtom<SelfType, SubType>::end()
{
    return static_cast<SelfType *>(this)->sub().end();
}


////////////////////////////////////////////////////////////////////////////////
// Filter Atoms
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
vector<Atom *> __NotAtom<SelfType, SubType>::filterAtoms(
    const unordered_set<string> &atomNameSet)
{
    vector<Atom *> atomPtrList;

    for (auto atomPtr: static_cast<SelfType *>(this)->getAtoms())
    {
        if (atomNameSet.count(atomPtr->name()))
        {
            atomPtrList.push_back(atomPtr);
        }
    }

    return atomPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// Get Atoms Coord
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
Matrix<double, Dynamic, 3> __NotAtom<SelfType, SubType>::getAtomsCoord()
{
    auto atomPtrList = static_cast<SelfType *>(this)->getAtoms();

    Matrix<double, Dynamic, 3> coordMatrix(atomPtrList.size(), 3);

    for (int idx = 0; idx < atomPtrList.size(); idx++)
    {
        coordMatrix.row(idx) = atomPtrList[idx]->coord();
    }

    return coordMatrix;
}


////////////////////////////////////////////////////////////////////////////////
// Filter Atoms Coord
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
Matrix<double, Dynamic, 3> __NotAtom<SelfType, SubType>::filterAtomsCoord(
    const unordered_set<string> &atomNameSet)
{
    auto atomPtrList = static_cast<SelfType *>(this)->filterAtoms();

    Matrix<double, Dynamic, 3> coordMatrix(atomPtrList.size(), 3);

    for (int idx = 0; idx < atomPtrList.size(); idx++)
    {
        if (atomNameSet.count(atomPtrList[idx]->name()))
        {
            coordMatrix.row(idx) = atomPtrList[idx]->coord();
        }
    }

    return coordMatrix;
}


////////////////////////////////////////////////////////////////////////////////
// Center
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
RowVector3d __NotAtom<SelfType, SubType>::center()
{
    return getAtomsCoord().colwise().mean();
}


////////////////////////////////////////////////////////////////////////////////
// Move Center
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::moveCenter()
{
    auto centerCoord = center();

    for (auto atomPtr: static_cast<SelfType *>(this)->getAtoms())
    {
        atomPtr->coord() -= centerCoord;
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Seq
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::seq()
{
    string seqStr;

    for (auto resPtr: static_cast<SelfType *>(this)->getResidues())
    {
        seqStr += RESIDUE_NAME_THREE_TO_ONE_MAP.at(resPtr->name());
    }

    return seqStr;
}


////////////////////////////////////////////////////////////////////////////////
// Fasta Str
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::fastaStr(const string &titleStr)
{
    return (format(">%s\n%s\n")                                               %
        (titleStr.empty() ? static_cast<SelfType *>(this)->name() : titleStr) %
        seq()
    ).str();
}


////////////////////////////////////////////////////////////////////////////////
// Dump Fasta
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::dumpFasta(const string &dumpFilePath,
    const string &fileMode, const string &titleStr)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    fprintf(fo, "%s", fastaStr(titleStr).c_str());

    fclose(fo);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Renum Residues
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::renumResidues(int startNum)
{
    for (auto resPtr: static_cast<SelfType *>(this)->getResidues())
    {
        resPtr->compNum(startNum++, "");
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Renum Atoms
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::renumAtoms(int startNum)
{
    for (auto atomPtr: static_cast<SelfType *>(this)->getAtoms())
    {
        atomPtr->num(startNum++);
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Append
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::append(SubType *subPtr, bool copyBool)
{
    if (copyBool)
    {
        subPtr = subPtr->copy();
    }

    subPtr->owner(static_cast<SelfType *>(this));

    static_cast<SelfType *>(this)->sub().push_back(subPtr);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Insert
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::insert(
    typename vector<SubType *>::iterator insertIter, SubType *subPtr, bool copyBool)
{
    if (copyBool)
    {
        subPtr = subPtr->copy();
    }

    subPtr->owner(static_cast<SelfType *>(this));

    static_cast<SelfType *>(this)->sub().insert(insertIter, subPtr);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Remove Alt
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::removeAlt()
{
    for (auto atomPtr: static_cast<SelfType *>(this)->getAtoms())
    {
        if (atomPtr->alt() == "A")
        {
            atomPtr->alt("");
        }
        else if (atomPtr->alt() != "")
        {
            atomPtr->remove();
        }
    }

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Dump Str
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::dumpStr()
{
    string pdbStr;

    for (auto atomPtr: static_cast<SelfType *>(this)->getAtoms())
    {
        pdbStr += atomPtr->dumpStr();
    }

    return pdbStr;
}


}  // End namespace PDBTools
