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
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "NotAtom.h"
#include "Predeclaration.h"
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
// Iterator
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
typename vector<SubType *>::iterator __NotAtom<SelfType, SubType>::begin()
{
    return static_cast<SelfType *>(this)->sub.begin();
}


template <typename SelfType, typename SubType>
typename vector<SubType *>::iterator __NotAtom<SelfType, SubType>::end()
{
    return static_cast<SelfType *>(this)->sub.end();
}


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


////////////////////////////////////////////////////////////////////////////////
// GetAtomsCoord
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
Matrix<double, Dynamic, 3> __NotAtom<SelfType, SubType>::GetAtomsCoord()
{
    auto atomPtrList = static_cast<SelfType *>(this)->GetAtoms();
    Matrix<double, Dynamic, 3> coordMatrix(atomPtrList.size(), 3);

    for (int idx = 0; idx < atomPtrList.size(); idx++)
    {
        coordMatrix.row(idx) = atomPtrList[idx]->coord;
    }

    return coordMatrix;
}


////////////////////////////////////////////////////////////////////////////////
// FilterAtomsCoord
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
Matrix<double, Dynamic, 3> __NotAtom<SelfType, SubType>::FilterAtomsCoord(
    const unordered_set<string> &atomNameSet)
{
    auto atomPtrList = static_cast<SelfType *>(this)->GetAtoms();
    Matrix<double, Dynamic, 3> coordMatrix(atomPtrList.size(), 3);

    for (int idx = 0; idx < atomPtrList.size(); idx++)
    {
        if (atomNameSet.count(atomPtrList[idx]->name))
        {
            coordMatrix.row(idx) = atomPtrList[idx]->coord;
        }
    }

    return coordMatrix;
}


////////////////////////////////////////////////////////////////////////////////
// Dumps
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::Dumps()
{
    string dumpStr;

    for (auto atomPtr: static_cast<SelfType *>(this)->GetAtoms())
    {
        dumpStr += atomPtr->Dumps();
    }

    return dumpStr;
}


////////////////////////////////////////////////////////////////////////////////
// center
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
RowVector3d __NotAtom<SelfType, SubType>::center()
{
    return GetAtomsCoord().colwise().mean();
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
string __NotAtom<SelfType, SubType>::seq()
{
    string seqStr;

    for (auto resPtr: static_cast<SelfType *>(this)->GetResidues())
    {
        seqStr += RESIDUE_NAME_THREE_TO_ONE_MAP.at(resPtr->name);
    }

    return seqStr;
}


////////////////////////////////////////////////////////////////////////////////
// fasta
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
string __NotAtom<SelfType, SubType>::fasta(const string &titleStr)
{
    return (format(">%s\n%s\n") %
        (titleStr.empty() ? static_cast<SelfType *>(this)->name : titleStr) %
        seq()).str();
}


////////////////////////////////////////////////////////////////////////////////
// DumpFasta
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::DumpFasta(const string &dumpFilePath,
    const string &titleStr, const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());
    fprintf(fo, "%s", fasta(titleStr).c_str());
    fclose(fo);

    return static_cast<SelfType *>(this);
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
// Append
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::Append(SubType *subPtr, bool copyBool)
{
    if (copyBool)
    {
        subPtr = subPtr->Copy();
    }

    subPtr->owner = static_cast<SelfType *>(this);
    static_cast<SelfType *>(this)->sub.push_back(subPtr);

    return static_cast<SelfType *>(this);
}


////////////////////////////////////////////////////////////////////////////////
// Insert
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
SelfType *__NotAtom<SelfType, SubType>::Insert(
    typename vector<SubType *>::iterator insertIter, SubType *subPtr, bool copyBool)
{
    if (copyBool)
    {
        subPtr = subPtr->Copy();
    }

    subPtr->owner = static_cast<SelfType *>(this);
    static_cast<SelfType *>(this)->sub.insert(insertIter, subPtr);

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
