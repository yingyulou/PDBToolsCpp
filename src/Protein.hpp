/*
    Protein.hpp
    ===========
        Class Protein implementation.
*/

#ifndef __PDBTOOLS_PROTEIN_HPP
#define __PDBTOOLS_PROTEIN_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include "Protein.h"
#include "Chain.h"
#include "Residue.h"
#include "Atom.h"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::vector;
using std::unordered_map;


////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////

Protein::Protein(const string &proteinID, int modelNum):
    name(proteinID), model(modelNum) {}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::Copy() const
{
    auto copyProPtr = new Protein(name);

    for (auto chainPtr: sub)
    {
        auto copyChainPtr = chainPtr->Copy();
        copyChainPtr->owner = copyProPtr;
        copyProPtr->sub.push_back(copyChainPtr);
    }

    return copyProPtr;
}


////////////////////////////////////////////////////////////////////////////////
// GetResidues
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> Protein::GetResidues()
{
    vector<Residue *> resPtrList;

    for (auto chainPtr: sub)
    {
        for (auto resPtr: chainPtr->sub)
        {
            resPtrList.push_back(resPtr);
        }
    }

    return resPtrList;
}


vector<const Residue *> Protein::GetResidues() const
{
    vector<const Residue *> resPtrList;

    for (auto chainPtr: sub)
    {
        for (auto resPtr: chainPtr->sub)
        {
            resPtrList.push_back(resPtr);
        }
    }

    return resPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// GetAtoms
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Protein::GetAtoms()
{
    vector<Atom *> atomPtrList;

    for (auto chainPtr: sub)
    {
        for (auto resPtr: chainPtr->sub)
        {
            for (auto atomPtr: resPtr->sub)
            {
                atomPtrList.push_back(atomPtr);
            }
        }
    }

    return atomPtrList;
}


vector<const Atom *> Protein::GetAtoms() const
{
    vector<const Atom *> atomPtrList;

    for (auto chainPtr: sub)
    {
        for (auto resPtr: chainPtr->sub)
        {
            for (auto atomPtr: resPtr->sub)
            {
                atomPtrList.push_back(atomPtr);
            }
        }
    }

    return atomPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// subMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, Chain *> Protein::subMap()
{
    unordered_map<string, Chain *> chainPtrMap;

    for (auto chainPtr: sub) chainPtrMap.emplace(chainPtr->name, chainPtr);

    return chainPtrMap;
}


unordered_map<string, const Chain *> Protein::subMap() const
{
    unordered_map<string, const Chain *> chainPtrMap;

    for (auto chainPtr: sub) chainPtrMap.emplace(chainPtr->name, chainPtr);

    return chainPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////

Protein::~Protein()
{
    for (auto subPtr: sub)
    {
        delete subPtr;
    }
}


}  // End namespace PDBTools


#endif  // __PDBTOOLS_PROTEIN_HPP
