/*
    Protein.hpp
    ===========
        Class Protein implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <boost/format.hpp>
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
using boost::format;


////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////

Protein::Protein(const string &proteinID, int modelNum):
    name (proteinID),
    model(modelNum) {}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Protein::str() const
{
    return (format("<Protein object: %s (Model: %d), at 0x%p>") %
        name                                                    %
        model                                                   %
        this
    ).str();
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::copy()
{
    auto copyProPtr = new Protein(name);

    for (auto chainPtr: sub)
    {
        auto copyChainPtr = chainPtr->copy();
        copyChainPtr->owner = copyProPtr;
        copyProPtr->sub.push_back(copyChainPtr);
    }

    return copyProPtr;
}


////////////////////////////////////////////////////////////////////////////////
// GetResidues
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> Protein::getResidues()
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


////////////////////////////////////////////////////////////////////////////////
// GetAtoms
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Protein::getAtoms()
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


////////////////////////////////////////////////////////////////////////////////
// subMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, Chain *> Protein::subMap()
{
    unordered_map<string, Chain *> chainPtrMap;

    for (auto chainPtr: sub)
    {
        chainPtrMap.emplace(chainPtr->name, chainPtr);
    }

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
