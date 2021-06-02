/*
    Chain.hpp
    =========
        Class Chain implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <boost/format.hpp>
#include "Chain.h"
#include "Protein.h"
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

Chain::Chain(const string &chainName, Protein *chainOwner):
    name(chainName), owner(chainOwner)
{
    if (chainOwner)
    {
        chainOwner->sub.push_back(this);
    }
}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Chain::str() const
{
    return (format("<Chain object: %s, at 0x%p>") % name % this).str();
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Chain *Chain::Copy()
{
    auto copyChainPtr = new Chain(name);

    for (auto resPtr: sub)
    {
        auto copyResPtr = resPtr->Copy();
        copyResPtr->owner = copyChainPtr;
        copyChainPtr->sub.push_back(copyResPtr);
    }

    return copyChainPtr;
}


////////////////////////////////////////////////////////////////////////////////
// GetResidues
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> Chain::GetResidues()
{
    return sub;
}


////////////////////////////////////////////////////////////////////////////////
// GetAtoms
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Chain::GetAtoms()
{
    vector<Atom *> atomPtrList;

    for (auto resPtr: sub)
    {
        for (auto atomPtr: resPtr->sub)
        {
            atomPtrList.push_back(atomPtr);
        }
    }

    return atomPtrList;
}


////////////////////////////////////////////////////////////////////////////////
// subMap
////////////////////////////////////////////////////////////////////////////////

unordered_map<string, Residue *> Chain::subMap()
{
    unordered_map<string, Residue *> resPtrMap;

    for (auto resPtr: sub)
    {
        resPtrMap.emplace(resPtr->compNum(), resPtr);
    }

    return resPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////

Chain::~Chain()
{
    for (auto subPtr: sub)
    {
        delete subPtr;
    }
}


}  // End namespace PDBTools
