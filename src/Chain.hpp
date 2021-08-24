/*
    Chain.hpp
    =========
        Class Chain implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>
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

Chain::Chain(const string &name, Protein *owner):
    __name (name),
    __owner(owner)
{
    if (owner)
    {
        owner->sub().push_back(this);
    }
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __name
////////////////////////////////////////////////////////////////////////////////

string &Chain::name()
{
    return __name;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __owner
////////////////////////////////////////////////////////////////////////////////

Protein *Chain::owner()
{
    return __owner;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __sub
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> &Chain::sub()
{
    return __sub;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __name
////////////////////////////////////////////////////////////////////////////////

Chain *Chain::name(const string &val)
{
    __name = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __owner
////////////////////////////////////////////////////////////////////////////////

Chain *Chain::owner(Protein *val)
{
    __owner = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __sub
////////////////////////////////////////////////////////////////////////////////

Chain *Chain::sub(const vector<Residue *> &val)
{
    __sub = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Chain *Chain::copy()
{
    auto copyChainPtr = new Chain(__name);

    for (auto resPtr: __sub)
    {
        auto copyResPtr = resPtr->copy();

        copyResPtr->owner(copyChainPtr);

        copyChainPtr->__sub.push_back(copyResPtr);
    }

    return copyChainPtr;
}


////////////////////////////////////////////////////////////////////////////////
// GetResidues
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> Chain::getResidues()
{
    return __sub;
}


////////////////////////////////////////////////////////////////////////////////
// GetAtoms
////////////////////////////////////////////////////////////////////////////////

vector<Atom *> Chain::getAtoms()
{
    vector<Atom *> atomPtrList;

    for (auto resPtr: __sub)
    {
        for (auto atomPtr: resPtr->sub())
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

    for (auto resPtr: __sub)
    {
        resPtrMap.emplace(resPtr->compNum(), resPtr);
    }

    return resPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// Dump
////////////////////////////////////////////////////////////////////////////////

Chain *Chain::dump(const string &dumpFilePath, const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto resPtr: __sub)
    {
        for (auto atomPtr: *resPtr)
        {
            if (isdigit(atomPtr->name()[0]) || atomPtr->name().size() == 4)
            {
                fprintf(fo, "ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n",
                    atomPtr->num(),
                    atomPtr->name().c_str(),
                    atomPtr->alt().c_str(),
                    resPtr->name().c_str(),
                    __name.c_str(),
                    resPtr->num(),
                    resPtr->ins().c_str(),
                    atomPtr->coord()[0],
                    atomPtr->coord()[1],
                    atomPtr->coord()[2],
                    atomPtr->occ().c_str(),
                    atomPtr->tempF().c_str(),
                    atomPtr->ele().c_str(),
                    atomPtr->chg().c_str()
                );
            }
            else
            {
                fprintf(fo, "ATOM  %5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n",
                    atomPtr->num(),
                    atomPtr->name().c_str(),
                    atomPtr->alt().c_str(),
                    resPtr->name().c_str(),
                    __name.c_str(),
                    resPtr->num(),
                    resPtr->ins().c_str(),
                    atomPtr->coord()[0],
                    atomPtr->coord()[1],
                    atomPtr->coord()[2],
                    atomPtr->occ().c_str(),
                    atomPtr->tempF().c_str(),
                    atomPtr->ele().c_str(),
                    atomPtr->chg().c_str()
                );
            }
        }
    }

    fclose(fo);

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////

Chain::~Chain()
{
    for (auto subPtr: __sub)
    {
        delete subPtr;
    }
}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Chain::__str() const
{
    return (format("<Chain object: %s, at 0x%p>") %
        __name                                    %
        this
    ).str();
}


}  // End namespace PDBTools
