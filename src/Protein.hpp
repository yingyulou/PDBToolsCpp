/*
    Protein.hpp
    ===========
        Class Protein implementation.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>
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

Protein::Protein(const string &name, int model):
    __name (name),
    __model(model) {}


////////////////////////////////////////////////////////////////////////////////
// Getter: __name
////////////////////////////////////////////////////////////////////////////////

string &Protein::name()
{
    return __name;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __model
////////////////////////////////////////////////////////////////////////////////

int Protein::model()
{
    return __model;
}


////////////////////////////////////////////////////////////////////////////////
// Getter: __sub
////////////////////////////////////////////////////////////////////////////////

vector<Chain *> &Protein::sub()
{
    return __sub;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __name
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::name(const string &val)
{
    __name = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __model
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::model(int val)
{
    __model = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Setter: __sub
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::sub(const vector<Chain *> &val)
{
    __sub = val;

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Protein::str() const
{
    return (format("<Protein object: %s (Model: %d), at 0x%p>") %
        __name                                                  %
        __model                                                 %
        this
    ).str();
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::copy()
{
    auto copyProPtr = new Protein(__name);

    for (auto chainPtr: __sub)
    {
        auto copyChainPtr = chainPtr->copy();

        copyChainPtr->owner(copyProPtr);

        copyProPtr->__sub.push_back(copyChainPtr);
    }

    return copyProPtr;
}


////////////////////////////////////////////////////////////////////////////////
// GetResidues
////////////////////////////////////////////////////////////////////////////////

vector<Residue *> Protein::getResidues()
{
    vector<Residue *> resPtrList;

    for (auto chainPtr: __sub)
    {
        for (auto resPtr: chainPtr->sub())
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

    for (auto chainPtr: __sub)
    {
        for (auto resPtr: chainPtr->sub())
        {
            for (auto atomPtr: resPtr->sub())
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

    for (auto chainPtr: __sub)
    {
        chainPtrMap.emplace(chainPtr->name(), chainPtr);
    }

    return chainPtrMap;
}


////////////////////////////////////////////////////////////////////////////////
// Dump
////////////////////////////////////////////////////////////////////////////////

Protein *Protein::dump(const string &dumpFilePath, const string &fileMode)
{
    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    for (auto chainPtr: __sub)
    {
        for (auto resPtr: *chainPtr)
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
                        chainPtr->name().c_str(),
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
                        chainPtr->name().c_str(),
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
    }

    fclose(fo);

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////

Protein::~Protein()
{
    for (auto subPtr: __sub)
    {
        delete subPtr;
    }
}


}  // End namespace PDBTools
