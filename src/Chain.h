/*
    Chain.h
    =======
        Class Chain header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "NotAtom.h"
#include "NotProtein.h"
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
using std::ostream;


////////////////////////////////////////////////////////////////////////////////
// Class Chain
////////////////////////////////////////////////////////////////////////////////

class Chain: public __NotAtom<Chain, Residue>,
    public __NotProtein<Chain, Protein>
{
    // Friend
    friend ostream &operator<<(ostream &os, const Chain &chainObj);


public:

    // Constructor
    explicit Chain(const string &name = "", Protein *owner = nullptr);


    // Getter: __name
    string &name();


    // Getter: __owner
    Protein *owner();


    // Getter: __sub
    vector<Residue *> &sub();


    // Setter: __name
    Chain *name(const string &val);


    // Setter: __owner
    Chain *owner(Protein *val);


    // Setter: __sub
    Chain *sub(const vector<Residue *> &val);


    // Copy
    Chain *copy();


    // GetResidues
    vector<Residue *> getResidues();


    // GetAtoms
    vector<Atom *> getAtoms();


    // subMap
    unordered_map<string, Residue *> subMap();


    // Dump
    Chain *dump(const string &dumpFilePath, const string &fileMode = "w");


    // Destructor
    ~Chain();


private:

    // Attribute
    string __name;
    Protein *__owner;
    vector<Residue *> __sub;


    // str
    string __str() const;
};


}  // End namespace PDBTools
