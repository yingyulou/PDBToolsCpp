/*
    Protein.h
    =========
        Class Protein header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "NotAtom.h"
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
using std::ostream;


////////////////////////////////////////////////////////////////////////////////
// Class Protein
////////////////////////////////////////////////////////////////////////////////

class Protein: public __NotAtom<Protein, Chain>
{
    // Friend
    friend ostream &operator<<(ostream &os, const Protein &proObj);


public:

    // Constructor
    explicit Protein(const string &name = "", int model = 0);


    // Getter: __name
    string &name();


    // Getter: __model
    int model();


    // Getter: __sub
    vector<Chain *> &sub();


    // Setter: __name
    Protein *name(const string &val);


    // Setter: __model
    Protein *model(int val);


    // Setter: __sub
    Protein *sub(const vector<Chain *> &val);


    // Copy
    Protein *copy();


    // Get Residues
    vector<Residue *> getResidues();


    // Get Atoms
    vector<Atom *> getAtoms();


    // subMap
    unordered_map<string, Chain *> subMap();


    // Dump
    Protein *dump(const string &dumpFilePath, const string &fileMode = "w");


    // Destructor
    ~Protein();


private:

    // Attribute
    string __name;
    int __model;
    vector<Chain *> __sub;


    // str
    string __str() const;
};


}  // End namespace PDBTools
