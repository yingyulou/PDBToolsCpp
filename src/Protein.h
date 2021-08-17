/*
    Protein.h
    =========
        Class Protein header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
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


////////////////////////////////////////////////////////////////////////////////
// Class Protein
////////////////////////////////////////////////////////////////////////////////

class Protein: public __NotAtom<Protein, Chain>
{
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


    // str
    string str() const;


    // Copy
    Protein *copy();


    // GetResidues
    vector<Residue *> getResidues();


    // GetAtoms
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
};


}  // End namespace PDBTools
