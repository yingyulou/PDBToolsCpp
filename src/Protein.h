/*
    Protein.h
    =========
        Class Protein header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "StructBase.h"
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

class Protein: public __StructBase<Protein>, public __NotAtom<Protein, Chain>
{
public:

    // Attribute
    string name;
    int model;
    vector<Chain *> sub;


    // Constructor
    explicit Protein(const string &proteinID = "", int modelNum = 0);


    // str
    string str() const;


    // Copy
    Protein *Copy();


    // GetResidues
    vector<Residue *> GetResidues();


    // GetAtoms
    vector<Atom *> GetAtoms();


    // subMap
    unordered_map<string, Chain *> subMap();


    // Destructor
    ~Protein();
};


}  // End namespace PDBTools
