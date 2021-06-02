/*
    Chain.h
    =======
        Class Chain header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "BaseStruct.h"
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


////////////////////////////////////////////////////////////////////////////////
// Class Chain
////////////////////////////////////////////////////////////////////////////////

class Chain: public __BaseStruct<Chain>, public __NotAtom<Chain, Residue>,
    public __NotProtein<Chain, Protein>
{
public:

    // Attribute
    string name;
    Protein *owner;
    vector<Residue *> sub;


    // Constructor
    explicit Chain(const string &chainName = "", Protein *chainOwner = nullptr);


    // str
    string str() const;


    // Copy
    Chain *Copy();


    // GetResidues
    vector<Residue *> GetResidues();


    // GetAtoms
    vector<Atom *> GetAtoms();


    // subMap
    unordered_map<string, Residue *> subMap();


    // Destructor
    ~Chain();
};


}  // End namespace PDBTools
