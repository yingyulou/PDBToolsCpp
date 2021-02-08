/*
    Chain.h
    =======
        Class Chain header.
*/

#ifndef __PDBTOOLS_CHAIN_H
#define __PDBTOOLS_CHAIN_H

#include <string>
#include <vector>
#include <unordered_map>
#include "StructBase.h"
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

class Chain: public __StructBase<Chain>, public __NotAtom<Chain, Residue>,
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
    Chain *Copy() const;


    // GetResidues
    vector<Residue *> GetResidues() { return sub; }
    vector<const Residue *> GetResidues() const;


    // GetAtoms
    vector<Atom *> GetAtoms();
    vector<const Atom *> GetAtoms() const;


    // subMap
    unordered_map<string, Residue *> subMap();
    unordered_map<string, const Residue *> subMap() const;


    // Destructor
    ~Chain();
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_CHAIN_H
