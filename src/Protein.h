/*
    Protein.h
    =========
        Class Protein header.
*/

#ifndef __PDBTOOLS_PROTEIN_H
#define __PDBTOOLS_PROTEIN_H

#include <string>
#include <vector>
#include <unordered_map>
#include "BaseStruct.h"
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

class Protein: public __BaseStruct<Protein>, public __NotAtom<Protein, Chain>
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


#endif  // __PDBTOOLS_PROTEIN_H
