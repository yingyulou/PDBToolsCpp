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
#include <boost/format.hpp>
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
using boost::format;


////////////////////////////////////////////////////////////////////////////////
// Class Protein
////////////////////////////////////////////////////////////////////////////////

class Protein: public __StructBase<Protein>, public __NotAtom<Protein, Chain>
{
public:

    // Attribute
    string name;
    vector<Chain *> sub;


    // Constructor
    explicit Protein(const string &proteinID = "");


    // str
    string str() const
    {
        return (format("<Protein object: %s, at 0x%p>") % name % this).str();
    }


    // Copy
    Protein *Copy() const;


    // GetResidues
    vector<Residue *> GetResidues();
    vector<const Residue *> GetResidues() const;


    // GetAtoms
    vector<Atom *> GetAtoms();
    vector<const Atom *> GetAtoms() const;


    // subMap
    unordered_map<string, Chain *> subMap();
    unordered_map<string, const Chain *> subMap() const;


    // Destructor
    ~Protein();
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_PROTEIN_H
