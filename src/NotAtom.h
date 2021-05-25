/*
    NotAtom.h
    =========
        Class __NotAtom header.
*/

#ifndef __PDBTOOLS_NOT_ATOM_H
#define __PDBTOOLS_NOT_ATOM_H

#include <string>
#include <vector>
#include <unordered_set>
#include <initializer_list>
#include <Eigen/Dense>
#include "Struct.h"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::vector;
using std::unordered_set;
using std::initializer_list;
using Eigen::RowVector3d;
using Eigen::Matrix;
using Eigen::Dynamic;


////////////////////////////////////////////////////////////////////////////////
// Class __NotAtom
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
class __NotAtom
{
public:

    // Iterator
    typename vector<SubType *>::iterator begin();
    typename vector<SubType *>::iterator end();


    // FilterAtoms
    vector<Atom *> FilterAtoms(
        const unordered_set<string> &atomNameSet = {"CA"});


    // GetAtomsCoord
    Matrix<double, Dynamic, 3> GetAtomsCoord();


    // FilterAtomsCoord
    Matrix<double, Dynamic, 3> FilterAtomsCoord(
        const unordered_set<string> &atomNameSet = {"CA"});


    // Dumps
    string Dumps();


    // center
    RowVector3d center();


    // MoveCenter
    SelfType *MoveCenter();


    // seq
    string seq();


    // fasta
    string fasta(const string &titleStr = "");


    // DumpFasta
    SelfType *DumpFasta(const string &dumpFilePath,
        const string &titleStr = "", const string &fileMode = "w");


    // Renum Residues
    SelfType *RenumResidues(int startNum = 1);


    // Renum Atoms
    SelfType *RenumAtoms(int startNum = 1);


    // Append
    SelfType *Append(SubType *subPtr, bool copyBool = true);


    // Insert
    SelfType *Insert(typename vector<SubType *>::iterator insertIter,
        SubType *subPtr, bool copyBool = true);


    // RemoveAlt
    SelfType *RemoveAlt();
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_NOT_ATOM_H
