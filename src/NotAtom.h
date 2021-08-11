/*
    NotAtom.h
    =========
        Class __NotAtom header.
*/

#pragma once

#include <string>
#include <vector>
#include <unordered_set>
#include <initializer_list>
#include <Eigen/Dense>
#include "Predeclaration.h"

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
    vector<Atom *> filterAtoms(
        const unordered_set<string> &atomNameSet = {"CA"});


    // GetAtomsCoord
    Matrix<double, Dynamic, 3> getAtomsCoord();


    // FilterAtomsCoord
    Matrix<double, Dynamic, 3> filterAtomsCoord(
        const unordered_set<string> &atomNameSet = {"CA"});


    // center
    RowVector3d center();


    // MoveCenter
    SelfType *moveCenter();


    // seq
    string seq();


    // fasta
    string fasta(const string &titleStr = "");


    // DumpFasta
    SelfType *dumpFasta(const string &dumpFilePath,
        const string &titleStr = "", const string &fileMode = "w");


    // Renum Residues
    SelfType *renumResidues(int startNum = 1);


    // Renum Atoms
    SelfType *renumAtoms(int startNum = 1);


    // Append
    SelfType *append(SubType *subPtr, bool copyBool = true);


    // Insert
    SelfType *insert(typename vector<SubType *>::iterator insertIter,
        SubType *subPtr, bool copyBool = true);


    // RemoveAlt
    SelfType *removeAlt();


    // Dumps
    string dumps();
};


}  // End namespace PDBTools
