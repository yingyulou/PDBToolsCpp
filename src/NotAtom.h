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
#include <Eigen>
#include "Predecl.h"

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

    // Begin
    typename vector<SubType *>::iterator begin();


    // End
    typename vector<SubType *>::iterator end();


    // Filter Atoms
    vector<Atom *> filterAtoms(
        const unordered_set<string> &atomNameSet = {"CA"});


    // Get Atoms Coord
    Matrix<double, Dynamic, 3> getAtomsCoord();


    // Filter Atoms Coord
    Matrix<double, Dynamic, 3> filterAtomsCoord(
        const unordered_set<string> &atomNameSet = {"CA"});


    // Center
    RowVector3d center();


    // Move Center
    SelfType *moveCenter();


    // Seq
    string seq();


    // Fasta Str
    string fastaStr(const string &titleStr = "");


    // Dump Fasta
    SelfType *dumpFasta(const string &dumpFilePath,
        const string &fileMode = "w", const string &titleStr = "");


    // Renum Residues
    SelfType *renumResidues(int startNum = 1);


    // Renum Atoms
    SelfType *renumAtoms(int startNum = 1);


    // Append
    SelfType *append(SubType *subPtr, bool copyBool = true);


    // Insert
    SelfType *insert(typename vector<SubType *>::iterator insertIter,
        SubType *subPtr, bool copyBool = true);


    // Remove Alt
    SelfType *removeAlt();


    // Dump Str
    string dumpStr();
};


}  // End namespace PDBTools
