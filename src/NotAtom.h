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
using Eigen::RowVector3d;


////////////////////////////////////////////////////////////////////////////////
// Class __NotAtom
////////////////////////////////////////////////////////////////////////////////

template <typename SelfType, typename SubType>
class __NotAtom
{
public:

    // Iterator typedef
    typedef typename vector<SubType *>::iterator iterator;
    typedef typename vector<SubType *>::const_iterator const_iterator;
    typedef typename vector<SubType *>::reverse_iterator reverse_iterator;
    typedef typename vector<SubType *>::const_reverse_iterator const_reverse_iterator;


    // Iterator
    iterator begin()
    {
        return static_cast<SelfType *>(this)->sub.begin();
    }


    const_iterator begin() const
    {
        return static_cast<const SelfType *>(this)->sub.begin();
    }


    iterator end()
    {
        return static_cast<SelfType *>(this)->sub.end();
    }


    const_iterator end() const
    {
        return static_cast<const SelfType *>(this)->sub.end();
    }


    const_iterator cbegin() const
    {
        return static_cast<const SelfType *>(this)->sub.cbegin();
    }


    const_iterator cend() const
    {
        return static_cast<const SelfType *>(this)->sub.cend();
    }


    reverse_iterator rbegin()
    {
        return static_cast<SelfType *>(this)->sub.rbegin();
    }


    const_reverse_iterator rbegin() const
    {
        return static_cast<const SelfType *>(this)->sub.rbegin();
    }


    reverse_iterator rend()
    {
        return static_cast<SelfType *>(this)->sub.rend();
    }


    const_reverse_iterator rend() const
    {
        return static_cast<const SelfType *>(this)->sub.rend();
    }


    const_reverse_iterator crbegin() const
    {
        return static_cast<const SelfType *>(this)->sub.crbegin();
    }


    const_reverse_iterator crend() const
    {
        return static_cast<const SelfType *>(this)->sub.crend();
    }


    // FilterAtoms
    vector<Atom *> FilterAtoms(const unordered_set<string> &atomNameSet = {"CA"});
    vector<const Atom *> FilterAtoms(const unordered_set<string> &atomNameSet = {"CA"}) const;


    // GetAtomsCoord
    vector<RowVector3d *> GetAtomsCoord();
    vector<const RowVector3d *> GetAtomsCoord() const;


    // FilterAtomsCoord
    vector<RowVector3d *> FilterAtomsCoord(const unordered_set<string> &atomNameSet = {"CA"});

    vector<const RowVector3d *> FilterAtomsCoord(
        const unordered_set<string> &atomNameSet = {"CA"}) const;


    // Dumps
    string Dumps() const;


    // center
    RowVector3d center() const;


    // MoveCenter
    SelfType *MoveCenter();


    // seq
    string seq() const;


    // fasta
    string fasta(const string &titleStr = "") const
    {
        return ">" + (titleStr == "" ? static_cast<const SelfType *>(this)->name :
            titleStr) + "\n" + seq() + "\n";
    }


    // DumpFasta
    SelfType *DumpFasta(const string &dumpFilePath,
        const string &titleStr = "", const string &fileMode = "w");

    const SelfType *DumpFasta(const string &dumpFilePath,
        const string &titleStr = "", const string &fileMode = "w") const;


    // Renum Residues
    SelfType *RenumResidues(int startNum = 1);


    // Renum Atoms
    SelfType *RenumAtoms(int startNum = 1);


    // Append
    template <typename T>
    SelfType *Append(const vector<T *> &subPtrList);


    // Insert
    template <typename T>
    SelfType *Insert(iterator insertIter, const vector<T *> &subPtrList);


    // Move
    SelfType *Move(const vector<SubType *> &subPtrList);


    // MoveInsert
    SelfType *MoveInsert(iterator insertIter, const vector<SubType *> &subPtrList);


    // RemoveAlt
    SelfType *RemoveAlt();
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_NOT_ATOM_H
