/*
    Atom.h
    ======
        Class Atom header.
*/

#pragma once

#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "NotProtein.h"
#include "Residue.h"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::ostream;
using Eigen::RowVector3d;


////////////////////////////////////////////////////////////////////////////////
// Class Atom
////////////////////////////////////////////////////////////////////////////////

class Atom: public __NotProtein<Atom, Residue>
{
    // Friend
    friend ostream &operator<<(ostream &os, const Atom &atomObj);


public:

    // Constructor
    explicit Atom(const string &name = "", int num = 0,
        const RowVector3d &coord = RowVector3d::Zero(),
        const string &alt = "", const string &occ = "",
        const string &tempF = "", const string &ele = "",
        const string &chg = "", Residue *owner = nullptr);


    // Getter: __name
    string &name();


    // Getter: __num
    int num();


    // Getter: __coord
    RowVector3d &coord();


    // Getter: __alt
    string &alt();


    // Getter: __occ
    string &occ();


    // Getter: __tempF
    string &tempF();


    // Getter: __ele
    string &ele();


    // Getter: __chg
    string &chg();


    // Getter: __owner
    Residue *owner();


    // Setter: __name
    Atom *name(const string &val);


    // Setter: __num
    Atom *num(int val);


    // Setter: __coord
    Atom *coord(const RowVector3d &val);


    // Setter: __alt
    Atom *alt(const string &val);


    // Setter: __occ
    Atom *occ(const string &val);


    // Setter: __tempF
    Atom *tempF(const string &val);


    // Setter: __ele
    Atom *ele(const string &val);


    // Setter: __chg
    Atom *chg(const string &val);


    // Setter: __owner
    Atom *owner(Residue *val);


    // Copy
    Atom *copy();


    // operator-
    double operator-(const Atom &rhs) const;


    // Dump
    Atom *dump(const string &dumpFilePath, const string &fileMode = "w");


    // Dump Str
    string dumpStr();


private:

    // Attribute
    string __name;
    int __num;
    RowVector3d __coord;
    string __alt;
    string __occ;
    string __tempF;
    string __ele;
    string __chg;
    Residue *__owner;

    // str
    string __str() const;
};


}  // End namespace PDBTools
