/*
    Atom.h
    ======
        Class Atom header.
*/

#ifndef __PDBTOOLS_ATOM_H
#define __PDBTOOLS_ATOM_H

#include <string>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "StructBase.h"
#include "NotProtein.h"
#include "Residue.h"

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using boost::format;
using Eigen::RowVector3d;


////////////////////////////////////////////////////////////////////////////////
// Class Atom
////////////////////////////////////////////////////////////////////////////////

class Atom: public __StructBase<Atom>, public __NotProtein<Atom, Residue>
{
public:

    // Attribute
    string name;
    int num;
    RowVector3d coord;
    string alt;
    string occ;
    string tempF;
    string ele;
    string chg;
    Residue *owner;


    // Constructor
    explicit Atom(const string &atomName = "", int atomNum = 0,
        const RowVector3d &atomCoord = RowVector3d(0., 0., 0.),
        const string &atomAltLoc = "", const string &atomOccupancy = "",
        const string &atomTempFactor = "", const string &atomElement = "",
        const string &atomCharge = "", Residue *atomOwner = nullptr);


    // str
    string str() const
    {
        return (format("<Atom object: %d %s [%.3f, %.3f, %.3f], at 0x%p>") %
            num % name % coord[0] % coord[1] % coord[2] % this).str();
    }


    // Copy
    Atom *Copy() const
    {
        return new Atom(name, num, coord, alt, occ, tempF, ele, chg);
    }


    // Dumps
    string Dumps() const;


    // operator-
    double operator-(const Atom &rhs) const
    {
        return (coord - rhs.coord).norm();
    }
};


}  // End namespace PDBTools


#endif  // __PDBTOOLS_ATOM_H
