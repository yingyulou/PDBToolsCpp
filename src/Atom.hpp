/*
    Atom.hpp
    ========
        Class Atom implementation.
*/

#pragma once

#include <string>
#include <cstdio>
#include <cctype>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "Atom.h"
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
// Constructor
////////////////////////////////////////////////////////////////////////////////

Atom::Atom(const string &atomName, int atomNum, const RowVector3d &atomCoord,
    const string &atomAltLoc, const string &atomOccupancy,
    const string &atomTempFactor, const string &atomElement,
    const string &atomCharge, Residue *atomOwner):
    name (atomName),
    num  (atomNum),
    coord(atomCoord),
    alt  (atomAltLoc),
    occ  (atomOccupancy),
    tempF(atomTempFactor),
    ele  (atomElement),
    chg  (atomCharge),
    owner(atomOwner)
{
    if (atomOwner)
    {
        atomOwner->sub.push_back(this);
    }
}


////////////////////////////////////////////////////////////////////////////////
// str
////////////////////////////////////////////////////////////////////////////////

string Atom::str() const
{
    return (format("<Atom object: %d %s [%.3f, %.3f, %.3f], at 0x%p>") %
        num                                                            %
        name                                                           %
        coord[0]                                                       %
        coord[1]                                                       %
        coord[2]                                                       %
        this
    ).str();
}


////////////////////////////////////////////////////////////////////////////////
// Copy
////////////////////////////////////////////////////////////////////////////////

Atom *Atom::copy()
{
    return new Atom(name, num, coord, alt, occ, tempF, ele, chg);
}


////////////////////////////////////////////////////////////////////////////////
// operator-
////////////////////////////////////////////////////////////////////////////////

double Atom::operator-(const Atom &rhs) const
{
    return (coord - rhs.coord).norm();
}


////////////////////////////////////////////////////////////////////////////////
// Dump
////////////////////////////////////////////////////////////////////////////////

Atom *Atom::dump(const string &dumpFilePath, const string &fileMode)
{
    string chainName, resName, resIns, dumpStr;
    int resNum = 0;

    if (owner)
    {
        resName = owner->name;
        resNum  = owner->num;
        resIns  = owner->ins;

        if (owner->owner)
        {
            chainName = owner->owner->name;
        }
    }

    FILE *fo = fopen(dumpFilePath.c_str(), fileMode.c_str());

    if (isdigit(name[0]) || name.size() == 4)
    {
        fprintf(fo, "ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n",
            num,
            name.c_str(),
            alt.c_str(),
            resName.c_str(),
            chainName.c_str(),
            resNum,
            resIns.c_str(),
            coord[0],
            coord[1],
            coord[2],
            occ.c_str(),
            tempF.c_str(),
            ele.c_str(),
            chg.c_str()
        );
    }
    else
    {
        fprintf(fo, "ATOM  %5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n",
            num,
            name.c_str(),
            alt.c_str(),
            resName.c_str(),
            chainName.c_str(),
            resNum,
            resIns.c_str(),
            coord[0],
            coord[1],
            coord[2],
            occ.c_str(),
            tempF.c_str(),
            ele.c_str(),
            chg.c_str()
        );
    }

    fclose(fo);

    return this;
}


////////////////////////////////////////////////////////////////////////////////
// Dumps
////////////////////////////////////////////////////////////////////////////////

string Atom::dumps()
{
    string chainName, resName, resIns, dumpStr, formatStr;
    int resNum = 0;

    if (owner)
    {
        resName = owner->name;
        resNum  = owner->num;
        resIns  = owner->ins;

        if (owner->owner)
        {
            chainName = owner->owner->name;
        }
    }

    if (isdigit(name[0]) || name.size() == 4)
    {
        formatStr = "ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n";
    }
    else
    {
        formatStr = "ATOM  %5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s          %2s%2s\n";
    }

    dumpStr = (format(formatStr) %
        num                      %
        name                     %
        alt                      %
        resName                  %
        chainName                %
        resNum                   %
        resIns                   %
        coord[0]                 %
        coord[1]                 %
        coord[2]                 %
        occ                      %
        tempF                    %
        ele                      %
        chg
    ).str();

    return dumpStr;
}


}  // End namespace PDBTools
