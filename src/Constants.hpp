/*
    Constants.hpp
    =============
        Constants define.
*/

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <regex>

namespace PDBTools
{

////////////////////////////////////////////////////////////////////////////////
// Using
////////////////////////////////////////////////////////////////////////////////

using std::string;
using std::unordered_map;
using std::vector;
using std::regex;


////////////////////////////////////////////////////////////////////////////////
// H Atom Name Regex
////////////////////////////////////////////////////////////////////////////////

const regex __H_RE(R"(^\d*H)");


////////////////////////////////////////////////////////////////////////////////
// CompNum Regex
////////////////////////////////////////////////////////////////////////////////

const regex __COMP_NUM_RE(R"((\d+)([A-Za-z]?))");


////////////////////////////////////////////////////////////////////////////////
// Dihedral Enum
////////////////////////////////////////////////////////////////////////////////

enum class DIH
{
    PHI = 0,
    L   = 0,
    PSI = 1,
    R   = 1,
};


////////////////////////////////////////////////////////////////////////////////
// Side Enum
////////////////////////////////////////////////////////////////////////////////

enum class SIDE
{
    N = 0,
    L = 0,
    C = 1,
    R = 1,
};


////////////////////////////////////////////////////////////////////////////////
// Residue Name 3 Letters => 1 Letter
////////////////////////////////////////////////////////////////////////////////

const unordered_map<string, string> RESIDUE_NAME_THREE_TO_ONE_MAP
{
    {"ALA", "A"},
    {"ARG", "R"},
    {"ASN", "N"},
    {"ASP", "D"},
    {"CYS", "C"},
    {"GLN", "Q"},
    {"GLU", "E"},
    {"GLY", "G"},
    {"HIS", "H"},
    {"ILE", "I"},
    {"LEU", "L"},
    {"LYS", "K"},
    {"MET", "M"},
    {"PHE", "F"},
    {"PRO", "P"},
    {"SER", "S"},
    {"THR", "T"},
    {"TRP", "W"},
    {"TYR", "Y"},
    {"VAL", "V"},
    {"UNK", "X"},
};


////////////////////////////////////////////////////////////////////////////////
// Residue Name 1 Letter => 3 Letters
////////////////////////////////////////////////////////////////////////////////

const unordered_map<string, string> RESIDUE_NAME_ONE_TO_THREE_MAP
{
    {"A", "ALA"},
    {"R", "ARG"},
    {"N", "ASN"},
    {"D", "ASP"},
    {"C", "CYS"},
    {"Q", "GLN"},
    {"E", "GLU"},
    {"G", "GLY"},
    {"H", "HIS"},
    {"I", "ILE"},
    {"L", "LEU"},
    {"K", "LYS"},
    {"M", "MET"},
    {"F", "PHE"},
    {"P", "PRO"},
    {"S", "SER"},
    {"T", "THR"},
    {"W", "TRP"},
    {"Y", "TYR"},
    {"V", "VAL"},
    {"X", "UNK"},
};


////////////////////////////////////////////////////////////////////////////////
// Residue Side Chain Rotation Atoms Name
////////////////////////////////////////////////////////////////////////////////

const unordered_map<string, vector<vector<string>>>
    __RESIDUE_SIDE_CHAIN_ROTATION_ATOMS_NAME_MAP
{
    {"ALA", {}},

    {"ARG", {
        {"N", "CA", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
        {"CA", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
        {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
        {"CG", "CD", "NE", "CZ", "NH1", "NH2"},
    }},

    {"ASN", {
        {"N", "CA", "CB", "CG", "OD1", "ND2"},
        {"CA", "CB", "CG", "OD1", "ND2"},
    }},

    {"ASP", {
        {"N", "CA", "CB", "CG", "OD1", "OD2"},
        {"CA", "CB", "CG", "OD1", "OD2"},
    }},

    {"CYS", {
        {"N", "CA", "CB", "SG"},
    }},

    {"GLN", {
        {"N", "CA", "CB", "CG", "CD", "OE1", "NE2"},
        {"CA", "CB", "CG", "CD", "OE1", "NE2"},
        {"CB", "CG", "CD", "OE1", "NE2"},
    }},

    {"GLU", {
        {"N", "CA", "CB", "CG", "CD", "OE1", "OE2"},
        {"CA", "CB", "CG", "CD", "OE1", "OE2"},
        {"CB", "CG", "CD", "OE1", "OE2"},
    }},

    {"GLY", {}},

    {"HIS", {
        {"N", "CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        {"CA", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    }},

    {"ILE", {
        {"N", "CA", "CB", "CG1", "CG2", "CD1"},
        {"CA", "CB", "CG1", "CD1"},
    }},

    {"LEU", {
        {"N", "CA", "CB", "CG", "CD1", "CD2"},
        {"CA", "CB", "CG", "CD1", "CD2"},
    }},

    {"LYS", {
        {"N", "CA", "CB", "CG", "CD", "CE", "NZ"},
        {"CA", "CB", "CG", "CD", "CE", "NZ"},
        {"CB", "CG", "CD", "CE", "NZ"},
        {"CG", "CD", "CE", "NZ"},
    }},

    {"MET", {
        {"N", "CA", "CB", "CG", "SD", "CE"},
        {"CA", "CB", "CG", "SD", "CE"},
        {"CB", "CG", "SD", "CE"},
    }},

    {"PHE", {
        {"N", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
        {"CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    }},

    {"PRO", {
        {"N", "CA", "CB", "CG", "CD"},
        {"CA", "CB", "CG", "CD"},
    }},

    {"SER", {
        {"N", "CA", "CB", "OG"},
    }},

    {"THR", {
        {"N", "CA", "CB", "OG1", "CG2"},
    }},

    {"TRP", {
        {"N", "CA", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
        {"CA", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    }},

    {"TYR", {
        {"N", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
        {"CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
        {"N", "CA", "CB", "OG1", "CG2"},
    }},

    {"VAL", {
        {"N", "CA", "CB", "CG1", "CG2"},
    }},
};


}  // End namespace PDBTools
