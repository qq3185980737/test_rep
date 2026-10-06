// add_hydrogen_atoms.h — C++ translation of MOPAC 2016
// "add_hydrogen_atoms.F90" (public entry points only).
//
// The file defines several subprograms; only the ones called from outside
// this file are declared here. The rest (h_type, add_a_*_hydrogen_atom,
// aromatic, aromatic_5, near_a_metal) are file-local helpers.
#pragma once

#include <vector>

// Adds hydrogen atoms to a system, intended for converting PDB files into
// input files suitable for MOPAC. Operates entirely on module state.
void add_hydrogen_atoms();

// Detects and prints unusually short hydrogen bonds.
void bridge_H();

// Re-sets "breaks" after RESEQ / ADD_H moves residues.
void reset_breaks();

// Public wrapper: places a hydrogen atom on a 4-coordinate atom icc
// (SP3 geometry). Forwarded to the file-local add_a_sp3_hydrogen_atom.
// Called from geochk.F90 ("site").
void add_a_sp3_hydrogen_atom_ext(int icc, int nb_icc, int nc_icc, int nd_icc,
                                 double bond_length,
                                 const std::vector<int>& metals, int nmetals);
