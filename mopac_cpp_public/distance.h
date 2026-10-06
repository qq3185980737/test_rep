// distance.h — C++ translation of distance/angle/torsion functions from
// H_bond_correction_bits.F90 (lines 267-299).
#pragma once
#include <vector>

// Interatomic distance between atoms a and b, with lattice translation
// (id != 0) taken into account.
double distance(int a, int b);

// Bond angle at b between a and c (radians).
double angle(int a, int b, int c);

// Dihedral angle between planes (i,j,k) and (j,k,l) (radians).
double torsion(int i, int j, int k, int l);
