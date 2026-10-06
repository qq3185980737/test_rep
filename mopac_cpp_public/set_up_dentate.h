// set_up_dentate.h — C++ translation of "set_up_dentate.F90".
#pragma once
#include <vector>

// set_up_dentate: works out which atoms are connected (within 1.1x the sum
// of covalent radii).  On exit: nbonds(i) = number of attached atoms;
// ibonds(1..nbonds(i), i) = attached atom numbers.
void set_up_dentate();

// nsp2_correction: molecular-mechanics correction for nitrogen atoms with
// exactly three ligands (PM6/PM7 only).
double nsp2_correction();

// nsp2_atom_correction: penalty for non-planarity about central atom n.
// vectors(1..3, atom) Cartesian coordinates (1-based).
double nsp2_atom_correction(const std::vector<std::vector<double>>& vectors,
                            int n, int i, int j, int k);

// C_triple_bond_C: stabilization for acetylenic bonds (PM6/PM7 only).
double C_triple_bond_C();

// Si_O_H_Correction: bending perturbation for Si-O-H structures.
double Si_O_H_Correction();

// Si_O_H_bond_correction: Gaussian penalty when Si-O > 1.7 A and O-H > 1.0 A
// and the Si-O-H angle deviates from 125 degrees.
double Si_O_H_bond_correction(const std::vector<std::vector<double>>& coord,
                              int Si, int O, int H);
