// lewis.h — C++ translation of lewis.F90 (MOPAC 2016).
#pragma once
#include <string>

// Build the Lewis structure: raw topography from set_up_dentate is tidied by
// filters (max bonds, hydrogen-bond checks, nitro-group O-O bridge, CVB).
// On exit nbonds/ibonds hold the final connectivity (max 4 bonds per atom).
void lewis(bool use_cvs);

// Remove the longest bond from atom i (used to cap Fe/Ni/Pd/Pt at 5 bonds).
void remove_bond(int i);

// Check hydrogen bonding: nearest non-H neighbour, remove bridge bonds,
// report bad positions (ibad count).  Only active for MOZYME runs.
void check_h(int& ibad);

// Process the CVB keyword: add/delete bonds defined as CVB(j,l) or
// CVB("label","label"); verifies geometry and atom labels.
void check_CVS(bool let);

// Convert a PDB/Jmol atom label in "text" (1-based position j_in after the
// opening quote) to an atom number.  On success text is rewritten with the
// number.  m receives the compressed label length (or the Jmol ']' position).
void txt_to_atom_no(std::string& text, int j_in, bool let, int& m);
