// atomrs.h — C++ translation of MOPAC 2016 "atomrs.F90".
#pragma once
#include <string>
#include <vector>

// atomrs: identify a protein residue starting at iatom, walk the bond graph,
// determine the residue type, and label atoms via txtatm.
// Arrays keep Fortran 1-based indexing (index 0 is padding).
void atomrs(std::vector<int>& lused, std::vector<bool>& ioptl, int& ires,
            int n1, int io, int uni_res, bool first_res);

// peptide_n(l): true if atom l looks like a backbone peptide nitrogen
// (three bonds: two C, one H, with one C attached to a terminal O).
bool peptide_n(int l);
