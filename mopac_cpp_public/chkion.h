// chkion.h — C++ translation of MOPAC 2016 "chkion.F90".
#pragma once
#include <string>
#include <vector>

// chkion: determine which atoms are ionized from Lewis structure.
// atom_charge(i): user-specified charge override (' ', '+', '-', '0').
void chkion(std::vector<int>& ox_calc, int& n_lone_pairs,
            const std::vector<char>& atom_charge);
