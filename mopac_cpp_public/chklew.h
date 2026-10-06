// chklew.h — C++ translation of MOPAC 2016 "chklew.F90".
#pragma once
#include <vector>

// add_Lewis_element: record a Lewis-structure element and update bookkeeping.
// element_type is incremented for sigma/lone/pi bond types.
void add_Lewis_element(int atom_i, int atom_j, int charge, int& element_type);

// ring5 / arom / arom2: five-membered-ring and aromaticity helpers.
void ring5(int i, const std::vector<int>& mb, std::vector<int>& ir5);
bool arom(int ii, int jj, const std::vector<int>& mpii);
bool arom2(int ii, int jj, const std::vector<int>& mpii);

// chklew: build the Lewis structure (sigma, lone pairs, pi bonds).
void chklew(std::vector<int>& mb, std::vector<int>& numbon, int& l,
            int large, bool debug);
