// deri22.h — C++ translation of MOPAC 2016 "deri22.F90".
#pragma once
#include <vector>

// Build (D-A)*B supervector ab from density derivative b, plus the
// CI-active Fock diagonal blocks fci. b/ab/fci are column-major supervectors
// over the relaxation basis; bcol/abcol/fcicol select the 1-based column
// (F90 b(1,j) / ab(1,j) / fci(1,j)).
void deri22(const std::vector<std::vector<double>>& c,
            std::vector<double>& b, int bcol,
            std::vector<std::vector<double>>& work,
            std::vector<double>& foc2,
            std::vector<double>& ab, int abcol, int minear,
            std::vector<double>& fci, int fcicol,
            std::vector<double>& w,
            const std::vector<double>& diag,
            const std::vector<double>& scalar, int ninear);
