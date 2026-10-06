// fock1.h — C++ translation of MOPAC 2016 "fock1.F90".
#pragma once
#include <vector>
void fock1(std::vector<double>& Fmat, const std::vector<double>& ptot,
           const std::vector<double>& PAa, const std::vector<double>& PBb);
