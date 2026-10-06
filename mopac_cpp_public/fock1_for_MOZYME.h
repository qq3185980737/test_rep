// fock1_for_MOZYME.h — C++ translation of MOPAC 2016.
#pragma once
#include <vector>
void fock1_for_MOZYME(std::vector<double>& f, const std::vector<double>& ptot,
                      const std::vector<std::vector<double>>& w,
                      int& kr, int iab, int ilim);
