// buildf.h — C++ translation of MOPAC 2016 "buildf.F90".
#pragma once
#include <vector>

// buildf: combine a partial Fock matrix with the one-electron matrix h,
// then hand off to fock2z.
//   mode=-1: f = partf - h
//   mode= 0: f = h
//   mode= 1: f = partf + h
void buildf(std::vector<double>& f, const std::vector<double>& partf, int mode);
