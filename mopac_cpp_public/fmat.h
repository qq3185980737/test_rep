// fmat.h — C++ translation of MOPAC 2016 "fmat.F90".
#pragma once
#include <vector>

// 2nd-order central difference of gradient -> Hessian element.
double fmat_second(double grold, double grad, double delta, double fact);
// 4th-order correction (precis).
double fmat_fourth(double grold, double grad, double g2old, double g2rad,
                   double delta, double fact);
void fmat(std::vector<double>& fmatrx, int& nreal, double tscf, double tder,
          std::vector<double>& deldip, double heat,
          std::vector<double>& evecs, bool ts);
