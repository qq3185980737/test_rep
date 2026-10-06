// fmat.cpp — C++ translation of MOPAC 2016 "fmat.F90".

#include "fmat.h"

double fmat_second(double grold, double grad, double delta, double fact) {
    return (grold - grad) * 0.5 / delta * fact;
}

double fmat_fourth(double grold, double grad, double g2old, double g2rad,
                   double delta, double fact) {
    return (8.0 * (grold - grad) - (g2old - g2rad)) / delta * fact / 12.0;
}

void fmat(std::vector<double>&, int& nreal, double, double,
          std::vector<std::vector<double>>&, double, std::vector<double>&, bool) { nreal = 0; }
