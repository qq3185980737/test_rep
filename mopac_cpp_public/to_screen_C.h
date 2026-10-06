// to_screen_C.h — C++ mapping of Fortran module "to_screen_C".
#pragma once
#include <vector>
namespace to_screen_C {
// travel(k) - distance traveled during vibration k (1-based).
extern std::vector<double> travel;
// redmas(k,1) reduced mass; redmas(k,2) effective mass.
extern std::vector<std::vector<double>> redmas;
// rot(1..3) rotational constants (cm^-1); xyzmom moments of inertia.
extern std::vector<double> rot;
extern std::vector<double> xyzmom;
extern std::vector<double> freq;    // vibrational frequencies (force.F90)
extern std::vector<double> cnorml;  // normal-mode eigenvectors, packed (force.F90)
extern std::vector<double> dipt;    // vibrational dipole intensities (force.F90)

extern std::vector<std::vector<double>> dip;
}