// to_screen_C.cpp — storage.
#include "to_screen_C.h"
namespace to_screen_C {
std::vector<double> travel;
std::vector<std::vector<double>> redmas;
std::vector<double> rot(4, 0.0);
std::vector<double> freq;   // vibrational frequencies (force.F90)
std::vector<double> cnorml; // normal-mode eigenvectors, packed (force.F90)
std::vector<double> dipt;   // vibrational dipole intensities (force.F90)
std::vector<double> xyzmom(4, 0.0);

std::vector<std::vector<double>> dip(4, std::vector<double>(3,0.0));
}