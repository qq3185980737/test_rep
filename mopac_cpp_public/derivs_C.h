// derivs_C.h — C++ mapping of Fortran module "derivs_C" (subset used).
#pragma once
#include <vector>
namespace derivs_C {
extern std::vector<double> b;        // derivative scratch
extern std::vector<double> ab;       // (alpha,beta) derivative scratch
extern std::vector<double> fb_ci;    // CI Fock (beta) alias
extern std::vector<double> fmat;     // f matrix scratch
extern std::vector<double> hmat;     // h matrix scratch
extern std::vector<double> wmat;     // w matrix scratch
extern std::vector<double> aidref;   // reference aid array
extern std::vector<double> work2;    // 2-D work array, flat
}
