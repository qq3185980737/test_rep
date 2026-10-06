// ef_C.h — globals from ef.C (extended for setup_mopac_arrays).
#pragma once
#include <vector>
namespace ef_C {
extern int nstep, negreq, iprnt, ef_mode;
extern double ddx, xlamd, xlamd0, skal, rmin, rmax, omin; extern double x0, x1, x2; extern int iloop; extern std::vector<std::vector<double>> alparm;
// EF optimizer workspaces (setup_mopac_arrays).
extern std::vector<std::vector<double>> hess;   // Hessian (3*nvar x 3*nvar)
extern std::vector<std::vector<double>> bmat;   // ef B-matrix (alias of ef_C bmat)
extern std::vector<double> u;                   // unitary transform scratch
extern std::vector<double> oldhss, oldu;        // old Hessian / old u
extern std::vector<std::vector<double>> pmat;   // p matrix
extern std::vector<std::vector<double>> uc;     // u * c
extern std::vector<double> hessc;               // Hessian * c
extern std::vector<std::vector<double>> oldf;   // old f matrix
extern std::vector<std::vector<double>> d;      // d matrix
extern std::vector<std::vector<double>> vmode;  // vibrational modes
}
