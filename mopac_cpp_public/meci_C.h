// meci_C.h — C++ translation of the Fortran module "meci_C" (MOPAC 2016).
//
// Original module (excerpt of the members used so far):
//   module meci_C
//     USE vast_kind_param, ONLY:  double
//     ...
//     real(double), dimension(:), allocatable :: occa  ! Initial M.O. occupancy
//   end module meci_C
//
// Only the members actually referenced by already-translated sources are
// declared here; the header will be extended as further .F90 files are
// ported. The namespace name matches the Fortran module name.
#pragma once

#include <vector>

namespace meci_C {

// occa — Initial M.O. occupancy: number of electrons in each M.O. of the
// SCF in the active space. Fortran: real(double), dimension(:), allocatable.
extern std::vector<double> occa;

// is, iiloop, jloop — control integers used by the MECI configuration-loop
// bookkeeping (Fortran: integer :: is, iiloop, jloop). Always set by the
// surrounding MECI driver before use.
extern int is;
extern int iiloop;
extern int jloop;

// ispqr — Fortran 2-D array ispqr(lab, nmeci+1). It keeps its Fortran
// (1-based) indexing convention: row 0 and column 0 are padding, and
// Fortran accesses such as ispqr(iiloop, is) / ispqr(j, l+1) are written
// verbatim as ispqr[iiloop][is] / ispqr[j][l + 1]. The allocation (with
// the -99999 sentinel) happens in the not-yet-ported MECI driver.
extern std::vector<std::vector<int>> ispqr;
extern int lab;
extern int nmos;
extern int nstate;
extern int nbo[4];
extern int nelec;
extern std::vector<int> nalmat;
extern std::vector<std::vector<int>> microa;
extern std::vector<std::vector<int>> microb;
extern std::vector<double> conf;
extern int nmeci;   // 22 (Fortran parameter)
extern std::vector<double> dijkl; extern std::vector<double> xy;
// vectci: CI eigenvectors, column-major (lab rows x nstate cols).
extern std::vector<double> vectci;
// deltap: density correction workspace (norbs x nmos).
extern std::vector<std::vector<double>> deltap;
// meci driver state.
extern std::vector<std::vector<double>> rjkaa;
extern std::vector<std::vector<double>> rjkab;
extern std::vector<double> eig;
extern std::vector<int> ispin;
extern int maxci;
extern std::vector<double> cdiag;
extern double cdiagi;
extern double cif1, cif2;
extern int k, dummy;
extern int labsiz;
extern int root_requested;
extern int msdel;
extern std::vector<double> spin;
extern std::vector<int> nfa;
// eiga: active-space eigenvalues (length nmos). cimat: packed CI matrix.
extern std::vector<double> eiga;
extern std::vector<double> cimat;
}  // namespace meci_C
