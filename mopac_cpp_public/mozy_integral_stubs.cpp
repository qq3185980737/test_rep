// mozy_integral_stubs.cpp — extern "C" stubs for the closed-source Fortran
// integral kernels that the C++ ports of h1elec / rotate / outer2 call.
//
// Status (2026-09-28):
//  * rotatd / elenuc / reppd — real F90 sources exist (mndod.F90) and are
//    translated in mndod.cpp; nddo_to_point lives in solrot.F90 and is
//    translated in solrot.cpp.  These extern "C" forwarders adapt the C++
//    1-based vector conventions to the Fortran-calling drivers (rotate.cpp /
//    h1elec.cpp / outer2.cpp), which use 0-based flat buffers.
//  * diat_ — diat.F90 is fully translated (diat.cpp: diat driver + diat2
//    sp-path + ss general overlap).  This forwarder builds a 1-based padded
//    view of the Fortran 9x9 column-major block and calls the C++ driver.
//  * trunk / to_point — internal subroutines of solrot.F90 / mndod.F90, now
//    translated inline in outer1.cpp (trunk_cpp / to_point_cpp).  The extern
//    "C" forwarders have been removed; no caller references them any more.
//  * fillij / memory_error — real C++ translations exist (fillij.cpp,
//    memory_error.cpp); the duplicates that used to live here were removed.
#include "parameters_C.h"
#include "diat.h"
#include <vector>

extern "C" void diat_(int* ni, int* nj, double* xjuc, double* smat) {
  // Diatomic one-electron integrals (s/p/d block).
  std::vector<double> xj(xjuc, xjuc + 3);
  std::vector<std::vector<double>> di(10, std::vector<double>(10, 0.0));
  diat(*ni, *nj, xj, di);
  for (int i = 1; i <= 9; ++i)
    for (int j = 1; j <= 9; ++j) smat[(j - 1) * 9 + (i - 1)] = di[i][j];
}

extern void rotatd(int, int, const std::vector<double>&, const std::vector<double>&,
                std::vector<double>&, int&, double&);
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj,
                        double* w, int* kr, double* enuc) {
  // Forward to the real C++ translation (mndod.cpp, mndod.F90 rotatd).
  // C++ writes 1-based w[1..kr-1]; the Fortran caller reads 0-based w[0..kr-2].
  // C++ rotatd uses 1-based ci[1..3]/cj[1..3]; Fortran caller passes 0-based
  // 3-vectors.  Pad to size 4 so ci[3]/cj[3] stay in bounds.
  std::vector<double> civ(4, 0.0), cjv(4, 0.0);
  civ[1] = xi[0]; civ[2] = xi[1]; civ[3] = xi[2];
  cjv[1] = xj[0]; cjv[2] = xj[1]; cjv[3] = xj[2];
  std::vector<double> wv(2026, 0.0);
  int krv = *kr;
  double enucv = *enuc;
  rotatd(*ni, *nj, civ, cjv, wv, krv, enucv);
  for (int k = 0; k < krv - 1; ++k) w[k] = wv[k + 1];
  *kr = krv;
  *enuc = enucv;
}

extern void elenuc(int, int, int, int, std::vector<double>&);
extern "C" void elenuc_(int* ni, int* nj, int* li, int* lj, double* en) {
  // Forward to the real C++ translation (mndod.cpp, mndod.F90 elenuc).
  // Caller passes (ilow,iup,jlow,jup) as (ni,nj,li,lj): ia/ib=ni..nj, ja/jb=li..lj.
  // elenuc writes h[m] with m=(i*(i-1))/2+j, max m = jb*(jb+1)/2 = lj*(lj+1)/2.
  int n = (*lj) * (*lj + 1) / 2;
  std::vector<double> hv(172, 0.0);
  elenuc(*ni, *nj, *li, *lj, hv);
  for (int k = 0; k < n; ++k) en[k] = hv[k + 1];
}

extern void nddo_to_point(double*, double*, double*, double&, double, int, int);
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double* rij, int* ni, int* nj) {
  // Forward to the real C++ translation (solrot.cpp, solrot.F90 nddo_to_point).
  // C++ arrays are 1-based; the Fortran caller uses 0-based flat buffers.
  nddo_to_point(w - 1, e1b - 1, e2a - 1, *enuc, *rij, *ni, *nj);
}

extern void reppd(int, int, double, std::vector<double>&, double&);
extern "C" void reppd_(int* ni, int* nj, double* rij, double* ri, double* gab) {
  // Forward to the real C++ translation (mndod.cpp, mndod.F90 reppd).
  // C++ ri is 1-based (22 entries); Fortran caller reads ri 0-based.
  std::vector<double> rv(23, 0.0);
  double gabv = 0.0;
  reppd(*ni, *nj, *rij, rv, gabv);
  for (int k = 0; k < 22; ++k) ri[k] = rv[k + 1];
  *gab = gabv;
}
