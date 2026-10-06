// esp_C.h — C++ translation of MOPAC 2016 "esp_C.F90" (ESP fit module data).
#pragma once

#include <vector>

namespace esp_C {
// Scalars (Fortran module globals)
extern int nesp, idip, iz, ipx, isc, is_esp, icd, ipe, npr, ic, ip, ncc;
extern double dens, scale, cf, rms, rrms, dx, dy, dz, den;

// fv(0:8,821)
extern double fv[9][822];
// dex(-1:96) -> dex[i + 1] == dex(i)
extern double dex[98];
// tf(0:2)
extern double tf[3];
// fac(0:7)
extern double fac[8];
// ixn..jzn(10), 1-based -> [11]
extern int ixn[11], iyn[11], izn[11], jxn[11], jyn[11], jzn[11];

// Allocatables -> std::vector (1-D)
extern std::vector<int> ind, itemp, ird, indc;
extern std::vector<double> cc, ex, temp, rad, es, b_esp, esp_array, cesp,
    cespml, al, qesp, td, dx_array, dy_array, dz_array, ptd, pexs, pce,
    pexpn, pewcx, pewcy, pewcz, pf0, pf1, pf2, pmtd, fc, qsc, espx, espp;
// Allocatables -> std::vector (2-D)
extern std::vector<std::vector<int>> iam;
extern std::vector<std::vector<double>> cespm, cen, potpc, co, potpt, ovl,
    cespm2, a2, exs, ce, expn, ewcx, ewcy, ewcz, f0, f1, u_esp, rnai, rnai1,
    rnai2, espi, exsr;
extern std::vector<std::vector<bool>> cequiv;
}  // namespace esp_C
