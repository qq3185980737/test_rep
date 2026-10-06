// overlaps_C.h — C++ mapping of Fortran module "overlaps_C".
#pragma once
namespace overlaps_C {
constexpr int N_AO = 60;
// STO-6G coefficients and exponents, dimension(60,6), 1-based.
extern double ccc[N_AO + 1][7];
extern double zzz[N_AO + 1][7];
// STO-6G master tables (setupg.F90): allz/allc(6,6,2), 1-based [nqn][i][l].
extern double allz[7][7][3];
extern double allc[7][7][3];
// fact(0..17): factorials n!
extern double fact[18];
extern double cutof1, cutof2;
}
