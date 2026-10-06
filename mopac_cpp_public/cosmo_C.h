// cosmo_C.h — C++ translation of MOPAC 2016 "cosmo_C.F90".
#pragma once
#include <vector>

namespace cosmo_C {
extern bool iseps, noeps, useps, lpka;
extern int nspa, nps, nps2, nden, lenabc, nppa, amat_dim, isude_dim, nipc, ioldcv;
extern int n0[3];
extern std::vector<int> iatsp, nar_csm, nsetf, cosmo_i, ipiden, idenat, nset;
extern std::vector<std::vector<int>> isude, nn;
extern double fepsi, rds, disex2, ediel, solv_energy, area, cosvol;
extern double cif1, cif2, cdiagi, fnsq, rsolv;
extern double tm[4][4][1083];
extern double dirsm[5][1083], dirvec[5][1083];
extern std::vector<double> amat, cmat, gden, qscat, arat, srad, abcmat, qden, bh, cdiag;
extern std::vector<std::vector<double>> bmat, phinet, qscnet, qdenet, cosurf, sude, xsp, cxy;
}
