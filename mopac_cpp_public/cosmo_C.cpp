// cosmo_C.cpp — C++ translation of MOPAC 2016 "cosmo_C.F90".
#include "cosmo_C.h"

namespace cosmo_C {
bool iseps = false, noeps = false, useps = false, lpka = false;
int nspa = 0, nps = 0, nps2 = 0, nden = 0, lenabc = 0, nppa = 1082;
int amat_dim = 0, isude_dim = 0, nipc = 0, ioldcv = 0;
int n0[3] = {0, 0, 0};
std::vector<int> iatsp, nar_csm, nsetf, cosmo_i, ipiden, idenat, nset;
std::vector<std::vector<int>> isude, nn;
double fepsi = 0.0, rds = 0.0, disex2 = 0.0, ediel = 0.0, solv_energy = 0.0;
double area = 0.0, cosvol = 0.0;
double cif1 = 0.0, cif2 = 0.0, cdiagi = 0.0, fnsq = 0.0, rsolv = 0.0;
double tm[4][4][1083] = {};
double dirsm[5][1083] = {}, dirvec[5][1083] = {};
std::vector<double> amat, cmat, gden, qscat, arat, srad, abcmat, qden, bh, cdiag;
std::vector<std::vector<double>> bmat, phinet, qscnet, qdenet, cosurf, sude, xsp, cxy;
}
