// local_for_MOZYME.cpp — C++ translation of MOPAC 2016 "local_for_MOZYME.F90".
// Occupied/virtual LMO localisation (maximise sum of (psi)^4 by pairwise rotation).
#include "local_for_MOZYME.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "chanel_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
using namespace MOZYME_C;

namespace molkst_C {
extern int norbs, numat, nelecs, natoms;
}
namespace MOZYME_C {
}
namespace chanel_C {
extern int iw;
}

extern void memory_error(const std::string&);
extern void mopend(const std::string&);

// Global MOZYME linkage.

// localize_for_MOZYME: rotate pairs of LMOs to maximise (PSI)^4.
static void localize_for_MOZYME(
    std::vector<double>& c, int n289, const std::vector<int>& ic, int n267,
    const std::vector<int>& nc, const std::vector<int>& ncstrt,
    int nmos_loc, const std::vector<int>& iorbs,
    std::vector<double>& psi1, std::vector<double>& psi2,
    std::vector<double>& axiiii, std::vector<int>& nf, std::vector<int>& nl,
    std::vector<std::vector<int>>& ioc, const std::vector<int>& nnc_loc,
    double& totij) {
    int i = molkst_C::natoms;  // dummy use
    (void)i;
    // Build all the DII.
    int l = 0;
    for (int k = 1; k <= nmos_loc; ++k) {
        l = nnc_loc[k];
        int m = ncstrt[k];
        axiiii[k] = 0.0;
        for (int ii = 1; ii <= nc[k]; ++ii) {
            l = l + 1;
            int i1 = ic[l];
            double dii = 0.0;
            for (int j = 1; j <= iorbs[i1]; ++j) {
                m = m + 1;
                dii = dii + c[m] * c[m];
            }
            axiiii[k] = axiiii[k] + dii * dii;
        }
    }
    double total = 0.0;
    totij = 0.0;
    for (int iloop = 1; iloop <= nmos_loc; ++iloop) {
        int i5 = ncstrt[iloop];
        int i8 = nnc_loc[iloop];
        double xiiii = axiiii[iloop];
        for (int jloop = 1; jloop <= nmos_loc; ++jloop) {
            int j8 = nnc_loc[jloop];
            int j5 = ncstrt[jloop];
            if (jloop == iloop) continue;
            bool common_atom = false;
            for (int ia = 1; ia <= 2; ++ia)
                for (int jb = 1; jb <= 2; ++jb)
                    if (ic[ia + i8] == ic[jb + j8]) common_atom = true;
            if (!common_atom) continue;
            int ij = 0, ijorb = 0, il = 0;
            for (int ia = 1; ia <= nc[iloop]; ++ia) {
                int i1 = ic[ia + i8];
                int jl = 0;
                for (int jb = 1; jb <= nc[jloop]; ++jb) {
                    int j1 = ic[jb + j8];
                    if (i1 == j1) {
                        ij = ij + 1;
                        nf[ij] = ijorb + 1;
                        nl[ij] = ijorb + iorbs[i1];
                        ioc[1][ij] = il;
                        ioc[2][ij] = jl;
                        int jl1 = jl, il1 = il;
                        for (int k = 1; k <= iorbs[i1]; ++k) {
                            ijorb = ijorb + 1;
                            il1 = il1 + 1;
                            jl1 = jl1 + 1;
                            psi1[ijorb] = c[il1 + i5];
                            psi2[ijorb] = c[jl1 + j5];
                        }
                    }
                    jl = jl + iorbs[j1];
                }
                il = il + iorbs[i1];
            }
            double xijjj = 0.0, xjiii = 0.0, xijij = 0.0, xiijj = 0.0;
            for (int k1 = 1; k1 <= ij; ++k1) {
                double dij = 0.0, dii = 0.0, djj = 0.0;
                for (int k = nf[k1]; k <= nl[k1]; ++k) {
                    dij = dij + psi1[k] * psi2[k];
                    dii = dii + psi1[k] * psi1[k];
                    djj = djj + psi2[k] * psi2[k];
                }
                xijjj = xijjj + dij * djj;
                xjiii = xjiii + dij * dii;
                xijij = xijij + dij * dij;
                xiijj = xiijj + dii * djj;
            }
            if (xiijj >= 0.001) {
                double xjjjj = axiiii[jloop];
                double aij = xijij - (xiiii + xjjjj - 2.0 * xiijj) / 4.0;
                double bij = xjiii - xijjj;
                double ca = std::sqrt(aij * aij + bij * bij);
                double sa = aij + ca;
                if (sa > 1.e-14) {
                    ca = (1.0 + std::sqrt((1.0 - aij / ca) / 2.0)) / 2.0;
                    sa = std::sqrt(1.0 - ca);
                    ca = std::sqrt(ca);
                    totij = totij + sa;
                    int ii = 0;
                    for (int k = 1; k <= ij; ++k) {
                        int il2 = 0;
                        for (int i2 = nf[k]; i2 <= nl[k]; ++i2) {
                            il2 = il2 + 1;
                            ii = ii + 1;
                            c[ioc[1][k] + il2 + i5] = ca * psi1[ii] + sa * psi2[ii];
                            c[ioc[2][k] + il2 + j5] = -sa * psi1[ii] + ca * psi2[ii];
                        }
                    }
                }
            }
        }
        total = total + xiiii;
    }
}

void local_for_MOZYME(const char* type) {
    int nocc = molkst_C::nelecs / 2;
    int nvir = molkst_C::norbs - nocc;
    std::vector<double> psi1(molkst_C::norbs + 1), psi2(molkst_C::norbs + 1),
        axiiii(molkst_C::norbs + 1);
    std::vector<int> nf(molkst_C::numat + 1), nl(molkst_C::numat + 1);
    std::vector<std::vector<int>> ioc(3, std::vector<int>(molkst_C::numat + 1));
    std::string t(type);
    double totij = 0.0;
    if (t == "OCCUPIED") {
        for (int i = 1; i <= 100; ++i) {
            localize_for_MOZYME(cocc, cocc_dim, icocc, icocc_dim, ncf, ncocc,
                                nocc, MOZYME_C::iorbs, psi1, psi2, axiiii, nf, nl, ioc, nncf, totij);
            if (totij < 1.e-10) break;
        }
    } else if (t == "VIRTUAL") {
        for (int i = 1; i <= 100; ++i) {
            localize_for_MOZYME(cvir, cvir_dim, icvir, icvir_dim, nce, ncvir,
                                nvir, MOZYME_C::iorbs, psi1, psi2, axiiii, nf, nl, ioc, nnce, totij);
            if (totij < 1.e-10) break;
        }
    } else {
        std::fprintf(stdout, " Error\n");
        mopend("Error in LOCAL");
    }
}
