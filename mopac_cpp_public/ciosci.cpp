// ciosci.cpp — C++ translation of MOPAC 2016 "ciosci.F90".
// matout is an external stub.

#include "ciosci.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace funcon_C;
using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;
using namespace parameters_C;

namespace {
void matout(const std::vector<std::vector<double>>&,
            const std::vector<std::vector<double>>&, int, int, int) {}

const int nspqn[108] = {
    0,
    1,1, 2,2,2,2,2,2,2,2,        // 2*1, 8*2
    3,3,3,3,3,3,3,3,             // 8*3
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,  // 18*4
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,  // 18*5
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,  // 32*6
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  // 21*0
};
}

void ciosci(const std::vector<std::vector<double>>& vects, int lroot,
            std::vector<std::vector<double>>& oscil,
            const std::vector<double>& conf_in) {
    bool debug = keywrd.find("CIOSCI") != std::string::npos;
    (void)debug;

    // factorial table up to 14.
    std::vector<int> ll(15, 1);
    ll[1] = 1;
    for (int i = 2; i <= 14; ++i) ll[i] = ll[i - 1] * i;

    std::vector<double> sppol(numat + 1, 0.0), pdpol(numat + 1, 0.0);
    for (int i = 1; i <= numat; ++i) sppol[i] = a0 * dd[nat[i]];

    for (int i = 1; i <= numat; ++i) {
        if (nlast[i] - nfirst[i] > 5) {
            int j = nat[i];
            int nd = nspqn[j];
            int npp = nd;
            if (j > 20) npp--;
            pdpol[i] = a0 * ll[npp + nd + 1] * std::pow(2.0, npp + nd + 1) *
                       std::pow(zp[j], npp + 0.5) * std::pow(zd[j], nd + 0.5) /
                       (std::sqrt(5.0) * std::pow(zp[j] + zd[j], npp + nd + 2) *
                        std::sqrt((double)(ll[2 * npp] * ll[2 * nd])));
        }
    }

    oscil.assign(4, std::vector<double>(lab + 1, 0.0));

    // (norbs, nmos) working arrays, 1-based.
    std::vector<std::vector<double>> vect1(norbs + 1, std::vector<double>(nmos + 1, 0.0));
    std::vector<std::vector<double>> vect2(norbs + 1, std::vector<double>(nmos + 1, 0.0));
    std::vector<std::vector<double>> t2(nmos + 1, std::vector<double>(nmos + 1, 0.0));
    std::vector<std::vector<double>> t4(lab + 1, std::vector<double>(lab + 1, 0.0));
    std::vector<double> work(lab + 1, 0.0);
    std::vector<int> lll(15, 0);

    for (int loop = 1; loop <= 3; ++loop) {
        for (int iloop = 1; iloop <= nmos; ++iloop) {
            for (int iatom = 1; iatom <= numat; ++iatom) {
                for (int i = nfirst[iatom]; i <= nlast[iatom]; ++i) {
                    vect1[i][iloop] = vects[i][iloop];
                    vect2[i][iloop] = vects[i][iloop] * coord[loop-1][iatom];
                    if (nat[iatom] != 1) {
                        if (i == nfirst[iatom]) {
                            vect2[i][iloop] += sppol[iatom] * vects[i + loop][iloop];
                        } else if (i - nfirst[iatom] == loop) {
                            vect2[i][iloop] += sppol[iatom] * vects[i - loop][iloop];
                        }
                        if (i == nfirst[iatom] + 8) {
                            int k = nfirst[iatom], l = k + 3;
                            double oneos3 = 1.0 / std::sqrt(3.0);
                            if (loop == 1) {
                                vect2[k+1][iloop] += pdpol[iatom] * (vects[l+1][iloop] - oneos3 * vects[l+3][iloop]);
                                vect2[k+2][iloop] += pdpol[iatom] * vects[l+5][iloop];
                                vect2[k+3][iloop] += pdpol[iatom] * vects[l+2][iloop];
                                vect2[l+1][iloop] += pdpol[iatom] * vects[k+1][iloop];
                                vect2[l+2][iloop] += pdpol[iatom] * vects[k+3][iloop];
                                vect2[l+3][iloop] -= oneos3 * pdpol[iatom] * vects[k+1][iloop];
                                vect2[l+5][iloop] += pdpol[iatom] * vects[k+2][iloop];
                            } else if (loop == 2) {
                                vect2[k+1][iloop] += pdpol[iatom] * vects[l+5][iloop];
                                vect2[k+2][iloop] += pdpol[iatom] * (-vects[l+1][iloop] - oneos3 * vects[l+3][iloop]);
                                vect2[k+3][iloop] += pdpol[iatom] * vects[l+4][iloop];
                                vect2[l+1][iloop] -= pdpol[iatom] * vects[k+2][iloop];
                                vect2[l+3][iloop] -= oneos3 * pdpol[iatom] * vects[k+2][iloop];
                                vect2[l+4][iloop] += pdpol[iatom] * vects[k+3][iloop];
                                vect2[l+5][iloop] += pdpol[iatom] * vects[k+1][iloop];
                            } else {
                                vect2[k+1][iloop] += pdpol[iatom] * vects[l+2][iloop];
                                vect2[k+2][iloop] += pdpol[iatom] * vects[l+4][iloop];
                                vect2[k+3][iloop] += pdpol[iatom] * 2.0 / std::sqrt(3.0) * vects[l+3][iloop];
                                vect2[l+2][iloop] += pdpol[iatom] * vects[k+1][iloop];
                                vect2[l+3][iloop] += 2.0 / std::sqrt(3.0) * pdpol[iatom] * vects[k+3][iloop];
                                vect2[l+4][iloop] += pdpol[iatom] * vects[k+2][iloop];
                            }
                        }
                    }
                }
            }
        }

        for (int i = 1; i <= nmos; ++i) {
            t2[i][i] = 0.0;
            for (int j = 1; j < i; ++j) {
                double sum = 0.0;
                for (int k = 1; k <= norbs; ++k) sum += vect1[k][i] * vect2[k][j];
                t2[i][j] = sum;
                t2[j][i] = -sum;
            }
        }

        for (int i = 1; i <= lab; ++i) {
            for (int j = 1; j <= i; ++j) {
                int lc = 0, mc = 0;
                for (int k = 1; k <= nmos; ++k) {
                    if (microa[k][i] != microa[k][j]) { ++lc; lll[lc] = k; }
                    if (microb[k][i] != microb[k][j]) { ++mc; lll[mc] = k; }
                }
                if ((lc == 2 && mc == 0) || (mc == 2 && lc == 0)) {
                    t4[i][j] = t2[lll[1]][lll[2]];
                    int ij;
                    if (lc == 2) {
                        int ii2 = lll[1];
                        ij = microb[ii2][i] + microa[ii2][i];
                        for (int jj = ii2 + 1; jj <= lll[2] - 1; ++jj)
                            ij += microa[jj][i] + microb[jj][i];
                    } else {
                        int ii2 = lll[1];
                        ij = microb[ii2][i];
                        for (int jj = ii2 + 1; jj <= lll[2] - 1; ++jj)
                            ij += microa[jj][i] + microb[jj][i];
                        ij += microa[lll[2]][i];
                    }
                    if (ij % 2 == 1) t4[i][j] = -t4[i][j];
                    t4[j][i] = -t4[i][j];
                } else {
                    t4[i][j] = 0.0;
                    t4[j][i] = 0.0;
                }
            }
        }

        for (int ns = 1; ns <= nstate; ++ns) {
            int ii = (lroot + ns - 2) * lab;
            for (int istate = 1; istate <= lab; ++istate) {
                for (int j = 1; j <= lab; ++j) {
                    double sum = 0.0;
                    for (int k = 1; k <= lab; ++k) sum += conf_in[k + ii] * t4[j][k];
                    work[j] = sum;
                }
                double sum = 0.0;
                for (int k = 1; k <= lab; ++k) sum += work[k] * conf_in[k + (istate - 1) * lab];
                oscil[loop][istate] += sum * sum;
            }
        }
    }
}
