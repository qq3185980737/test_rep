// fillij.cpp — C++ translation of MOPAC 2016 "fillij.F90".
#include "fillij.h"
#include <cmath>
#include <algorithm>
#include <string>
#include "MOZYME_C.h"
#include "molkst_C.h"
#include "overlaps_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace overlaps_C;

extern double reada(const std::string&, int);

void fillij(bool count) {
    static int ix = 0;
    bool first;
    if (MOZYME_C::nijbo.empty()) {
        ix = 0;
        n2elec = 0;
        mpack = 0;
        first = true;
        MOZYME_C::morb = 4;
        if (ispd > 0) MOZYME_C::morb = 9;
    } else {
        first = false;
    }
    if (!count) {
        if (MOZYME_C::nijbo.empty()) {
            if (MOZYME_C::lijbo) {
                MOZYME_C::nijbo.assign(numat, std::vector<int>(numat, 0));
            }
            if (!MOZYME_C::lijbo) {
                MOZYME_C::iijj.assign(MOZYME_C::ij_dim, 0);
                MOZYME_C::ijall.assign(MOZYME_C::ij_dim, 0);
                MOZYME_C::iij.assign(natoms + 1, 0);
                MOZYME_C::numij.assign(natoms + 1, 0);
            }
        }
    } else {
        MOZYME_C::lijbo = true;
    }
    cutof1 = 10.0 * 10.0;
    if (cutofp < 100.0) cutof1 = cutofp * cutofp + 1e-4;
    if (numat < 30) cutof1 = 1e6 + 10;
    cutof2 = 9.9 * 9.9;
    if (cutofp < 100.0) cutof2 = cutofp * cutofp;
    if (numat < 30) cutof2 = 1e6;
    size_t pos;
    pos = keywrd.find(" CUTOFF=");
    if (pos != std::string::npos) {
        int ii = (int)pos + 8;
        double v = reada(keywrd, ii);
        cutof1 = v * v;
        cutof2 = cutof1 - 1e-1;
    } else {
        pos = keywrd.find(" CUTOF1=");
        if (pos != std::string::npos) {
            int ii = (int)pos + 8;
            double v = reada(keywrd, ii);
            cutof1 = v * v;
        }
        pos = keywrd.find(" CUTOF2=");
        if (pos != std::string::npos) {
            int ii = (int)pos + 8;
            double v = reada(keywrd, ii);
            cutof2 = v * v;
        }
    }
    pos = keywrd.find(" CUTOFS=");
    if (pos != std::string::npos) {
        int ii = (int)pos + 8;
        double v = reada(keywrd, ii);
        MOZYME_C::cutofs = v * v;
    } else {
        MOZYME_C::cutofs = 7.0 * 7.0;
    }
    if (cutof1 < cutof2) cutof1 = cutof2 + 1e-4;
    if (keywrd.find(" NODIRECT") != std::string::npos || id != 0) {
        MOZYME_C::direct = false;
        MOZYME_C::semidr = false;
    } else if (keywrd.find(" SEMIDIRECT") != std::string::npos) {
        MOZYME_C::direct = false;
        MOZYME_C::semidr = true;
    } else {
        MOZYME_C::direct = true;
        MOZYME_C::semidr = true;
    }
    int i = 0;
    for (int iloop = 1; iloop <= numat; ++iloop) {
        int io = MOZYME_C::iorbs[iloop];
        i = i + 1;
        double x1 = common_arrays_C::coord[0][i];
        double x2 = common_arrays_C::coord[1][i];
        double x3 = common_arrays_C::coord[2][i];
        int ii = (io * (io + 1)) / 2;
        if (first) n2elec = n2elec + ii * ii;
        if (!MOZYME_C::lijbo && !count) MOZYME_C::iij[i] = ix + 1;
        int j = 0;
        for (int jloop = 1; jloop <= iloop - 1; ++jloop) {
            int jo = MOZYME_C::iorbs[jloop];
            j = j + 1;
            double r;
            if (id == 0) {
                r = (x1 - common_arrays_C::coord[0][j]) * (x1 - common_arrays_C::coord[0][j]) +
                    (x2 - common_arrays_C::coord[1][j]) * (x2 - common_arrays_C::coord[1][j]) +
                    (x3 - common_arrays_C::coord[2][j]) * (x3 - common_arrays_C::coord[2][j]);
            } else {
                r = 1e10;
                double xj[3];
                for (int ip = -l1u; ip <= l1u; ++ip)
                  for (int jp = -l2u; jp <= l2u; ++jp)
                    for (int kp = -l3u; kp <= l3u; ++kp) {
                      for (int lp = 0; lp < 3; ++lp)
                        xj[lp] = common_arrays_C::coord[lp][j] +
                          common_arrays_C::tvec[lp][0] * ip +
                          common_arrays_C::tvec[lp][1] * jp +
                          common_arrays_C::tvec[lp][2] * kp;
                      double rr = (x1 - xj[0]) * (x1 - xj[0]) +
                                  (x2 - xj[1]) * (x2 - xj[1]) +
                                  (x3 - xj[2]) * (x3 - xj[2]);
                      r = std::min(r, rr);
                    }
            }
            if (r < cutof2) {
                ix = ix + 1;
                if (!count) {
                    if (MOZYME_C::lijbo) {
                        if (first || MOZYME_C::nijbo[i-1][j-1] < 0) {
                            MOZYME_C::nijbo[i-1][j-1] = mpack;
                            MOZYME_C::nijbo[j-1][i-1] = mpack;
                        } else {
                            mpack = mpack - io * jo;
                        }
                    } else {
                        MOZYME_C::iijj[ix] = mpack;
                        MOZYME_C::ijall[ix] = j;
                    }
                }
                mpack = mpack + io * jo;
                if (!MOZYME_C::direct) n2elec = n2elec + (jo * (jo + 1)) / 2 * ii;
            } else if (r < cutof1) {
                if (!MOZYME_C::semidr) {
                    if (io > 1) n2elec = n2elec + (jo > 1 ? 7 : 4);
                    else n2elec = n2elec + (jo > 1 ? 4 : 1);
                }
                if (MOZYME_C::lijbo && !count) {
                    if (first || MOZYME_C::nijbo[i-1][j-1] == -1) {
                        MOZYME_C::nijbo[i-1][j-1] = -2;
                        MOZYME_C::nijbo[j-1][i-1] = -2;
                    }
                }
            } else {
                if (!MOZYME_C::semidr) n2elec = n2elec + 1;
                if (MOZYME_C::lijbo && !count && first) {
                    MOZYME_C::nijbo[i-1][j-1] = -1;
                    MOZYME_C::nijbo[j-1][i-1] = -1;
                }
            }
        }
        ix = ix + 1;
        if (!count) {
            if (MOZYME_C::lijbo) {
                if (first) {
                    MOZYME_C::nijbo[i-1][i-1] = mpack;
                } else {
                    mpack = mpack - (io * (io + 1)) / 2;
                }
            } else {
                MOZYME_C::iijj[ix] = mpack;
                MOZYME_C::ijall[ix] = i;
                MOZYME_C::numij[i] = ix;
            }
        }
        if (id != 0) n2elec = n2elec + ((io * (io + 1)) / 2) * ((io * (io + 1)) / 2);
        mpack = mpack + (io * (io + 1)) / 2;
    }
    if (count) MOZYME_C::ij_dim = ix;
    if (n2elec < 2025) n2elec = 2025;
    if (first) {
        n2elec = n2elec + 10;
        if (MOZYME_C::direct && ispd == 0) n2elec = n2elec + 100;
        if (MOZYME_C::direct && ispd != 0) n2elec = n2elec + 2025;
    }
    if (id != 0) n2elec = n2elec * 2;
}
