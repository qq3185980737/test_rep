// dfock2.cpp — C++ translation of MOPAC 2016 "dfock2.F90".
// NDDO 2-electron 2-center repulsion derivative of the Fock matrix.
// ifact/i1fact/ptot2 live in common_arrays_C (shared with setup_mopac_arrays),
// matching F90.  All packed arrays 1-based.
#include "dfock2.h"
#include <algorithm>
#include <cmath>
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "jab.h"
#include "kab.h"
using namespace common_arrays_C;
using namespace molkst_C;

namespace {
// fockdorbs: general d-orbital two-center derivative block (same as fock2.F90).
void fockdorbs(int ia, int ib, int ja, int jb, double* f,
               const double* p, const double* ptot, const double* w,
               int& kr, const std::vector<int>& fact) {
    if (ia > ja) {
        for (int i = ia; i <= ib; ++i) {
            int ka = fact[i];
            double aa = 2.0;
            for (int j = ia; j <= i; ++j) {
                if (i == j) aa = 1.0;
                int kb = fact[j];
                int ij = ka + j;
                for (int k = ja; k <= jb; ++k) {
                    int kc = fact[k];
                    int ik = ka + k;
                    int jk = kb + k;
                    double bb = 2.0;
                    for (int l = ja; l <= k; ++l) {
                        if (k == l) bb = 1.0;
                        int il = ka + l;
                        int jl = kb + l;
                        int kl = kc + l;
                        kr = kr + 1;
                        double a = w[kr];
                        f[ij] += bb * a * ptot[kl];
                        f[kl] += aa * a * ptot[ij];
                        a = a * aa * bb * 0.25;
                        f[ik] -= a * p[jl];
                        f[il] -= a * p[jk];
                        f[jk] -= a * p[il];
                        f[jl] -= a * p[ik];
                    }
                }
            }
        }
    } else {
        int kref = kr;
        int nn = jb - ja + 1;
        nn = (nn * (nn + 1)) / 2;
        int n1 = 0;
        for (int i = ja; i <= jb; ++i) {
            int ka = fact[i];
            double aa = 2.0;
            for (int j = ja; j <= i; ++j) {
                n1 = n1 + 1;
                if (i == j) aa = 1.0;
                int kb = fact[j];
                int ij = ka + j;
                int n2 = 0;
                for (int k = ia; k <= ib; ++k) {
                    int kc = fact[k];
                    int ik = ka + k;
                    int jk = kb + k;
                    double bb = 2.0;
                    for (int l = ia; l <= k; ++l) {
                        n2 = n2 + 1;
                        if (k == l) bb = 1.0;
                        int il = ka + l;
                        int jl = kb + l;
                        int kl = kc + l;
                        kr = kr + 1;
                        double a = w[kref + (n2 - 1) * nn + n1];
                        f[ij] += bb * a * ptot[kl];
                        f[kl] += aa * a * ptot[ij];
                        a = a * aa * bb * 0.25;
                        f[ik] -= a * p[jl];
                        f[il] -= a * p[jk];
                        f[jk] -= a * p[il];
                        f[jl] -= a * p[ik];
                    }
                }
            }
        }
    }
}

int jindex[257];
int itype = 1;
int icalcn = 0;

} // namespace

void dfock2(double* f, const double* ptot, const double* p, const double* w,
            int numat, const int* nfirst, const int* nlast, int nati) {
    if (icalcn != molkst_C::numcal) {
        // jindex build needs ifact[1..4]; keep at least 5 slots even for
        // tiny test systems (F90 assumes norbs >= 4 in real molecules)
        int sz = std::max(molkst_C::norbs + 1, 5);
        common_arrays_C::ifact.assign(sz, 0);
        common_arrays_C::i1fact.assign(sz, 0);
        common_arrays_C::ptot2.assign(numat * 82, 0.0);
        icalcn = molkst_C::numcal;
        itype = 1;
    }
restart:
    switch (itype) {
    default:
        for (int i = 1; i <= molkst_C::norbs; ++i) {
            common_arrays_C::ifact[i] = (i * (i - 1)) / 2;
            common_arrays_C::i1fact[i] = common_arrays_C::ifact[i] + i;
        }
        {
            int m = 0;
            for (int i = 1; i <= 4; ++i)
                for (int j = 1; j <= 4; ++j) {
                    int ij = std::min(i, j);
                    int ji = i + j - ij;
                    for (int k = 1; k <= 4; ++k) {
                        int ik = std::min(i, k);
                        for (int l = 1; l <= 4; ++l) {
                            m = m + 1;
                            int kl = std::min(k, l);
                            int lk = k + l - kl;
                            jindex[m] = (common_arrays_C::ifact[ji] + ij) * 10 +
                                        common_arrays_C::ifact[lk] + kl - 10;
                        }
                    }
                }
        }
        itype = 3;
        goto restart;
    case 3:
        break;
    case 2:
        // MNDO (d-orbital-free) fast path: one elrep integral per atom pair
        {
            int kr = 0;
            int ii = nati;
            int ia = nfirst[ii], ib = nlast[ii];
            for (int jj = 1; jj <= numat; ++jj) {
                if (jj == ii) continue;
                kr = kr + 1;
                double elrep = w[kr];
                int ja = nfirst[jj], jb = nlast[jj];
                if (ja < ia) {
                    for (int i = ia; i <= ib; ++i) {
                        int ka = common_arrays_C::ifact[i];
                        int kk = ka + i;
                        for (int k = ja; k <= jb; ++k) {
                            int ll = common_arrays_C::i1fact[k];
                            int ik = ka + k;
                            f[kk] += ptot[ll] * elrep;
                            f[ll] += ptot[kk] * elrep;
                            f[ik] -= p[ik] * elrep;
                        }
                    }
                } else {
                    for (int i = ia; i <= ib; ++i) {
                        int ka = common_arrays_C::ifact[i];
                        int kk = ka + i;
                        for (int k = ja; k <= jb; ++k) {
                            int ll = common_arrays_C::i1fact[k];
                            int ik = ll + i - k;
                            f[kk] += ptot[ll] * elrep;
                            f[ll] += ptot[kk] * elrep;
                            f[ik] -= p[ik] * elrep;
                        }
                    }
                }
            }
        }
        return;
    }
    // case 3: main NDDO path
    {
        int kk = 0;
        for (int i = 1; i <= numat; ++i) {
            int ia = nfirst[i];
            int ib = nlast[i];
            int m = 0;
            for (int j = ia; j <= ib; ++j)
                for (int k = ia; k <= ib; ++k) {
                    m = m + 1;
                    int jk = std::min(j, k);
                    int kj = k + j - jk;
                    jk = jk + (kj * (kj - 1)) / 2;
                    common_arrays_C::ptot2[(i - 1) * 82 + (m - 1)] = ptot[jk];
                }
        }
        int ii = nati;
        int ia = nfirst[ii], ib = nlast[ii];
        for (int jj = 1; jj <= numat; ++jj) {
            if (ii == jj) continue;
            int ja = nfirst[jj], jb = nlast[jj];
            if (ib - ia < 0 || jb - ja < 0) continue;  // sparkle atom
            if (ib - ia >= 6 || jb - ja >= 6) {
                fockdorbs(ia, ib, ja, jb, f, p, ptot, w, kk,
                          common_arrays_C::ifact);
            } else if (ib - ia >= 3 && jb - ja >= 3) {
                // heavy - heavy: jab + kab, w block 100
                double pja[17] = {}, pjb[17] = {}, pk[17] = {};
                for (int mm = 1; mm <= 16; ++mm) {
                    pja[mm] = common_arrays_C::ptot2[(ii - 1) * 82 + (mm - 1)];
                    pjb[mm] = common_arrays_C::ptot2[(jj - 1) * 82 + (mm - 1)];
                }
                jab(ia, ja, pja, pjb, &w[kk + 1], f);
                int l = 0;
                if (ia > ja) {
                    for (int i = ia; i <= ib; ++i) {
                        for (int q = 0; q < jb - ja + 1; ++q)
                            pk[l + 1 + q] = p[common_arrays_C::ifact[i] + ja + q];
                        l = jb - ja + 1 + l;
                    }
                } else {
                    for (int i = ia; i <= ib; ++i) {
                        for (int q = 0; q < jb - ja + 1; ++q)
                            pk[l + 1 + q] = p[common_arrays_C::ifact[ja + q] + i];
                        l = jb - ja + 1 + l;
                    }
                }
                kab(ia, ja, pk, &w[kk + 1], f);
                kk = kk + 100;
            } else if (ib - ia >= 3) {
                // light-atom jj - heavy-atom ii (Coulomb + exchange)
                if (nati == 1) { fprintf(stderr, "[DF2] L-H jj=%d kk=%d jindex1..16:", jj, kk); for (int q = 1; q <= 16; ++q) fprintf(stderr, " %d", jindex[q]); fprintf(stderr, "\n"); fflush(stderr); }
                double sumdia = 0.0, sumoff = 0.0;
                int ll = common_arrays_C::i1fact[ja];
                int k = 0;
                for (int i = 0; i <= 3; ++i) {
                    int j1 = common_arrays_C::ifact[ia + i] + ia - 1;
                    for (int j = 0; j <= i - 1; ++j) {
                        k = k + 1;
                        j1 = j1 + 1;
                        f[j1] += ptot[ll] * w[kk + k];
                        sumoff += ptot[j1] * w[kk + k];
                    }
                    j1 = j1 + 1;
                    k = k + 1;
                    f[j1] += ptot[ll] * w[kk + k];
                    sumdia += ptot[j1] * w[kk + k];
                }
                f[ll] += sumoff * 2.0 + sumdia;
                int i1, j1;
                if (ia > ja) {
                    k = 0;
                    for (int i = ia; i <= ib; ++i) {
                        i1 = common_arrays_C::ifact[i] + ja;
                        double sum = 0.0;
                        for (int j = ia; j <= ib; ++j) {
                            k = k + 1;
                            j1 = common_arrays_C::ifact[j] + ja;
                            sum += p[j1] * w[kk + jindex[k]];
                        }
                        f[i1] -= sum;
                    }
                } else {
                    k = 0;
                    for (int i = ia; i <= ib; ++i) {
                        i1 = common_arrays_C::ifact[ja] + i;
                        double sum = 0.0;
                        for (int j = ia; j <= ib; ++j) {
                            k = k + 1;
                            j1 = common_arrays_C::ifact[ja] + j;
                            sum += p[j1] * w[kk + jindex[k]];
                            if (nati == 1 && jj == 2 && i == 2)
                                fprintf(stderr, "[DF2] EX i2 j=%d k=%d j1=%d p=%+.6f ji=%d w=%+.6f term=%+.6f\n",
                                        j, k, j1, p[j1], jindex[k], w[kk + jindex[k]], p[j1] * w[kk + jindex[k]]);
                        }
                        f[i1] -= sum;
                        if (nati == 1 && jj == 2 && i == 2)
                            fprintf(stderr, "[DF2] EX i2 i1=%d sum=%+.6f f[i1]-=%+.6f\n", i1, sum, f[i1]);
                    }
                }
                if (nati == 1 && jj == 2) {
                    fprintf(stderr, "[DF2] COUL ll=%d w1..10:", ll);
                    for (int q = 1; q <= 10; ++q) fprintf(stderr, " %+.6f", w[kk + q]);
                    fprintf(stderr, " ptot15=%+.6f\n", ptot[ll]); fflush(stderr);
                }
                kk = kk + 10;
            } else if (jb - ja >= 3) {
                // heavy-atom jj - light-atom ii
                double sumdia = 0.0, sumoff = 0.0;
                int ll = common_arrays_C::i1fact[ia];
                int k = 0;
                for (int i = 0; i <= 3; ++i) {
                    int j1 = common_arrays_C::ifact[ja + i] + ja - 1;
                    for (int j = 0; j <= i - 1; ++j) {
                        k = k + 1;
                        j1 = j1 + 1;
                        f[j1] += ptot[ll] * w[kk + k];
                        sumoff += ptot[j1] * w[kk + k];
                    }
                    j1 = j1 + 1;
                    k = k + 1;
                    f[j1] += ptot[ll] * w[kk + k];
                    sumdia += ptot[j1] * w[kk + k];
                }
                f[ll] += sumoff * 2.0 + sumdia;
                if (ia > ja) {
                    k = common_arrays_C::ifact[ia] + ja;
                    int j = 0;
                    for (int i = k; i <= k + 3; ++i) {
                        double sum = 0.0;
                        for (int l = k; l <= k + 3; ++l) {
                            j = j + 1;
                            sum += p[l] * w[kk + jindex[j]];
                        }
                        f[i] -= sum;
                    }
                } else {
                    int j = 0;
                    for (int kk2 = ja; kk2 <= ja + 3; ++kk2) {
                        int i = common_arrays_C::ifact[kk2] + ia;
                        double sum = 0.0;
                        for (int ll2 = ja; ll2 <= ja + 3; ++ll2) {
                            int l = common_arrays_C::ifact[ll2] + ia;
                            j = j + 1;
                            sum += p[l] * w[kk + jindex[j]];
                        }
                        f[i] -= sum;
                    }
                }
                kk = kk + 10;
            } else {
                // light - light
                int i1 = common_arrays_C::i1fact[ia];
                int j1 = common_arrays_C::i1fact[ja];
                f[i1] += ptot[j1] * w[kk + 1];
                f[j1] += ptot[i1] * w[kk + 1];
                int ij = (ia > ja) ? (i1 + ja - ia) : (j1 + ia - ja);
                f[ij] -= p[ij] * w[kk + 1];
                kk = kk + 1;
            }
        }
    }
    return;
}
