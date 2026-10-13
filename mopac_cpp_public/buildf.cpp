// buildf.cpp — C++ translation of MOPAC 2016 "buildf.F90"
//   + fock2z.F90 (fz2n fast two-center path + focd2z d-orbital path).
// 1-based semantics throughout: f, ptot, w, wk are 1-based vectors of size
// mpack+1; nijbo[i][j] holds the 0-based block start (as Fortran ijbo returns).

#include "buildf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <omp.h>

#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace MOZYME_C;
using namespace parameters_C;
using namespace funcon_C;
using namespace cosmo_C;

extern void chrge_for_MOZYME(const double* p, double* q, int numat, const int* iorbs);
extern void rotate(int ni, int nj, const double* xi, const double* xj, double* w,
                   int& kr, double* e1b, double* e2a, double& enuc);
extern void jab_for_MOZYME(int ia, int ja, const double* pja, const double* pjb,
                           const double* w, double* f1, double* f2);
extern void kab_for_MOZYME(int ia, int ja, const double* pk, const double* w, double* f);
extern void fock1_for_MOZYME(std::vector<double>& f, const std::vector<double>& ptot,
                             const std::vector<std::vector<double>>& w,
                             int& kr, int iab, int ilim);
extern void to_point(double rij, double& point, double& cnst);
namespace linear_cosmo { extern void addfckz(); }

// Fortran SAVE locals of fz2n / focd2z
static int fz2n_icalcn = 0;
static int fz2n_krinc = 0;
static int jindex_tab[257] = {0};  // 1-based 1..256

// ---------------------------------------------------------------------------
// focd2z — d-orbital two-center Fock (fock2z.F90:1254-1322)
//   fii/fjj/fij/pii/pjj/pij are 0-based segment pointers into f / ptot
//   (segment start = &f[i], &ptot[i]); wj/wk are 0-based integral segments.
// ---------------------------------------------------------------------------
static void focd2z_cpp(int iab, int jba, double* fii, double* fjj, double* fij,
                       const double* pii, const double* pjj, const double* pij,
                       const double* wj, const double* wk, bool diagonal, int& kr) {
    int loop = 0;
    for (int i = 1; i <= iab; ++i) {
        int ka = (i * (i - 1)) / 2;
        double aa = 2.0;
        for (int j = 1; j <= i; ++j) {
            if (i == j) aa = 1.0;
            int ij = ka + j;
            for (int k = 1; k <= jba; ++k) {
                int kc = (k * (k - 1)) / 2;
                double bb = 2.0;
                for (int l = 1; l <= k; ++l) {
                    if (k == l) bb = 1.0;
                    int kl = kc + l;
                    ++loop;
                    double a = w[loop - 1];                       // w(loop), 0-based
                    fii[ij - 1] += bb * a * pjj[kl - 1];
                    if (!diagonal) {
                        fjj[kl - 1] += aa * a * pii[ij - 1];
                        a = wk[loop - 1] * aa * bb * 0.125;
                        int ik = (i - 1) * jba + k;
                        int il = (i - 1) * jba + l;
                        int jk = (j - 1) * jba + k;
                        int jl = (j - 1) * jba + l;
                        fij[ik - 1] -= a * pij[jl - 1];
                        fij[il - 1] -= a * pij[jk - 1];
                        fij[jk - 1] -= a * pij[il - 1];
                        fij[jl - 1] -= a * pij[ik - 1];
                    }
                }
            }
        }
    }
    kr += loop;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// fz2n_pair - one atom pair (ii,jj) of fz2n (fock2z.F90:649-1252).
//   Writes into target Fock pointer fw (serial: global f; parallel:
//   per-thread private copy). kr/krinc are the sequential W-cursor and
//   focd2z loop counter (private per thread in the parallel path).
//   Per-thread buffers: tl_wjloc (rotation), tl_e1b/e2a, tl_pja/pjb, tl_fdummy.
// ---------------------------------------------------------------------------
static thread_local std::vector<double> tl_wjloc(2026, 0.0);
static thread_local double tl_e1b[46] = {0}, tl_e2a[46] = {0};
static thread_local double tl_pja[17] = {0}, tl_pjb[17] = {0};
static thread_local double tl_fdummy[172] = {0};

static void fz2n_pair(int ii, int jj, bool calci, bool calcj, double* fw,
                      int& kr, int& krinc,
                      const std::vector<double>& ptot,
                      const std::vector<std::vector<double>>& ptot2,
                      const std::vector<double>& qe) {
    int i, i1, ij, j, j1, ji, jk, jba, k, kj, kl, l, li, lii, lij, lj, ljj, m, ni, nj;
    double sum, sumdia, sumoff, ade, aee, da, dx, dy, dz, r, r2, ri2, ri5, rm,
        rp, w1, w2, w3, w4, w5, w6, w7, enuc, rij, point, constv;
    int iab = iorbs[ii];
                jba = iorbs[jj];
                if (nijbo[ii][jj] >= 0) {
                    if (direct && (calci || calcj)) {
                        double xi[3] = {coord[0][ii], coord[1][ii], coord[2][ii]};
                        double xj[3] = {coord[0][jj], coord[1][jj], coord[2][jj]};
                        rotate(nat[ii], nat[jj], xi, xj, tl_wjloc.data(), kr, tl_e1b, tl_e2a, enuc);
                    }
                    if (iab > 5 || jba > 5) {
                        if (!(calci || calcj)) {
                            if (!direct) kr = kr + (iab * (iab + 1)) / 2 * (jba * (jba + 1)) / 2;
                        } else {
                            i = nijbo[ii][ii] + 1;
                            j = nijbo[jj][jj] + 1;
                            ij = nijbo[ii][jj] + 1;
                            if (direct)
                                focd2z_cpp(iab, jba, &fw[i], &fw[j], &fw[ij], &ptot[i],
                                           &ptot[j], &ptot[ij], tl_wjloc.data(),
                                           tl_wjloc.data(), i == j, krinc);
                            else
                                focd2z_cpp(iab, jba, &fw[i], &fw[j], &fw[ij], &ptot[i],
                                           &ptot[j], &ptot[ij], &w[kr + 1],
                                           &wk[kr + 1], i == j, kr);
                        }
                    } else if (iab >= 3 && jba >= 3) {
                        //  HEAVY-ATOM - HEAVY-ATOM
                        if (calci || calcj) {
                            for (i = 1; i <= 16; ++i) {
                                tl_pja[i] = ptot2[ii][i];
                                tl_pjb[i] = ptot2[jj][i];
                            }
                            //  COULOMB TERMS
                            if (ii == jj) {
                                if (direct)
                                    jab_for_MOZYME(1, 1, tl_pja, tl_pjb, tl_wjloc.data(),
                                                   &fw[nijbo[ii][ii]], tl_fdummy);
                                else
                                    jab_for_MOZYME(1, 1, tl_pja, tl_pjb, &w[kr + 1],
                                                   &fw[nijbo[ii][ii]], tl_fdummy);
                            } else if (direct) {
                                jab_for_MOZYME(1, 1, tl_pja, tl_pjb, tl_wjloc.data(),
                                               &fw[nijbo[ii][ii]], &fw[nijbo[jj][jj]]);
                            } else {
                                jab_for_MOZYME(1, 1, tl_pja, tl_pjb, &w[kr + 1],
                                               &fw[nijbo[ii][ii]], &fw[nijbo[jj][jj]]);
                            }
                            //  EXCHANGE TERMS
                            l = nijbo[ii][jj] + 1;
                            if (ii != jj) {
                                if (direct)
                                    kab_for_MOZYME(0, 0, &ptot[nijbo[ii][jj]],
                                                   tl_wjloc.data(), &fw[nijbo[ii][jj]]);
                                else
                                    kab_for_MOZYME(0, 0, &ptot[nijbo[ii][jj]],
                                                   &wk[kr + 1], &fw[nijbo[ii][jj]]);
                            }
                        }
                        if (!direct) kr += 100;
                    } else if (iab >= 3 && jba == 1) {
                        //  LIGHT-ATOM - HEAVY-ATOM
                        if (calci || calcj) {
                            sumdia = 0.0;
                            sumoff = 0.0;
                            l = nijbo[ii][jj];
                            lj = nijbo[jj][jj] + 1;
                            li = nijbo[ii][ii];
                            k = 0;
                            if (direct) {
                                for (i = 1; i <= 4; ++i) {
                                    for (j = 1; j <= i - 1; ++j) {
                                        li = li + 1;
                                        k = k + 1;
                                        fw[li] += ptot[lj] * tl_wjloc[k - 1];
                                        sumoff += ptot[li] * tl_wjloc[k - 1];
                                    }
                                    li = li + 1;
                                    k = k + 1;
                                    fw[li] += ptot[lj] * tl_wjloc[k - 1];
                                    sumdia += ptot[li] * tl_wjloc[k - 1];
                                }
                                fw[lj] += sumoff * 2.0 + sumdia;
                                l = nijbo[ii][jj];
                                k = 0;
                                for (i = 1; i <= 4; ++i) {
                                    sum = 0.0;
                                    for (j = 1; j <= 4; ++j) {
                                        k = k + 1;
                                        sum += ptot[l + j] * tl_wjloc[jindex_tab[k] - 1];
                                    }
                                    fw[l + i] -= sum * 0.5;
                                }
                            } else {
                                for (i = 1; i <= 4; ++i) {
                                    for (j = 1; j <= i - 1; ++j) {
                                        li = li + 1;
                                        k = k + 1;
                                        fw[li] += ptot[lj] * w[kr + k];
                                        sumoff += ptot[li] * w[kr + k];
                                    }
                                    li = li + 1;
                                    k = k + 1;
                                    fw[li] += ptot[lj] * w[kr + k];
                                    sumdia += ptot[li] * w[kr + k];
                                }
                                fw[lj] += sumoff * 2.0 + sumdia;
                                l = nijbo[ii][jj];
                                k = 0;
                                for (i = 1; i <= 4; ++i) {
                                    sum = 0.0;
                                    for (j = 1; j <= 4; ++j) {
                                        k = k + 1;
                                        sum += ptot[l + j] * wk[kr + jindex_tab[k]];
                                    }
                                    fw[l + i] -= sum * 0.5;
                                }
                            }
                        }
                        if (!direct) kr += 10;
                    } else if (jba >= 3 && iab == 1) {
                        //  HEAVY-ATOM - LIGHT-ATOM
                        if (calci || calcj) {
                            sumdia = 0.0;
                            sumoff = 0.0;
                            l = nijbo[ii][jj];
                            lj = nijbo[jj][jj];
                            li = nijbo[ii][ii] + 1;
                            k = 0;
                            if (direct) {
                                for (i = 1; i <= 4; ++i) {
                                    for (j = 1; j <= i - 1; ++j) {
                                        k = k + 1;
                                        lj = lj + 1;
                                        fw[lj] += ptot[li] * tl_wjloc[k - 1];
                                        sumoff += ptot[lj] * tl_wjloc[k - 1];
                                    }
                                    lj = lj + 1;
                                    k = k + 1;
                                    fw[lj] += ptot[li] * tl_wjloc[k - 1];
                                    sumdia += ptot[lj] * tl_wjloc[k - 1];
                                }
                                fw[li] += sumoff * 2.0 + sumdia;
                                l = nijbo[ii][jj];
                                k = 0;
                                for (i = 1; i <= 4; ++i) {
                                    sum = 0.0;
                                    for (j = 1; j <= 4; ++j) {
                                        k = k + 1;
                                        sum += ptot[l + j] * tl_wjloc[jindex_tab[k] - 1];
                                    }
                                    fw[l + i] -= sum * 0.5;
                                }
                            } else {
                                for (i = 1; i <= 4; ++i) {
                                    for (j = 1; j <= i - 1; ++j) {
                                        k = k + 1;
                                        lj = lj + 1;
                                        fw[lj] += ptot[li] * w[kr + k];
                                        sumoff += ptot[lj] * w[kr + k];
                                    }
                                    lj = lj + 1;
                                    k = k + 1;
                                    fw[lj] += ptot[li] * w[kr + k];
                                    sumdia += ptot[lj] * w[kr + k];
                                }
                                fw[li] += sumoff * 2.0 + sumdia;
                                l = nijbo[ii][jj];
                                k = 0;
                                for (i = 1; i <= 4; ++i) {
                                    sum = 0.0;
                                    for (j = 1; j <= 4; ++j) {
                                        k = k + 1;
                                        sum += ptot[l + j] * wk[kr + jindex_tab[k]];
                                    }
                                    fw[l + i] -= sum * 0.5;
                                }
                            }
                        }
                        if (!direct) kr += 10;
                    } else if (iab == 1 && jba == 1) {
                        //  LIGHT-ATOM - LIGHT-ATOM
                        if (calci || calcj) {
                            lii = nijbo[ii][ii] + 1;
                            lij = nijbo[ii][jj] + 1;
                            ljj = nijbo[jj][jj] + 1;
                            if (direct) {
                                fw[lii] += ptot[ljj] * tl_wjloc[0];
                                if (ii != jj) {
                                    fw[ljj] += ptot[lii] * tl_wjloc[0];
                                    fw[lij] -= ptot[lij] * tl_wjloc[0] * 0.5;
                                }
                            } else {
                                fw[lii] += ptot[ljj] * w[kr + 1];
                                if (ii != jj) {
                                    fw[ljj] += ptot[lii] * w[kr + 1];
                                    fw[lij] -= ptot[lij] * wk[kr + 1] * 0.5;
                                }
                            }
                        }
                        if (!direct) kr += 1;
                    }
                } else {
                    //  Use point-charge approximation
                    i1 = nijbo[ii][ii];
                    j1 = nijbo[jj][jj];
                    if (iorbs[ii] * iorbs[jj] > 0) {
                        if (calci || calcj) {
                            if (semidr) {
                                dx = coord[0][ii] - coord[0][jj];
                                dy = coord[1][ii] - coord[1][jj];
                                dz = coord[2][ii] - coord[2][jj];
                                r2 = dx * dx + dy * dy + dz * dz;
                                ni = nat[ii];
                                nj = nat[jj];
                                aee = 0.5 / am[ni] + 0.5 / am[nj];
                                if (method_pm7) {
                                    rij = std::sqrt(r2);
                                    to_point(rij, point, constv);
                                    w1 = ev / std::sqrt(r2 / (a0 * a0) + aee * aee);
                                    w1 = w1 * constv + (1.0 - constv) * point;
                                } else {
                                    w1 = ev / std::sqrt(r2 / (a0 * a0) + aee * aee);
                                }
                            } else {
                                w1 = w[kr + 1];
                            }
                            //  MONOPOLE
                            for (i = 1; i <= iorbs[ii]; ++i) {
                                i1 = i1 + i;
                                fw[i1] += qe[jj] * w1;
                            }
                            for (j = 1; j <= iorbs[jj]; ++j) {
                                j1 = j1 + j;
                                fw[j1] += qe[ii] * w1;
                            }
                            i1 = nijbo[ii][ii];
                            j1 = nijbo[jj][jj];
                        }
                        if (!semidr) kr += 1;
                        //  DIPOLE, IF NEEDED
                        if (nijbo[ii][jj] == -2) {
                            if ((calci || calcj) && semidr) {
                                r = std::sqrt(r2);
                                dx = dx / r;
                                dy = dy / r;
                                dz = dz / r;
                                if (std::fabs(dz) > 0.99999999) dz = (dz > 0 ? 1.0 : -1.0);
                                if (iorbs[ii] > 1) {
                                    da = dd[ni];
                                    ade = 0.5 / ad[ni] + 0.5 / am[nj];
                                    rp = std::sqrt((r / a0 + da) * (r / a0 + da) + ade * ade);
                                    rm = std::sqrt((r / a0 - da) * (r / a0 - da) + ade * ade);
                                    ri2 = ev * (0.5 / rp - 0.5 / rm);
                                    if (method_pm7) {
                                        to_point(r, point, constv);
                                        ri2 = ri2 * constv;
                                    }
                                    w5 = ri2 * dx;
                                    w6 = ri2 * dy;
                                    w7 = ri2 * dz;
                                }
                                if (iorbs[jj] > 1) {
                                    da = dd[nj];
                                    ade = 0.5 / am[ni] + 0.5 / ad[nj];
                                    rp = std::sqrt((r / a0 + da) * (r / a0 + da) + ade * ade);
                                    rm = std::sqrt((r / a0 - da) * (r / a0 - da) + ade * ade);
                                    ri5 = -ev * (0.5 / rp - 0.5 / rm);
                                    if (method_pm7) {
                                        to_point(r, point, constv);
                                        ri5 = ri5 * constv;
                                    }
                                    w2 = ri5 * dx;
                                    w3 = ri5 * dy;
                                    w4 = ri5 * dz;
                                }
                            }
                            if (iorbs[ii] > 1) {
                                if (iorbs[jj] > 1) {
                                    if (calci || calcj) {
                                        if (!semidr) {
                                            w2 = w[kr + 1];
                                            w3 = w[kr + 2];
                                            w4 = w[kr + 3];
                                            w5 = w[kr + 4];
                                            w6 = w[kr + 5];
                                            w7 = w[kr + 6];
                                        }
                                        fw[j1 + 2] += qe[ii] * w2;
                                        fw[j1 + 4] += qe[ii] * w3;
                                        fw[j1 + 7] += qe[ii] * w4;
                                        fw[i1 + 2] += qe[jj] * w5;
                                        fw[i1 + 4] += qe[jj] * w6;
                                        fw[i1 + 7] += qe[jj] * w7;
                                        fw[j1 + 1] += ptot[i1 + 2] * w5 * 2;
                                        fw[j1 + 3] += ptot[i1 + 2] * w5 * 2;
                                        fw[j1 + 6] += ptot[i1 + 2] * w5 * 2;
                                        fw[j1 + 10] += ptot[i1 + 2] * w5 * 2;
                                        fw[j1 + 1] += ptot[i1 + 4] * w6 * 2;
                                        fw[j1 + 3] += ptot[i1 + 4] * w6 * 2;
                                        fw[j1 + 6] += ptot[i1 + 4] * w6 * 2;
                                        fw[j1 + 10] += ptot[i1 + 4] * w6 * 2;
                                        fw[j1 + 1] += ptot[i1 + 7] * w7 * 2;
                                        fw[j1 + 3] += ptot[i1 + 7] * w7 * 2;
                                        fw[j1 + 6] += ptot[i1 + 7] * w7 * 2;
                                        fw[j1 + 10] += ptot[i1 + 7] * w7 * 2;
                                        fw[i1 + 1] += ptot[j1 + 2] * w2 * 2;
                                        fw[i1 + 3] += ptot[j1 + 2] * w2 * 2;
                                        fw[i1 + 6] += ptot[j1 + 2] * w2 * 2;
                                        fw[i1 + 10] += ptot[j1 + 2] * w2 * 2;
                                        fw[i1 + 1] += ptot[j1 + 4] * w3 * 2;
                                        fw[i1 + 3] += ptot[j1 + 4] * w3 * 2;
                                        fw[i1 + 6] += ptot[j1 + 4] * w3 * 2;
                                        fw[i1 + 10] += ptot[j1 + 4] * w3 * 2;
                                        fw[i1 + 1] += ptot[j1 + 7] * w4 * 2;
                                        fw[i1 + 3] += ptot[j1 + 7] * w4 * 2;
                                        fw[i1 + 6] += ptot[j1 + 7] * w4 * 2;
                                        fw[i1 + 10] += ptot[j1 + 7] * w4 * 2;
                                    }
                                    if (!semidr) kr += 6;
                                } else {
                                    if (calci || calcj) {
                                        if (!semidr) {
                                            w5 = w[kr + 1];
                                            w6 = w[kr + 2];
                                            w7 = w[kr + 3];
                                        }
                                        fw[i1 + 2] += qe[jj] * w5;
                                        fw[i1 + 4] += qe[jj] * w6;
                                        fw[i1 + 7] += qe[jj] * w7;
                                        fw[j1 + 1] += ptot[i1 + 2] * w5 * 2;
                                        fw[j1 + 1] += ptot[i1 + 4] * w6 * 2;
                                        fw[j1 + 1] += ptot[i1 + 7] * w7 * 2;
                                    }
                                    if (!semidr) kr += 3;
                                }
                            } else if (iorbs[jj] > 1) {
                                if (calci || calcj) {
                                    if (!semidr) {
                                        w2 = w[kr + 1];
                                        w3 = w[kr + 2];
                                        w4 = w[kr + 3];
                                    }
                                    fw[j1 + 2] += qe[ii] * w2;
                                    fw[j1 + 4] += qe[ii] * w3;
                                    fw[j1 + 7] += qe[ii] * w4;
                                    fw[i1 + 1] += ptot[j1 + 2] * w2 * 2;
                                    fw[i1 + 1] += ptot[j1 + 4] * w3 * 2;
                                    fw[i1 + 1] += ptot[j1 + 7] * w4 * 2;
                                }
                                if (!semidr) kr += 3;
                            }
                        }
                    }
    }
}

// fz2n — fast two-center Fock (fock2z.F90:641-1253), lijbo=true path.
//   wjloc here is a 0-based 2026-element work buffer filled by rotate()
//   (rotate writes w[0..len-1]); jab/kab/focd2z read 0-based w.
// ---------------------------------------------------------------------------
static void fz2n_cpp(std::vector<double>& f, const std::vector<double>& ptot,
                     std::vector<double>& q, std::vector<double>& qe,
                     std::vector<std::vector<double>>& ptot2, int mode, int ione) {
    bool calci, calcj;
    int i, i1, iab, ii, iim1, ij, ilim, ired, j, j1, jba, ji, jj, jk, jred, k,
        kj, kl, kr, l, li, lii, lij, lj, ljj, lk, m, ni, nj;
    double sum, sumdia, sumoff, ade, aee, da, dx, dy, dz, r, r2, ri2, ri5, rm,
        rp, w1, w2, w3, w4, w5, w6, w7, enuc, rij, point, constv;
    double pja[17], pjb[17];     // 1-based 1..16
    double e1b[46], e2a[46];     // 1-based 1..45
    double fdummy[172] = {};     // 1-based 1..171
    std::vector<double> wjloc(2026, 0.0);  // 0-based
    std::vector<char> calci_tab(numat + 1, 0);
    std::vector<std::vector<char>> calcj_tab(numat + 1);

    // Recreate the original per-call zero-initialized work buffers: the
    // original fz2n_cpp allocated wjloc(2026, 0.0) on every call, while the
    // thread-local versions persist across calls, so stale data could be
    // read where a rotation wrote fewer elements than focd2z/jab read.
    std::fill(tl_wjloc.begin(), tl_wjloc.end(), 0.0);
    std::fill(tl_e1b, tl_e1b + 46, 0.0);
    std::fill(tl_e2a, tl_e2a + 46, 0.0);
    std::fill(tl_pja, tl_pja + 17, 0.0);
    std::fill(tl_pjb, tl_pjb + 17, 0.0);

    if (fz2n_icalcn != numcal) {
        fz2n_icalcn = numcal;
        m = 0;
        for (i = 1; i <= 4; ++i)
            for (j = 1; j <= 4; ++j) {
                ij = std::min(i, j);
                ji = i + j - ij;
                for (k = 1; k <= 4; ++k)
                    for (l = 1; l <= 4; ++l) {
                        m = m + 1;
                        kl = std::min(k, l);
                        lk = k + l - kl;
                        jindex_tab[m] = (ifact[ji] + ij) * 10 + ifact[lk] + kl - 10;
                    }
            }
    }
    //  Put P(A,A) density into array PTOT2
    chrge_for_MOZYME(ptot.data(), qe.data(), numat, iorbs.data());
    if (mode == -1) {
        for (int kk = 1; kk <= mpack; ++kk) f[kk] = -f[kk];
    }
    l = 0;
    for (ii = 1; ii <= numat; ++ii) {
        q[ii] = tore[nat[ii]] - qe[ii];
        i = nijbo[ii][ii];
        iab = iorbs[ii];
        m = 0;
        for (j = 1; j <= iab; ++j)
            for (k = 1; k <= iab; ++k) {
                m = m + 1;
                jk = std::min(j, k);
                kj = k + j - jk;
                ptot2[ii][m] = ptot[i + (kj * (kj - 1)) / 2 + jk];
            }
    }
    kr = 0;
    {
        // Precompute calci/calcj lookup tables: bitwise-identical values to
        // the original sequential ired/jred cursors, but without loop-carried
        // state, so the atom-pair loop can be parallelized.
        int ired = 1;
        for (int ii = 1; ii <= numat; ++ii) {
            bool ci = (kopt[ired] == ii);
            if (ci && ired < numat) ++ired;
            if (mode == 0) ci = true;
            calci_tab[ii] = ci ? 1 : 0;
            calcj_tab[ii].assign(ii, 0);      // jj in 1..ii-1
            int jred = 1;
            for (int jj = 1; jj <= ii - 1; ++jj) {
                bool cj = (kopt[jred] == jj);
                if (cj && jred < numat) ++jred;
                calcj_tab[ii][jj] = cj ? 1 : 0;
            }
        }
    }
    const bool fz2n_par = direct && std::getenv("MOPAC_FZ2N_OMP") != nullptr;
    if (fz2n_par) {
        // Parallel path: each thread accumulates into a private Fock copy,
        // then a deterministic serial reduction merges them (paper C2).
        int nt = omp_get_max_threads();
        std::vector<std::vector<double>> fpriv(nt, std::vector<double>(mpack + 1, 0.0));
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            double* fw = fpriv[tid].data();
            int kr_loc = 0, krinc_loc = 0;
            std::fill(tl_wjloc.begin(), tl_wjloc.end(), 0.0);
            std::fill(tl_e1b, tl_e1b + 46, 0.0);
            std::fill(tl_e2a, tl_e2a + 46, 0.0);
            std::fill(tl_pja, tl_pja + 17, 0.0);
            std::fill(tl_pjb, tl_pjb + 17, 0.0);
#pragma omp for schedule(dynamic)
            for (int ii = 1; ii <= numat; ++ii) {
                if (iorbs[ii] == 0) continue;
                const bool calci = calci_tab[ii] != 0;
                for (int jj = 1; jj <= ii - 1; ++jj)
                    fz2n_pair(ii, jj, calci, calcj_tab[ii][jj] != 0, fw,
                              kr_loc, krinc_loc, ptot, ptot2, qe);
            }
        }
        for (int t = 0; t < nt; ++t)
            for (int q = 1; q <= mpack; ++q) f[q] += fpriv[t][q];
    } else {
        // Serial path: original loop-carried order (kr advances atom-pair by
        // atom-pair; diagonal fock1 segments interleaved per atom).
        for (int ii = 1; ii <= numat; ++ii) {
            const bool calci = calci_tab[ii] != 0;
            int iab = iorbs[ii];
            if (iab != 0) {
                for (int jj = 1; jj <= ii - 1; ++jj)
                    fz2n_pair(ii, jj, calci, calcj_tab[ii][jj] != 0, f.data(),
                              kr, fz2n_krinc, ptot, ptot2, qe);
                if (!direct) {
                    int i = nijbo[ii][ii] + 1;
                    int ilim = (iab * (iab + 1)) / 2;
                    std::vector<double> fseg(ilim + 1), pseg(ilim + 1);
                    for (int t = 1; t <= ilim; ++t) {
                        fseg[t] = f[i + t - 1];
                        pseg[t] = ptot[i + t - 1];
                    }
                    std::vector<std::vector<double>> w2(ilim + 1,
                                                        std::vector<double>(ilim + 1, 0.0));
                    for (int row = 1; row <= ilim; ++row)
                        for (int col = 1; col <= ilim; ++col)
                            w2[row][col] = w[kr + (col - 1) * ilim + (row - 1) + 1];
                    fock1_for_MOZYME(fseg, pseg, w2, kr, iab, ilim);
                    for (int t = 1; t <= ilim; ++t) f[i + t - 1] = fseg[t];
                }
            }
        }
    }

    if (direct) {
        kr = 0;
        ired = 1;
        for (ii = 1; ii <= numat; ++ii) {
            iab = iorbs[ii];
            if (iab != 0) {
                i = nijbo[ii][ii] + 1;
                ilim = (iab * (iab + 1)) / 2;
                std::vector<double> fseg(ilim + 1), pseg(ilim + 1);
                for (int t = 1; t <= ilim; ++t) {
                    fseg[t] = f[i + t - 1];
                    pseg[t] = ptot[i + t - 1];
                }
                std::vector<std::vector<double>> w2(ilim + 1,
                                                    std::vector<double>(ilim + 1, 0.0));
                for (int row = 1; row <= ilim; ++row)
                    for (int col = 1; col <= ilim; ++col)
                        w2[row][col] = w[kr + (col - 1) * ilim + (row - 1) + 1];
                fock1_for_MOZYME(fseg, pseg, w2, kr, iab, ilim);
                for (int t = 1; t <= ilim; ++t) f[i + t - 1] = fseg[t];
            }
        }
    }
    if (mode == -1) {
        for (int kk = 1; kk <= mpack; ++kk) f[kk] = -f[kk];
    }
    //  The following routine adds the dielectric correction to F.
    if (useps && mode == 0) linear_cosmo::addfckz();
}

void buildf(std::vector<double>& F, const std::vector<double>& partf, int mode) {
    std::vector<double> q(numat + 1, 0.0), qe(numat + 1, 0.0);
    std::vector<std::vector<double>> ptot2(numat + 1, std::vector<double>(82, 0.0));

    if (mode == -1) {
        for (int k = 1; k <= mpack; ++k) F[k] = partf[k] - h[k];
    } else if (mode == 0) {
        for (int k = 1; k <= mpack; ++k) F[k] = h[k];
    } else {  // mode == 1
        for (int k = 1; k <= mpack; ++k) F[k] = partf[k] + h[k];
    }
    if (id == 0) {
        //  Discrete system: fock2z -> fz2n when lijbo (fast two-center).
        //  (fz2 full path is not ported; lijbo=true is the default MOZYME path.)
        fz2n_cpp(F, p, q, qe, ptot2, mode, 1);
    } else {
        //  A solid (polymer, layer, or 3-D solid)
        fz2n_cpp(F, p, q, qe, ptot2, mode, 0);
    }
}
