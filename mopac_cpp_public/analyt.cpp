// analyt.cpp — C++ translation of MOPAC 2016 "analyt.F90".
#pragma warning(disable: 4459)

#include "analyt.h"

#include <algorithm>
#include <cmath>

#include "analyt_C.h"
#include "chanel_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "overlaps_C.h"
#include "parameters_C.h"

using namespace analyt_C;
using namespace parameters_C;

namespace {
int icalcn = 0;
bool am1 = false;
}  // namespace

// ===========================================================================
void analyt(const std::vector<double>& psum, const std::vector<double>& palpha,
            const std::vector<double>& pbeta,
            const std::vector<std::vector<double>>& coord,
            const int nat[3], int jja, int jjd, int iia, int iid,
            double eng[4]) {
    if (icalcn != molkst_C::numcal) {
        icalcn = molkst_C::numcal;
        am1 = molkst_C::keywrd.find("AM1") != std::string::npos ||
              molkst_C::keywrd.find("PM3") != std::string::npos;
    }
    int jd = jjd - jja + 1;
    int ja = 1;
    int id = iid - iia + 1 + jd;
    int ia = jd + 1;

    double eaa[4] = {}, eab[4] = {}, enuc[4] = {};
    eng[1] = eng[2] = eng[3] = 0.0;

    int i = 2, j = 1;
    int ni = nat[i], nj = nat[j];
    int istart = nztype[ni] * 4 - 3;
    int jstart = nztype[nj] * 4 - 3;

    double r2 = (coord[0][i] - coord[0][j]) * (coord[0][i] - coord[0][j]) +
                (coord[1][i] - coord[1][j]) * (coord[1][i] - coord[1][j]) +
                (coord[2][i] - coord[2][j]) * (coord[2][i] - coord[2][j]);
    double rij = std::sqrt(r2);
    double r0 = rij / funcon_C::a0;
    double rr = r2 / (funcon_C::a0 * funcon_C::a0);

    for (int ix = 1; ix <= 3; ++ix) {
        double del1 = coord[ix-1][i] - coord[ix-1][j];
        double termaa = 0.0, termab = 0.0;
        int isp = 0, iol = 0;

        // First derivatives of overlap integrals.
        for (int k = ia; k <= id; ++k) {
            int ka = k - ia;
            int kg = istart + ka;
            for (int l = ja; l <= jd; ++l) {
                int la = l - ja;
                int lg = jstart + la;
                iol = iol + 1;
                ds[iol] = 0.0;
                int is;
                double del2 = 0.0, del3 = 0.0;
                if (ka == 0 && la == 0) {
                    if (std::fabs(del1) <= 1.0e-6) continue;
                    is = 1;
                } else if (ka == 0 && la > 0) {
                    is = 3;
                    if (ix == la) {
                        // fall through to call ders
                    } else if (std::fabs(del1) <= 1.0e-6) {
                        continue;
                    } else {
                        is = 2;
                        del2 = coord[la-1][i] - coord[la-1][j];
                    }
                } else if (ka > 0 && la == 0) {
                    is = 5;
                    if (ix == ka) {
                        // fall through
                    } else if (std::fabs(del1) <= 1.0e-6) {
                        continue;
                    } else {
                        is = 4;
                        del2 = coord[ka-1][i] - coord[ka-1][j];
                    }
                } else {
                    if (ka == la) {
                        is = 9;
                        if (ix == ka) {
                            // fall through
                        } else if (std::fabs(del1) <= 1.0e-6) {
                            continue;
                        } else {
                            is = 8;
                            del2 = coord[ka-1][i] - coord[ka-1][j];
                        }
                    } else if (ix != ka && ix != la) {
                        if (std::fabs(del1) <= 1.0e-6) continue;
                        is = 7;
                        del2 = coord[ka-1][i] - coord[ka-1][j];
                        del3 = coord[la-1][i] - coord[la-1][j];
                    } else {
                        del2 = coord[ka + la - ix - 1][i] - coord[ka + la - ix - 1][j];
                        is = 6;
                    }
                }
                ders(kg, lg, rr, del1, del2, del3, is, iol);
            }
        }
        if (ix == 1) {
            for (int i22 = 1; i22 <= 22; ++i22) g[i22] = g[i22];  // read(irot) g
        }
        delri(dg, ni, nj, r0, del1);
        delmol(coord, i, j, ni, nj, ia, id, ja, jd, ix, rij, del1, isp);

        // First derivative of nuclear repulsion.
        double termnc;
        if (rij < 1.0 && natorb[ni] * natorb[nj] == 0) {
            termnc = 0.0;
        } else {
            double c1 = tore[ni] * tore[nj];
            double dd;
            if (ni == 1 && (nj == 7 || nj == 8)) {
                double f3 = 1.0 + std::exp(-alp[1] * rij) + rij * std::exp(-alp[nj] * rij);
                dd = (dg[1] * f3 -
                      g[1] * (del1 / rij) *
                          (alp[1] * std::exp(-alp[1] * rij) +
                           (alp[nj] * rij - 1.0) * std::exp(-alp[nj] * rij))) *
                     c1;
            } else if ((ni == 7 || ni == 8) && nj == 1) {
                double f3 = 1.0 + std::exp(-alp[1] * rij) + rij * std::exp(-alp[ni] * rij);
                dd = (dg[1] * f3 -
                      g[1] * (del1 / rij) *
                          (alp[1] * std::exp(-alp[1] * rij) +
                           (alp[ni] * rij - 1.0) * std::exp(-alp[ni] * rij))) *
                     c1;
            } else {
                double part1 = dg[1] * c1;
                double part2 = -(g[1] * (del1 / rij) *
                                 (alp[ni] * std::exp(-alp[ni] * rij) +
                                  alp[nj] * std::exp(-alp[nj] * rij))) *
                               std::fabs(c1);
                double part3 = dg[1] *
                               (std::exp(-alp[ni] * rij) + std::exp(-alp[nj] * rij)) *
                               std::fabs(c1);
                dd = part1 + part2 + part3;
            }
            termnc = dd;
            if (am1) {
                double anam1 = 0.0;
                for (int ig = 1; ig <= 4; ++ig) {
                    if (std::fabs(guess1[ni][ig]) > 0.0)
                        anam1 += guess1[ni][ig] *
                                 (1.0 / (rij * rij) +
                                  2.0 * guess2[ni][ig] * (rij - guess3[ni][ig]) / rij) *
                                 std::exp(std::max(-30.0, -guess2[ni][ig] *
                                                              (rij - guess3[ni][ig]) *
                                                              (rij - guess3[ni][ig])));
                    if (std::fabs(guess1[nj][ig]) <= 0.0) continue;
                    anam1 += guess1[nj][ig] *
                             (1.0 / (rij * rij) +
                              2.0 * guess2[nj][ig] * (rij - guess3[nj][ig]) / rij) *
                             std::exp(std::max(-30.0, -guess2[nj][ig] *
                                                          (rij - guess3[nj][ig]) *
                                                          (rij - guess3[nj][ig])));
                }
                anam1 = anam1 * tore[ni] * tore[nj];
                termnc = termnc - anam1 * del1 / rij;
            }
        }

        double bi[5] = {}, bj[5] = {};
        bi[1] = betas[ni]; bi[2] = betap[ni]; bi[3] = bi[2]; bi[4] = bi[2];
        bj[1] = betas[nj]; bj[2] = betap[nj]; bj[3] = bj[2]; bj[4] = bj[2];

        iol = 0;
        for (int k = ia; k <= id; ++k) {
            if (jd - ja + 1 > 0) {
                double termk = bi[k - ia + 1];
                int base = ja + k * (k - 1) / 2;
                for (int l = 1; l <= jd - ja + 1; ++l)
                    termab += (termk + bj[l]) * psum[base + (l - 1)] * ds[iol + l];
                iol = jd - ja + 1 + iol;
            }
        }

        // Core-electron attraction: atom core i affecting AOs on j.
        isp = 0;
        for (int m = ja; m <= jd; ++m) {
            double bb = 1.0;
            for (int n = m; n <= jd; ++n) {
                int mn = m + (n * (n - 1)) / 2;
                isp = isp + 1;
                termab = termab - bb * tore[ni] * psum[mn] * dr[isp];
                bb = 2.0;
            }
        }
        // Atom core j affecting AOs on i.
        int k = std::max(jd - ja + 1, 1);
        k = (k * (k + 1)) / 2;
        isp = -k + 1;
        for (int m = ia; m <= id; ++m) {
            double bb = 1.0;
            for (int n = m; n <= id; ++n) {
                int mn = m + (n * (n - 1)) / 2;
                isp = isp + k;
                termab = termab - bb * tore[nj] * psum[mn] * dr[isp];
                bb = 2.0;
            }
        }
        isp = 0;

        // Coulomb and exchange.
        for (int kk2 = ia; kk2 <= id; ++kk2) {
            double aa = 1.0;
            int kk = (kk2 * (kk2 - 1)) / 2;
            for (int ll = kk2; ll <= id; ++ll) {
                int llb = (ll * (ll - 1)) / 2;
                for (int m = ja; m <= jd; ++m) {
                    double bb = 1.0;
                    for (int n = m; n <= jd; ++n) {
                        isp = isp + 1;
                        int kl = kk2 + llb;
                        int mn = m + (n * (n - 1)) / 2;
                        termaa = termaa + aa * bb * psum[kl] * psum[mn] * dr[isp];
                        int mk = m + kk, nk = n + kk, ml = m + llb, nl = n + llb;
                        termaa = termaa -
                                 0.5 * aa * bb *
                                     (palpha[mk] * palpha[nl] + palpha[nk] * palpha[ml] +
                                      pbeta[mk] * pbeta[nl] + pbeta[nk] * pbeta[ml]) *
                                     dr[isp];
                        bb = 2.0;
                    }
                }
            }
            aa = 2.0;
        }
        eaa[ix] += termaa;
        eab[ix] += termab;
        enuc[ix] += termnc;
    }
    eng[1] = eaa[1] + eab[1] + enuc[1];
    eng[2] = eaa[2] + eab[2] + enuc[2];
    eng[3] = eaa[3] + eab[3] + enuc[3];
    eng[1] = -eng[1] * funcon_C::fpc_9;
    eng[2] = -eng[2] * funcon_C::fpc_9;
    eng[3] = -eng[3] * funcon_C::fpc_9;
}

// ===========================================================================
void delmol(const std::vector<std::vector<double>>& coord, int i, int j,
            int ni, int nj, int ia, int id, int ja, int jd, int ix,
            double rij, double tomb, int& isp) {
    if (ni > 1 || nj > 1) rotat(coord, i, j, ix, 2, rij, tomb);
    int ib = std::max(ia, id);
    int jb = std::max(ja, jd);
    for (int k = ia; k <= ib; ++k) {
        int kk = k - ia;
        for (int l = k; l <= ib; ++l) {
            int ll = l - ia;
            for (int m = ja; m <= jb; ++m) {
                int mm = m - ja;
                for (int n = m; n <= jb; ++n) {
                    int nn = n - ja;
                    isp = isp + 1;
                    if (nn == 0) {
                        if (ll == 0) {
                            dr[isp] = dg[1];
                        } else if (kk == 0) {
                            dr[isp] = dg[2] * tx[ll] + g[2] * tdx[ll];
                        } else {
                            dr[isp] = dg[3] * tx[kk] * tx[ll] +
                                      g[3] * (tdx[kk] * tx[ll] + tx[kk] * tdx[ll]) +
                                      dg[4] * (ty[kk] * ty[ll] + tz[kk] * tz[ll]) +
                                      g[4] * (tdy[kk] * ty[ll] + ty[kk] * tdy[ll] +
                                              tdz[kk] * tz[ll] + tz[kk] * tdz[ll]);
                        }
                    } else if (mm == 0) {
                        if (ll == 0) {
                            dr[isp] = dg[5] * tx[nn] + g[5] * tdx[nn];
                        } else if (kk == 0) {
                            dr[isp] = dg[6] * tx[ll] * tx[nn] +
                                      g[6] * (tdx[ll] * tx[nn] + tx[ll] * tdx[nn]) +
                                      dg[7] * (ty[ll] * ty[nn] + tz[ll] * tz[nn]) +
                                      g[7] * (tdy[ll] * ty[nn] + ty[ll] * tdy[nn] +
                                              tdz[ll] * tz[nn] + tz[ll] * tdz[nn]);
                        } else {
                            dr[isp] = dg[8] * tx[kk] * tx[ll] * tx[nn] +
                                      g[8] * (tdx[kk] * tx[ll] * tx[nn] +
                                              tx[kk] * tdx[ll] * tx[nn] +
                                              tx[kk] * tx[ll] * tdx[nn]) +
                                      dg[9] * (ty[kk] * ty[ll] + tz[kk] * tz[ll]) * tx[nn] +
                                      g[9] * ((tdy[kk] * ty[ll] + ty[kk] * tdy[ll] +
                                               tdz[kk] * tz[ll] + tz[kk] * tdz[ll]) *
                                                  tx[nn] +
                                              (ty[kk] * ty[ll] + tz[kk] * tz[ll]) * tdx[nn]) +
                                      dg[10] * (tx[kk] * (ty[ll] * ty[nn] + tz[ll] * tz[nn]) +
                                                tx[ll] * (ty[kk] * ty[nn] + tz[kk] * tz[nn])) +
                                      g[10] * (tdx[kk] * (ty[ll] * ty[nn] + tz[ll] * tz[nn]) +
                                               tdx[ll] * (ty[kk] * ty[nn] + tz[kk] * tz[nn]) +
                                               tx[kk] * (tdy[ll] * ty[nn] + ty[ll] * tdy[nn] +
                                                         tdz[ll] * tz[nn] + tz[ll] * tdz[nn]) +
                                               tx[ll] * (tdy[kk] * ty[nn] + ty[kk] * tdy[nn] +
                                                         tdz[kk] * tz[nn] + tz[kk] * tdz[nn]));
                        }
                    } else if (ll == 0) {
                        dr[isp] = dg[11] * tx[mm] * tx[nn] +
                                  g[11] * (tdx[mm] * tx[nn] + tx[mm] * tdx[nn]) +
                                  dg[12] * (ty[mm] * ty[nn] + tz[mm] * tz[nn]) +
                                  g[12] * (tdy[mm] * ty[nn] + ty[mm] * tdy[nn] +
                                           tdz[mm] * tz[nn] + tz[mm] * tdz[nn]);
                    } else if (kk == 0) {
                        dr[isp] = dg[13] * tx[ll] * tx[mm] * tx[nn] +
                                  g[13] * (tdx[ll] * tx[mm] * tx[nn] +
                                           tx[ll] * tdx[mm] * tx[nn] +
                                           tx[ll] * tx[mm] * tdx[nn]) +
                                  dg[14] * tx[ll] * (ty[mm] * ty[nn] + tz[mm] * tz[nn]) +
                                  g[14] * (tdx[ll] * (ty[mm] * ty[nn] + tz[mm] * tz[nn]) +
                                           tx[ll] * (tdy[mm] * ty[nn] + ty[mm] * tdy[nn] +
                                                     tdz[mm] * tz[nn] + tz[mm] * tdz[nn])) +
                                  dg[15] * (ty[ll] * (ty[mm] * tx[nn] + ty[nn] * tx[mm]) +
                                            tz[ll] * (tz[mm] * tx[nn] + tz[nn] * tx[mm])) +
                                  g[15] * (tdy[ll] * (ty[mm] * tx[nn] + ty[nn] * tx[mm]) +
                                           tdz[ll] * (tz[mm] * tx[nn] + tz[nn] * tx[mm]) +
                                           ty[ll] * (tdy[mm] * tx[nn] + ty[mm] * tdx[nn] +
                                                     tdy[nn] * tx[mm] + ty[nn] * tdx[mm]) +
                                           tz[ll] * (tdz[mm] * tx[nn] + tz[mm] * tdx[nn] +
                                                     tdz[nn] * tx[mm] + tz[nn] * tdx[mm]));
                    } else {
                        dr[isp] = dg[16] * tx[kk] * tx[ll] * tx[mm] * tx[nn] +
                                  g[16] * (tdx[kk] * tx[ll] * tx[mm] * tx[nn] +
                                           tx[kk] * tdx[ll] * tx[mm] * tx[nn] +
                                           tx[kk] * tx[ll] * tdx[mm] * tx[nn] +
                                           tx[kk] * tx[ll] * tx[mm] * tdx[nn]) +
                                  dg[17] * (ty[kk] * ty[ll] + tz[kk] * tz[ll]) * tx[mm] * tx[nn] +
                                  g[17] * ((tdy[kk] * ty[ll] + ty[kk] * tdy[ll] +
                                            tdz[kk] * tz[ll] + tz[kk] * tdz[ll]) *
                                               tx[mm] * tx[nn] +
                                           (ty[kk] * ty[ll] + tz[kk] * tz[ll]) *
                                               (tdx[mm] * tx[nn] + tx[mm] * tdx[nn])) +
                                  dg[18] * tx[kk] * tx[ll] *
                                      (ty[mm] * ty[nn] + tz[mm] * tz[nn]) +
                                  g[18] * ((tdx[kk] * tx[ll] + tx[kk] * tdx[ll]) *
                                               (ty[mm] * ty[nn] + tz[mm] * tz[nn]) +
                                           tx[kk] * tx[ll] * (tdy[mm] * ty[nn] +
                                                             ty[mm] * tdy[nn] +
                                                             tdz[mm] * tz[nn] +
                                                             tz[mm] * tdz[nn]));
                        dr[isp] = dr[isp] + dg[19] * (ty[kk] * ty[ll] * ty[mm] * ty[nn] +
                                                      tz[kk] * tz[ll] * tz[mm] * tz[nn]) +
                                  g[19] * (tdy[kk] * ty[ll] * ty[mm] * ty[nn] +
                                           ty[kk] * tdy[ll] * ty[mm] * ty[nn] +
                                           ty[kk] * ty[ll] * tdy[mm] * ty[nn] +
                                           ty[kk] * ty[ll] * ty[mm] * tdy[nn] +
                                           tdz[kk] * tz[ll] * tz[mm] * tz[nn] +
                                           tz[kk] * tdz[ll] * tz[mm] * tz[nn] +
                                           tz[kk] * tz[ll] * tdz[mm] * tz[nn] +
                                           tz[kk] * tz[ll] * tz[mm] * tdz[nn]) +
                                  dg[20] * (tx[kk] * (tx[mm] * (ty[ll] * ty[nn] +
                                                                tz[ll] * tz[nn]) +
                                                      tx[nn] * (ty[ll] * ty[mm] +
                                                                tz[ll] * tz[mm])) +
                                            tx[ll] * (tx[mm] * (ty[kk] * ty[nn] +
                                                                tz[kk] * tz[nn]) +
                                                      tx[nn] * (ty[kk] * ty[mm] +
                                                                tz[kk] * tz[mm])));
                        double temp1 =
                            tdx[kk] * (tx[mm] * (ty[ll] * ty[nn] + tz[ll] * tz[nn]) +
                                       tx[nn] * (ty[ll] * ty[mm] + tz[ll] * tz[mm])) +
                            tdx[ll] * (tx[mm] * (ty[kk] * ty[nn] + tz[kk] * tz[nn]) +
                                       tx[nn] * (ty[kk] * ty[mm] + tz[kk] * tz[mm])) +
                            tx[kk] * (tdx[mm] * (ty[ll] * ty[nn] + tz[ll] * tz[nn]) +
                                      tdx[nn] * (ty[ll] * ty[mm] + tz[ll] * tz[mm])) +
                            tx[ll] * (tdx[mm] * (ty[kk] * ty[nn] + tz[kk] * tz[nn]) +
                                      tdx[nn] * (ty[kk] * ty[mm] + tz[kk] * tz[mm]));
                        double temp2 =
                            tx[kk] * (tx[mm] * (tdy[ll] * ty[nn] + ty[ll] * tdy[nn] +
                                               tdz[ll] * tz[nn] + tz[ll] * tdz[nn]) +
                                      tx[nn] * (tdy[ll] * ty[mm] + ty[ll] * tdy[mm] +
                                               tdz[ll] * tz[mm] + tz[ll] * tdz[mm])) +
                            tx[ll] * (tx[mm] * (tdy[kk] * ty[nn] + ty[kk] * tdy[nn] +
                                               tdz[kk] * tz[nn] + tz[kk] * tdz[nn]) +
                                      tx[nn] * (tdy[kk] * ty[mm] + ty[kk] * tdy[mm] +
                                               tdz[kk] * tz[mm] + tz[kk] * tdz[mm]));
                        dr[isp] = dr[isp] + g[20] * (temp1 + temp2);
                        dr[isp] = dr[isp] +
                                  dg[21] * (ty[kk] * ty[ll] * tz[mm] * tz[nn] +
                                            tz[kk] * tz[ll] * ty[mm] * ty[nn]) +
                                  g[21] * (tdy[kk] * ty[ll] * tz[mm] * tz[nn] +
                                           ty[kk] * tdy[ll] * tz[mm] * tz[nn] +
                                           ty[kk] * ty[ll] * tdz[mm] * tz[nn] +
                                           ty[kk] * ty[ll] * tz[mm] * tdz[nn] +
                                           tdz[kk] * tz[ll] * ty[mm] * ty[nn] +
                                           tz[kk] * tdz[ll] * ty[mm] * ty[nn] +
                                           tz[kk] * tz[ll] * tdy[mm] * ty[nn] +
                                           tz[kk] * tz[ll] * ty[mm] * tdy[nn]);
                        dr[isp] = dr[isp] +
                                  dg[22] * (ty[kk] * tz[ll] + tz[kk] * ty[ll]) *
                                      (ty[mm] * tz[nn] + tz[mm] * ty[nn]) +
                                  g[22] * ((tdy[kk] * tz[ll] + ty[kk] * tdz[ll] +
                                            tdz[kk] * ty[ll] + tz[kk] * tdy[ll]) *
                                               (ty[mm] * tz[nn] + tz[mm] * ty[nn]) +
                                           (ty[kk] * tz[ll] + tz[kk] * ty[ll]) *
                                               (tdy[mm] * tz[nn] + ty[mm] * tdz[nn] +
                                                tdz[mm] * ty[nn] + tz[mm] * tdy[nn]));
                    }
                }
            }
        }
    }
}

// ===========================================================================
void delri(double dg[23], int ni, int nj, double rr, double del1) {
    double term = (funcon_C::ev * del1) / (rr * funcon_C::a0 * funcon_C::a0);
    double da = dd[ni], db = dd[nj];
    double qa = qq[ni], qb = qq[nj];
    double aee = 0.25 * (1.0 / am[ni] + 1.0 / am[nj]) * (1.0 / am[ni] + 1.0 / am[nj]);
    double ee = -rr / std::pow(std::sqrt(rr * rr + aee), 3);
    dg[1] = term * ee;
    double dze = 0.0, qzze = 0.0, qxxe = 0.0;
    if (natorb[ni] <= 2 && natorb[nj] <= 2) return;
    if (natorb[ni] > 2) {
        double ade = 0.25 * (1.0 / ad[ni] + 1.0 / am[nj]) * (1.0 / ad[ni] + 1.0 / am[nj]);
        double aqe = 0.25 * (1.0 / aq[ni] + 1.0 / am[nj]) * (1.0 / aq[ni] + 1.0 / am[nj]);
        dze = (rr + da) / std::pow(std::sqrt((rr + da) * (rr + da) + ade), 3) -
                     (rr - da) / std::pow(std::sqrt((rr - da) * (rr - da) + ade), 3);
        qzze = -(rr + 2.0 * qa) / std::pow(std::sqrt((rr + 2.0 * qa) * (rr + 2.0 * qa) + aqe), 3) -
                      (rr - 2.0 * qa) / std::pow(std::sqrt((rr - 2.0 * qa) * (rr - 2.0 * qa) + aqe), 3) +
                      (2.0 * rr) / std::pow(std::sqrt(rr * rr + aqe), 3);
        qxxe = -(2.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qa * qa + aqe), 3) +
                      (2.0 * rr) / std::pow(std::sqrt(rr * rr + aqe), 3);
        dg[2] = -(term * dze) / 2.0;
        dg[3] = term * (ee + qzze / 4.0);
        dg[4] = term * (ee + qxxe / 4.0);
        if (natorb[nj] <= 2) return;
    }
    double aed = 0.25 * (1.0 / am[ni] + 1.0 / ad[nj]) * (1.0 / am[ni] + 1.0 / ad[nj]);
    double aeq = 0.25 * (1.0 / am[ni] + 1.0 / aq[nj]) * (1.0 / am[ni] + 1.0 / aq[nj]);
    double edz = (rr - db) / std::pow(std::sqrt((rr - db) * (rr - db) + aed), 3) -
                 (rr + db) / std::pow(std::sqrt((rr + db) * (rr + db) + aed), 3);
    double eqzz = -(rr - 2.0 * qb) / std::pow(std::sqrt((rr - 2.0 * qb) * (rr - 2.0 * qb) + aeq), 3) -
                  (rr + 2.0 * qb) / std::pow(std::sqrt((rr + 2.0 * qb) * (rr + 2.0 * qb) + aeq), 3) +
                  (2.0 * rr) / std::pow(std::sqrt(rr * rr + aeq), 3);
    double eqxx = -(2.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qb * qb + aeq), 3) +
                  (2.0 * rr) / std::pow(std::sqrt(rr * rr + aeq), 3);
    dg[5] = -(term * edz) / 2.0;
    dg[11] = term * (ee + eqzz / 4.0);
    dg[12] = term * (ee + eqxx / 4.0);
    if (natorb[ni] <= 2) return;

    double add = 0.25 * (1.0 / ad[ni] + 1.0 / ad[nj]) * (1.0 / ad[ni] + 1.0 / ad[nj]);
    double adq = 0.25 * (1.0 / ad[ni] + 1.0 / aq[nj]) * (1.0 / ad[ni] + 1.0 / aq[nj]);
    double aqd = 0.25 * (1.0 / aq[ni] + 1.0 / ad[nj]) * (1.0 / aq[ni] + 1.0 / ad[nj]);
    double aqq = 0.25 * (1.0 / aq[ni] + 1.0 / aq[nj]) * (1.0 / aq[ni] + 1.0 / aq[nj]);
    double dxdx = -(2.0 * rr) / std::pow(std::sqrt(rr * rr + (da - db) * (da - db) + add), 3) +
                  (2.0 * rr) / std::pow(std::sqrt(rr * rr + (da + db) * (da + db) + add), 3);
    double dzdz = -(rr + da - db) / std::pow(std::sqrt((rr + da - db) * (rr + da - db) + add), 3) -
                  (rr - da + db) / std::pow(std::sqrt((rr - da + db) * (rr - da + db) + add), 3) +
                  (rr - da - db) / std::pow(std::sqrt((rr - da - db) * (rr - da - db) + add), 3) +
                  (rr + da + db) / std::pow(std::sqrt((rr + da + db) * (rr + da + db) + add), 3);
    double dzqxx = 2.0 * (rr + da) / std::pow(std::sqrt((rr + da) * (rr + da) + 4.0 * qb * qb + adq), 3) -
                   2.0 * (rr - da) / std::pow(std::sqrt((rr - da) * (rr - da) + 4.0 * qb * qb + adq), 3) -
                   2.0 * (rr + da) / std::pow(std::sqrt((rr + da) * (rr + da) + adq), 3) +
                   2.0 * (rr - da) / std::pow(std::sqrt((rr - da) * (rr - da) + adq), 3);
    double qxxdz = 2.0 * (rr - db) / std::pow(std::sqrt((rr - db) * (rr - db) + 4.0 * qa * qa + aqd), 3) -
                   2.0 * (rr + db) / std::pow(std::sqrt((rr + db) * (rr + db) + 4.0 * qa * qa + aqd), 3) -
                   2.0 * (rr - db) / std::pow(std::sqrt((rr - db) * (rr - db) + aqd), 3) +
                   2.0 * (rr + db) / std::pow(std::sqrt((rr + db) * (rr + db) + aqd), 3);
    double dzqzz = (rr + da - 2.0 * qb) / std::pow(std::sqrt((rr + da - 2.0 * qb) * (rr + da - 2.0 * qb) + adq), 3) -
                   (rr - da - 2.0 * qb) / std::pow(std::sqrt((rr - da - 2.0 * qb) * (rr - da - 2.0 * qb) + adq), 3) +
                   (rr + da + 2.0 * qb) / std::pow(std::sqrt((rr + da + 2.0 * qb) * (rr + da + 2.0 * qb) + adq), 3) -
                   (rr - da + 2.0 * qb) / std::pow(std::sqrt((rr - da + 2.0 * qb) * (rr - da + 2.0 * qb) + adq), 3) +
                   2.0 * (rr - da) / std::pow(std::sqrt((rr - da) * (rr - da) + adq), 3) -
                   2.0 * (rr + da) / std::pow(std::sqrt((rr + da) * (rr + da) + adq), 3);
    double qzzdz = (rr + 2.0 * qa - db) / std::pow(std::sqrt((rr + 2.0 * qa - db) * (rr + 2.0 * qa - db) + aqd), 3) -
                   (rr + 2.0 * qa + db) / std::pow(std::sqrt((rr + 2.0 * qa + db) * (rr + 2.0 * qa + db) + aqd), 3) +
                   (rr - 2.0 * qa - db) / std::pow(std::sqrt((rr - 2.0 * qa - db) * (rr - 2.0 * qa - db) + aqd), 3) -
                   (rr - 2.0 * qa + db) / std::pow(std::sqrt((rr - 2.0 * qa + db) * (rr - 2.0 * qa + db) + aqd), 3) -
                   2.0 * (rr - db) / std::pow(std::sqrt((rr - db) * (rr - db) + aqd), 3) +
                   2.0 * (rr + db) / std::pow(std::sqrt((rr + db) * (rr + db) + aqd), 3);
    double qxxqxx = -(2.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * (qa - qb) * (qa - qb) + aqq), 3) -
                    (2.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * (qa + qb) * (qa + qb) + aqq), 3) +
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qa * qa + aqq), 3) +
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qb * qb + aqq), 3) -
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + aqq), 3);
    double qxxqyy = -(4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qa * qa + 4.0 * qb * qb + aqq), 3) +
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qa * qa + aqq), 3) +
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qb * qb + aqq), 3) -
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + aqq), 3);
    double qxxqzz = -2.0 * (rr - 2.0 * qb) / std::pow(std::sqrt((rr - 2.0 * qb) * (rr - 2.0 * qb) + 4.0 * qa * qa + aqq), 3) -
                   2.0 * (rr + 2.0 * qb) / std::pow(std::sqrt((rr + 2.0 * qb) * (rr + 2.0 * qb) + 4.0 * qa * qa + aqq), 3) +
                   2.0 * (rr - 2.0 * qb) / std::pow(std::sqrt((rr - 2.0 * qb) * (rr - 2.0 * qb) + aqq), 3) +
                   2.0 * (rr + 2.0 * qb) / std::pow(std::sqrt((rr + 2.0 * qb) * (rr + 2.0 * qb) + aqq), 3) +
                   (4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qa * qa + aqq), 3) -
                   (4.0 * rr) / std::pow(std::sqrt(rr * rr + aqq), 3);
    double qzzqxx = -2.0 * (rr + 2.0 * qa) / std::pow(std::sqrt((rr + 2.0 * qa) * (rr + 2.0 * qa) + 4.0 * qb * qb + aqq), 3) -
                   2.0 * (rr - 2.0 * qa) / std::pow(std::sqrt((rr - 2.0 * qa) * (rr - 2.0 * qa) + 4.0 * qb * qb + aqq), 3) +
                   2.0 * (rr + 2.0 * qa) / std::pow(std::sqrt((rr + 2.0 * qa) * (rr + 2.0 * qa) + aqq), 3) +
                   2.0 * (rr - 2.0 * qa) / std::pow(std::sqrt((rr - 2.0 * qa) * (rr - 2.0 * qa) + aqq), 3) +
                   (4.0 * rr) / std::pow(std::sqrt(rr * rr + 4.0 * qb * qb + aqq), 3) -
                   (4.0 * rr) / std::pow(std::sqrt(rr * rr + aqq), 3);
    double qzzqzz = -(rr + 2.0 * qa - 2.0 * qb) / std::pow(std::sqrt((rr + 2.0 * qa - 2.0 * qb) * (rr + 2.0 * qa - 2.0 * qb) + aqq), 3) -
                    (rr + 2.0 * qa + 2.0 * qb) / std::pow(std::sqrt((rr + 2.0 * qa + 2.0 * qb) * (rr + 2.0 * qa + 2.0 * qb) + aqq), 3) -
                    (rr - 2.0 * qa - 2.0 * qb) / std::pow(std::sqrt((rr - 2.0 * qa - 2.0 * qb) * (rr - 2.0 * qa - 2.0 * qb) + aqq), 3) -
                    (rr - 2.0 * qa + 2.0 * qb) / std::pow(std::sqrt((rr - 2.0 * qa + 2.0 * qb) * (rr - 2.0 * qa + 2.0 * qb) + aqq), 3) +
                    2.0 * (rr - 2.0 * qa) / std::pow(std::sqrt((rr - 2.0 * qa) * (rr - 2.0 * qa) + aqq), 3) +
                    2.0 * (rr + 2.0 * qa) / std::pow(std::sqrt((rr + 2.0 * qa) * (rr + 2.0 * qa) + aqq), 3) +
                    2.0 * (rr - 2.0 * qb) / std::pow(std::sqrt((rr - 2.0 * qb) * (rr - 2.0 * qb) + aqq), 3) +
                    2.0 * (rr + 2.0 * qb) / std::pow(std::sqrt((rr + 2.0 * qb) * (rr + 2.0 * qb) + aqq), 3) -
                    (4.0 * rr) / std::pow(std::sqrt(rr * rr + aqq), 3);
    double dxqxz = 2.0 * (rr - qb) / std::pow(std::sqrt((rr - qb) * (rr - qb) + (da - qb) * (da - qb) + adq), 3) -
                   2.0 * (rr + qb) / std::pow(std::sqrt((rr + qb) * (rr + qb) + (da - qb) * (da - qb) + adq), 3) -
                   2.0 * (rr - qb) / std::pow(std::sqrt((rr - qb) * (rr - qb) + (da + qb) * (da + qb) + adq), 3) +
                   2.0 * (rr + qb) / std::pow(std::sqrt((rr + qb) * (rr + qb) + (da + qb) * (da + qb) + adq), 3);
    double qxzdx = 2.0 * (rr + qa) / std::pow(std::sqrt((rr + qa) * (rr + qa) + (qa - db) * (qa - db) + aqd), 3) -
                   2.0 * (rr - qa) / std::pow(std::sqrt((rr - qa) * (rr - qa) + (qa - db) * (qa - db) + aqd), 3) -
                   2.0 * (rr + qa) / std::pow(std::sqrt((rr + qa) * (rr + qa) + (qa + db) * (qa + db) + aqd), 3) +
                   2.0 * (rr - qa) / std::pow(std::sqrt((rr - qa) * (rr - qa) + (qa + db) * (qa + db) + aqd), 3);
    double qxzqxz = -2.0 * (rr + qa - qb) / std::pow(std::sqrt((rr + qa - qb) * (rr + qa - qb) + (qa - qb) * (qa - qb) + aqq), 3) +
                    2.0 * (rr + qa + qb) / std::pow(std::sqrt((rr + qa + qb) * (rr + qa + qb) + (qa - qb) * (qa - qb) + aqq), 3) +
                    2.0 * (rr - qa - qb) / std::pow(std::sqrt((rr - qa - qb) * (rr - qa - qb) + (qa - qb) * (qa - qb) + aqq), 3) -
                    2.0 * (rr - qa + qb) / std::pow(std::sqrt((rr - qa + qb) * (rr - qa + qb) + (qa - qb) * (qa - qb) + aqq), 3) +
                    2.0 * (rr + qa - qb) / std::pow(std::sqrt((rr + qa - qb) * (rr + qa - qb) + (qa + qb) * (qa + qb) + aqq), 3) -
                    2.0 * (rr + qa + qb) / std::pow(std::sqrt((rr + qa + qb) * (rr + qa + qb) + (qa + qb) * (qa + qb) + aqq), 3) -
                    2.0 * (rr - qa - qb) / std::pow(std::sqrt((rr - qa - qb) * (rr - qa - qb) + (qa + qb) * (qa + qb) + aqq), 3) +
                    2.0 * (rr - qa + qb) / std::pow(std::sqrt((rr - qa + qb) * (rr - qa + qb) + (qa + qb) * (qa + qb) + aqq), 3);
    dg[6] = (term * dzdz) / 4.0;
    dg[7] = (term * dxdx) / 4.0;
    dg[8] = -term * (edz / 2.0 + qzzdz / 8.0);
    dg[9] = -term * (edz / 2.0 + qxxdz / 8.0);
    dg[10] = -(term * qxzdx) / 8.0;
    dg[13] = -term * (dze / 2.0 + dzqzz / 8.0);
    dg[14] = -term * (dze / 2.0 + dzqxx / 8.0);
    dg[15] = -(term * dxqxz) / 8.0;
    dg[16] = term * (ee + eqzz / 4.0 + qzze / 4.0 + qzzqzz / 16.0);
    dg[17] = term * (ee + eqzz / 4.0 + qxxe / 4.0 + qxxqzz / 16.0);
    dg[18] = term * (ee + eqxx / 4.0 + qzze / 4.0 + qzzqxx / 16.0);
    dg[19] = term * (ee + eqxx / 4.0 + qxxe / 4.0 + qxxqxx / 16.0);
    dg[20] = (term * qxzqxz) / 16.0;
    dg[21] = term * (ee + eqxx / 4.0 + qxxe / 4.0 + qxxqyy / 16.0);
    dg[22] = term * (qxxqxx - qxxqyy) / 32.0;
}

// ===========================================================================
void rotat(const std::vector<std::vector<double>>& coord, int i, int j, int ix,
           int idx, double rij, double del1) {
    double xd = coord[0][i] - coord[0][j];
    double yd = coord[1][i] - coord[1][j];
    double zd = coord[2][i] - coord[2][j];
    double rxy = std::sqrt(xd * xd + yd * yd);
    double ryz = std::sqrt(yd * yd + zd * zd);
    double rzx = std::sqrt(zd * zd + xd * xd);
    tx[1] = tx[2] = tx[3] = 0.0;
    ty[1] = ty[2] = ty[3] = 0.0;
    tz[1] = tz[2] = tz[3] = 0.0;
    tdx[1] = tdx[2] = tdx[3] = 0.0;
    tdy[1] = tdy[2] = tdy[3] = 0.0;
    tdz[1] = tdz[2] = tdz[3] = 0.0;
    if (rxy < 1.0e-4) {
        tx[3] = 1.0;
        if (zd < 0.0) tx[3] = -1.0;
        ty[2] = 1.0;
        tz[1] = tx[3];
        if (idx == 1) return;
        if (ix == 1) tdx[1] = 1.0 / rij;
        if (ix == 2) tdx[2] = 1.0 / rij;
        if (ix == 1) tdz[3] = -1.0 / rij;
        if (ix == 2) tdy[3] = -tx[3] / rij;
    } else if (ryz < 1.0e-4) {
        tx[1] = 1.0;
        if (xd < 0.0) tx[1] = -1.0;
        ty[2] = tx[1];
        tz[3] = 1.0;
        if (idx == 1) return;
        if (ix == 2) tdx[2] = 1.0 / rij;
        if (ix == 3) tdx[3] = 1.0 / rij;
        if (ix == 2) tdy[1] = -1.0 / rij;
        if (ix == 3) tdz[1] = -tx[1] / rij;
    } else if (rzx < 1.0e-4) {
        tx[2] = 1.0;
        if (yd < 0.0) tx[2] = -1.0;
        ty[1] = -tx[2];
        tz[3] = 1.0;
        if (idx == 1) return;
        if (ix == 1) tdx[1] = 1.0 / rij;
        if (ix == 3) tdx[3] = 1.0 / rij;
        if (ix == 1) tdy[2] = 1.0 / rij;
        if (ix == 3) tdz[2] = -tx[2] / rij;
    } else {
        tx[1] = xd / rij; tx[2] = yd / rij; tx[3] = zd / rij;
        tz[3] = rxy / rij;
        ty[1] = -tx[2] * (tx[1] >= 0 ? 1.0 : -1.0) / tz[3];
        ty[2] = std::fabs(tx[1] / tz[3]);
        ty[3] = 0.0;
        tz[1] = -tx[1] * tx[3] / tz[3];
        tz[2] = -tx[2] * tx[3] / tz[3];
        if (idx == 1) return;
        double term = del1 / (rij * rij);
        if (ix == 1) {
            tdx[1] = 1.0 / rij - tx[1] * term;
            tdx[2] = -tx[2] * term;
            tdx[3] = -tx[3] * term;
            tdz[3] = tx[1] / rxy - tz[3] * term;
        } else if (ix == 2) {
            tdx[1] = -tx[1] * term;
            tdx[2] = 1.0 / rij - tx[2] * term;
            tdx[3] = -tx[3] * term;
            tdz[3] = tx[2] / rxy - tz[3] * term;
        } else {
            tdx[1] = -tx[1] * term;
            tdx[2] = -tx[2] * term;
            tdx[3] = 1.0 / rij - tx[3] * term;
            tdz[3] = -tz[3] * term;
        }
        double tmp = tdz[3] / (tz[3] * tz[3]);
        tdy[1] = (-tdx[2] / tz[3]) + tx[2] * tmp;
        if (tx[1] < 0.0) tdy[1] = -tdy[1];
        tdy[2] = tdx[1] / tz[3] - tx[1] * tmp;
        if (tx[1] < 0.0) tdy[2] = -tdy[2];
        tdy[3] = 0.0;
        tdz[1] = tx[1] * tx[3] * tmp - (tx[3] * tdx[1] + tx[1] * tdx[1]) / tz[3];
        tdz[2] = tx[2] * tx[3] * tmp - (tx[3] * tdx[2] + tx[2] * tdx[3]) / tz[3];
    }
}

// ===========================================================================
void ders(int m, int n, double rr, double del1, double del2, double del3,
          int is, int iol) {
    using overlaps_C::ccc;
    using overlaps_C::zzz;
    double ss[7][7] = {};
    for (int i = 1; i <= 6; ++i) {
        for (int j = 1; j <= 6; ++j) {
            double apb = zzz[m][i] * zzz[n][j];
            double amb = zzz[m][i] + zzz[n][j];
            double adb = apb / amb;
            double adr = std::min(adb * rr, 35.0);
            double abn;
            switch (is) {
                default:
                    abn = -2.0 * adb * del1 / (funcon_C::a0 * funcon_C::a0);
                    break;
                case 2:
                    abn = -4.0 * adb * adb * del1 * del2 /
                          (std::sqrt(zzz[n][j]) * funcon_C::a0 * funcon_C::a0 * funcon_C::a0);
                    break;
                case 3:
                    abn = (2.0 * adb / (std::sqrt(zzz[n][j]) * funcon_C::a0)) *
                          (1.0 - 2.0 * adb * del1 * del1 /
                              (funcon_C::a0 * funcon_C::a0));
                    break;
                case 4:
                    abn = 4.0 * adb * adb * del1 * del2 /
                          (std::sqrt(zzz[m][i]) * funcon_C::a0 * funcon_C::a0 * funcon_C::a0);
                    break;
                case 5:
                    abn = -(2.0 * adb / (std::sqrt(zzz[m][i]) * funcon_C::a0)) *
                          (1.0 - 2.0 * adb * del1 * del1 /
                              (funcon_C::a0 * funcon_C::a0));
                    break;
                case 6:
                    abn = -(4.0 * adb * adb * del2 /
                            (std::sqrt(apb) * funcon_C::a0 * funcon_C::a0)) *
                          (1.0 - 2.0 * adb * del1 * del1 /
                              (funcon_C::a0 * funcon_C::a0));
                    break;
                case 7:
                    abn = 8.0 * adb * adb * adb * del1 * del2 * del3 /
                          (std::sqrt(apb) * funcon_C::a0 * funcon_C::a0 * funcon_C::a0 * funcon_C::a0);
                    break;
                case 8:
                    abn = -(8.0 * adb * adb * del1 /
                            (std::sqrt(apb) * funcon_C::a0 * funcon_C::a0)) *
                          (0.5 - adb * del2 * del2 / (funcon_C::a0 * funcon_C::a0));
                    break;
                case 9:
                    abn = -(8.0 * adb * adb * del1 /
                            (std::sqrt(apb) * funcon_C::a0 * funcon_C::a0)) *
                          (1.5 - adb * del1 * del1 / (funcon_C::a0 * funcon_C::a0));
                    break;
            }
            ss[i][j] = std::sqrt(std::pow(2.0 * std::sqrt(apb) / amb, 3)) *
                       std::exp(-adr) * abn;
        }
    }
    for (int i = 1; i <= 6; ++i) {
        double s = 0.0;
        for (int j = 1; j <= 6; ++j) s += ss[i][j] * ccc[m][i] * ccc[n][j];
        ds[iol] += s;
    }
}
