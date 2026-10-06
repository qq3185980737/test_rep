// dftd3_bits.cpp — C++ translation of MOPAC 2016 "dftd3_bits.F90" + "gdisp.F90".
#include "dftd3_bits.h"
#include "d3_r0ab_data.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// ---------------------------------------------------------------- anagrdc6
void anagrdc6(int, int, int,
              const std::vector<double>& cn,
              const std::vector<std::vector<double>>& dcn2,
              const std::vector<std::vector<std::vector<double>>>& dcn3,
              const std::vector<int>& nat,
              const std::vector<int>& mxc,
              int iat, int jat, int kat,
              const c6ab_t& c6ab,
              std::vector<double>& anag) {
    const double k3 = -4.0;
    std::vector<double> dterm2(4, 0.0), dterm3(4, 0.0);
    if (iat == kat) {
        for (int k = 1; k <= 3; ++k) {
            dterm2[k] = dcn2[k][iat];
            dterm3[k] = dcn3[k][iat][jat];
        }
    } else {
        for (int k = 1; k <= 3; ++k) {
            dterm2[k] = dcn3[k][kat][iat];
            dterm3[k] = dcn3[k][kat][jat];
        }
    }
    double zaehler = 0.0, nenner = 0.0;
    std::vector<double> dzaehler(4, 0.0), dnenner(4, 0.0);
    for (int i = 1; i <= mxc[nat[iat]]; ++i)
        for (int j = 1; j <= mxc[nat[jat]]; ++j) {
            double term3 = c6ab[nat[iat]][nat[jat]][i][j][2] - cn[jat];
            double term2 = c6ab[nat[iat]][nat[jat]][i][j][1] - cn[iat];
            double term1 = std::exp(k3 * (term2 * term2 + term3 * term3));
            zaehler += c6ab[nat[iat]][nat[jat]][i][j][0] * term1;
            nenner += term1;
            double term4 = term1 * k3 * 2.0;
            for (int k = 1; k <= 3; ++k) {
                dzaehler[k] += c6ab[nat[iat]][nat[jat]][i][j][0] * term4 * (term2 * dterm2[k] + term3 * dterm3[k]);
                dnenner[k] += term4 * (term2 * dterm2[k] + term3 * dterm3[k]);
            }
        }
    if (nenner > 0) {
        double term4 = 1.0 / (nenner * nenner);
        for (int k = 1; k <= 3; ++k)
            anag[k] = (dzaehler[k] * nenner - dnenner[k] * zaehler) * term4;
    } else {
        for (int k = 1; k <= 3; ++k) anag[k] = 0.0;
    }
}

// ------------------------------------------------------------------- d3limit
void d3limit(int iat, int jat, int& iadr, int& jadr) {
    iadr = 1; jadr = 1;
    while (iat > 100) { iat -= 100; ++iadr; }
    while (jat > 100) { jat -= 100; ++jadr; }
}

// ---------------------------------------------------------------------- lin
int lin(int i1, int i2) {
    int idum1 = std::max(i1, i2);
    int idum2 = std::min(i1, i2);
    return idum2 + idum1 * (idum1 - 1) / 2;
}

// --------------------------------------------------------------------- esym
std::string esym(int i) {
    static const char* ELEMNT[95] = {"", "h ", "he",
        "li", "be", "b ", "c ", "n ", "o ", "f ", "ne",
        "na", "mg", "al", "si", "p ", "s ", "cl", "ar",
        "k ", "ca", "sc", "ti", "v ", "cr", "mn", "fe", "co", "ni", "cu",
        "zn", "ga", "ge", "as", "se", "br", "kr",
        "rb", "sr", "y ", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag",
        "cd", "in", "sn", "sb", "te", "i ", "xe",
        "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy",
        "ho", "er", "tm", "yb", "lu", "hf", "ta", "w ", "re", "os", "ir", "pt",
        "au", "hg", "tl", "pb", "bi", "po", "at", "rn",
        "fr", "ra", "ac", "th", "pa", "u ", "np", "pu"};
    if (i >= 1 && i <= 94) return std::string(ELEMNT[i]);
    return std::string("  ");
}

// --------------------------------------------------------------------- eabh
double eabh(int n, int A, int B, int H,
            const std::vector<std::vector<double>>& xyz,
            double shortcut, double cab) {
    double rab2 = (xyz[1][A]-xyz[1][B])*(xyz[1][A]-xyz[1][B])
                + (xyz[2][A]-xyz[2][B])*(xyz[2][A]-xyz[2][B])
                + (xyz[3][A]-xyz[3][B])*(xyz[3][A]-xyz[3][B]);
    double rab = std::sqrt(rab2);
    double d2ij = (xyz[1][A]-xyz[1][H])*(xyz[1][A]-xyz[1][H])
                + (xyz[2][A]-xyz[2][H])*(xyz[2][A]-xyz[2][H])
                + (xyz[3][A]-xyz[3][H])*(xyz[3][A]-xyz[3][H]);
    double d2jk = (xyz[1][H]-xyz[1][B])*(xyz[1][H]-xyz[1][B])
                + (xyz[2][H]-xyz[2][B])*(xyz[2][H]-xyz[2][B])
                + (xyz[3][H]-xyz[3][B])*(xyz[3][H]-xyz[3][B]);
    double d2ik = rab2;
    double xy = std::sqrt(d2ij * d2jk + 1e-14);
    double cosabh = 0.5 * (d2ij + d2jk - d2ik) / xy;
    double aterm = 0.125 * std::pow(cosabh - 1.0, 4);
    double xm = 0.5 * (xyz[1][A] + xyz[1][B]);
    double ym = 0.5 * (xyz[2][A] + xyz[2][B]);
    double zm = 0.5 * (xyz[3][A] + xyz[3][B]);
    double rhm = std::sqrt((xyz[1][H]-xm)*(xyz[1][H]-xm)
                         + (xyz[2][H]-ym)*(xyz[2][H]-ym)
                         + (xyz[3][H]-zm)*(xyz[3][H]-zm));
    const double alp = 20.0, longcut = 7.50388;
    double dampm = 1.0 - 1.0 / (1.0 + std::exp(-alp * (rhm / rab - 1.0)));
    double damps = 1.0 / (1.0 + std::exp(-alp * (rab / shortcut - 1.0)));
    double dampl = 1.0 - 1.0 / (1.0 + std::exp(-alp * (rab / longcut - 1.0)));
    (void)n;
    return -dampl * damps * dampm * aterm * cab / (rab2 * rab2);
}

// -------------------------------------------------------------------- hbpar
int hbpar(int elem) {
    if (elem == 7)  return 1;
    if (elem == 8)  return 2;
    if (elem == 9)  return 3;
    if (elem == 15) return 4;
    if (elem == 16) return 5;
    if (elem == 17) return 6;
    return 0;
}

// -------------------------------------------------------------------- getc6
void getc6(int maxc, int max_elem, const c6ab_t& c6ab, const std::vector<int>& mxc,
           int iat, int jat, double nci, double ncj, double& c6) {
    (void)maxc; (void)max_elem;
    double c6mem = -1e99;
    double rsum = 0.0, csum = 0.0;
    c6 = 0.0;
    for (int i = 1; i <= mxc[iat]; ++i)
        for (int j = 1; j <= mxc[jat]; ++j) {
            double c6v = c6ab[iat][jat][i][j][1];
            if (c6v > 0.0) {
                c6mem = c6v;
                double cn1 = c6ab[iat][jat][i][j][2];
                double cn2 = c6ab[iat][jat][i][j][3];
                double r = (cn1 - nci) * (cn1 - nci) + (cn2 - ncj) * (cn2 - ncj);
                double tmp1 = std::exp(-4.0 * r);
                rsum += tmp1;
                csum += tmp1 * c6v;
            }
        }
    c6 = (rsum > 0) ? csum / rsum : c6mem;
}

// -------------------------------------------------------------------- edisp
void edisp(int max_elem, int maxc, int n, const std::vector<std::vector<double>>& xyz,
           const std::vector<int>& nat, const c6ab_t& c6ab, const std::vector<int>& mxc,
           const std::vector<double>& r2r4, const std::vector<std::vector<double>>& r0ab,
           const std::vector<double>& rcov, double rs6, double rs8, double alp6, double alp8,
           double& e6, double& e8) {
    using funcon_C::a0;
    double rthr = 15.0;
    double rthr2 = (rthr / a0) * (rthr / a0);
    e6 = 0.0; e8 = 0.0;
    std::vector<double> cn(n + 1, 0.0);
    ncoord(n, rcov, nat, xyz, cn);
    for (int iat = 1; iat <= n - 1; ++iat)
        for (int jat = iat + 1; jat <= n; ++jat) {
            double dx = xyz[1][iat] - xyz[1][jat];
            double dy = xyz[2][iat] - xyz[2][jat];
            double dz = xyz[3][iat] - xyz[3][jat];
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > rthr2) continue;
            double r = std::sqrt(r2);
            double rr = r0ab[nat[jat]][nat[iat]] / r;
            double tmp = rs6 * rr;
            double damp6 = 1.0 / (1.0 + 6.0 * std::pow(tmp, alp6));
            tmp = rs8 * rr;
            double damp8 = 1.0 / (1.0 + 6.0 * std::pow(tmp, alp8));
            double c6 = 0.0;
            getc6(maxc, max_elem, c6ab, mxc, nat[iat], nat[jat], cn[iat], cn[jat], c6);
            double r6 = r2 * r2 * r2;
            double r8 = r6 * r2;
            double c8 = 3.0 * c6 * r2r4[nat[iat]] * r2r4[nat[jat]];
            e6 += c6 * damp6 / r6;
            e8 += c8 * damp8 / r8;
        }
}

// ------------------------------------------------------------------ hbsimple
namespace {
// Diagnostic printer for hydrogen bonds; full version lives in
// H_bond_correction_bits.F90 (prt_hbonds). Only invoked with " PRT " + " DISP(".
void prt_hbonds_stub(int B, int H, int A, double energy) {
    std::printf("HBOND: A=%d H=%d B=%d  E=%12.4f kcal/mol\n", A, H, B, energy);
    (void)B;
}
}

void hbsimple(int n, const std::vector<int>& at, std::vector<std::vector<double>>& xyz,
              double hbscale, double& energy, bool l_grad, std::vector<std::vector<double>>& g) {
    using molkst_C::N_Hbonds;
    using molkst_C::keywrd;
    double scalehb[7] = {0, hbscale, hbscale, hbscale, hbscale, hbscale, hbscale};
    bool prt = (keywrd.find(" PRT ") != std::string::npos && keywrd.find(" DISP(") != std::string::npos);
    double thr = 250.0;
    double r0ab[7][7] = {{0}};
    r0ab[1][1]=4.95580619390096; r0ab[2][1]=4.69521322049756; r0ab[2][2]=4.68973278168166;
    r0ab[3][1]=4.51361038303964; r0ab[3][2]=4.44293461883933; r0ab[3][3]=4.34561357746278;
    r0ab[4][1]=5.75705004325974; r0ab[4][2]=5.42861568913110; r0ab[4][3]=5.22773805284199;
    r0ab[4][4]=6.61725321392116;
    r0ab[5][1]=5.49040984445458; r0ab[5][2]=5.44335574472965; r0ab[5][3]=5.16462109512098;
    r0ab[5][4]=6.27011084754718; r0ab[5][5]=6.25631558644032;
    r0ab[6][1]=5.28216218036620; r0ab[6][2]=5.18862031695641; r0ab[6][3]=5.07977251245376;
    r0ab[6][4]=6.03124949908998; r0ab[6][5]=5.95698288505716; r0ab[6][6]=5.86684309279986;
    for (int i = 1; i <= 6; ++i)
        for (int j = 1; j <= i; ++j)
            r0ab[j][i] = r0ab[i][j];
    const double shortcutscale = 1.0;
    for (int i = 1; i <= 6; ++i)
        for (int j = 1; j <= 6; ++j)
            r0ab[i][j] *= shortcutscale;
    int nhb = 0;
    for (int i = 1; i <= n - 1; ++i)
        if (at[i] == 7 || at[i] == 8 || at[i] == 9 || at[i] == 15 || at[i] == 16 || at[i] == 17)
            for (int j = i + 1; j <= n; ++j)
                if (at[j] == 7 || at[j] == 8 || at[j] == 9 || at[j] == 15 || at[j] == 16 || at[j] == 17) {
                    double rab = (xyz[1][i]-xyz[1][j])*(xyz[1][i]-xyz[1][j])
                               + (xyz[2][i]-xyz[2][j])*(xyz[2][i]-xyz[2][j])
                               + (xyz[3][i]-xyz[3][j])*(xyz[3][i]-xyz[3][j]);
                    if (rab < thr)
                        for (int k = 1; k <= n; ++k)
                            if (at[k] == 1) ++nhb;
                }
    std::vector<std::vector<double>> hbs(4, std::vector<double>(std::max(nhb, 1), 0.0));
    nhb = 0;
    for (int i = 1; i <= n - 1; ++i)
        if (at[i] == 7 || at[i] == 8 || at[i] == 9 || at[i] == 15 || at[i] == 16 || at[i] == 17)
            for (int j = i + 1; j <= n; ++j)
                if (at[j] == 7 || at[j] == 8 || at[j] == 9 || at[j] == 15 || at[j] == 16 || at[j] == 17) {
                    double rab = (xyz[1][i]-xyz[1][j])*(xyz[1][i]-xyz[1][j])
                               + (xyz[2][i]-xyz[2][j])*(xyz[2][i]-xyz[2][j])
                               + (xyz[3][i]-xyz[3][j])*(xyz[3][i]-xyz[3][j]);
                    if (rab < thr)
                        for (int k = 1; k <= n; ++k)
                            if (at[k] == 1) {
                                ++nhb;
                                hbs[1][nhb] = (double)i;
                                hbs[2][nhb] = (double)j;
                                hbs[3][nhb] = (double)k;
                            }
                }
    N_Hbonds = 0;
    energy = 0.0;
    for (int m = 1; m <= nhb; ++m) {
        int A = (int)std::lround(hbs[1][m]);
        int B = (int)std::lround(hbs[2][m]);
        int H = (int)std::lround(hbs[3][m]);
        int typa = hbpar(at[A]);
        int typb = hbpar(at[B]);
        double shortcut = r0ab[typa][typb];
        double cab = 0.50 * (scalehb[typa] + scalehb[typb]);
        double e = eabh(n, A, B, H, xyz, shortcut, cab);
        if (e < -1.0 / 627.51) ++N_Hbonds;
        if (prt) prt_hbonds_stub(B, H, A, e * 627.51);
        energy += e;
    }
    if (l_grad) {
        for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= n; ++j) g[i][j] = 0.0;
        double step = 1e-5;
        for (int m = 1; m <= nhb; ++m) {
            int A = (int)std::lround(hbs[1][m]);
            int B = (int)std::lround(hbs[2][m]);
            int H = (int)std::lround(hbs[3][m]);
            int typa = hbpar(at[A]);
            int typb = hbpar(at[B]);
            double cab = 0.50 * (scalehb[typa] + scalehb[typb]);
            double shortcut = r0ab[hbpar(at[A])][hbpar(at[B])];
            for (int j = 1; j <= 3; ++j) {
                xyz[j][A] += step;
                double er = eabh(n, A, B, H, xyz, shortcut, cab);
                xyz[j][A] -= step * 2.0;
                double el = eabh(n, A, B, H, xyz, shortcut, cab);
                xyz[j][A] += step;
                g[j][A] += (er - el) / (2.0 * step);
            }
            for (int j = 1; j <= 3; ++j) {
                xyz[j][B] += step;
                double er = eabh(n, A, B, H, xyz, shortcut, cab);
                xyz[j][B] -= step * 2.0;
                double el = eabh(n, A, B, H, xyz, shortcut, cab);
                xyz[j][B] += step;
                g[j][B] += (er - el) / (2.0 * step);
            }
            for (int j = 1; j <= 3; ++j) {
                xyz[j][H] += step;
                double er = eabh(n, A, B, H, xyz, shortcut, cab);
                xyz[j][H] -= step * 2.0;
                double el = eabh(n, A, B, H, xyz, shortcut, cab);
                xyz[j][H] += step;
                g[j][H] += (er - el) / (2.0 * step);
            }
        }
    }
}

// ------------------------------------------------------------------ setr0ab
void setr0ab(int max_elem, double autoang, std::vector<std::vector<double>>& r) {
    int k = 0;
    for (int i = 1; i <= max_elem; ++i)
        for (int j = 1; j <= i; ++j) {
            ++k;
            double v = d3_r0ab_data[k - 1] / autoang;
            r[i][j] = v;
            r[j][i] = v;
        }
}
