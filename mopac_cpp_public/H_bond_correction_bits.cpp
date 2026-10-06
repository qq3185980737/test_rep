// H_bond_correction_bits.cpp — C++ translation.
#include "H_bond_correction_bits.h"
extern void mopend(const char*);

#include <algorithm>
#include <cmath>

#include "bangle.h"
#include "common_arrays_C.h"
#include "dihed.h"
#include "elemts_C.h"
#include "molkst_C.h"

#include <string>

double truncation(double R, double limit, double spread) {
    double a = limit - spread;
    double b = limit + spread;
    if (R < b) {
        if (R < a) return limit;
        return limit + (limit - a) / ((a - b) * (a - b)) * (R - a) * (R - a);
    }
    return R;
}

double distance_hb(int a, int b) {
    using namespace common_arrays_C;
    if (molkst_C::id == 0) {
        double dx = coord[0][a] - coord[0][b];
        double dy = coord[1][a] - coord[1][b];
        double dz = coord[2][a] - coord[2][b];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    double best = 1.0e6;
    for (int ik = -molkst_C::l1u; ik <= molkst_C::l1u; ++ik)
        for (int jk = -molkst_C::l2u; jk <= molkst_C::l2u; ++jk)
            for (int kl = -molkst_C::l3u; kl <= molkst_C::l3u; ++kl) {
                double c1 = coord[0][a] + tvec[1][1] * ik + tvec[1][2] * jk + tvec[1][3] * kl;
                double c2 = coord[1][a] + tvec[2][1] * ik + tvec[2][2] * jk + tvec[2][3] * kl;
                double c3 = coord[2][a] + tvec[3][1] * ik + tvec[3][2] * jk + tvec[3][3] * kl;
                double dx = c1 - coord[0][b];
                double dy = c2 - coord[1][b];
                double dz = c3 - coord[2][b];
                best = std::min(best, dx * dx + dy * dy + dz * dz);
            }
    return std::sqrt(best);
}

double angle_hb(int a, int b, int c) {
    double ang = 0.0;
    bangle(common_arrays_C::coord, a, b, c, ang);
    return ang;
}

double torsion_hb(int i, int j, int k, int l) {
    double tor = 0.0;
    dihed(common_arrays_C::coord, i, j, k, l, tor);
    return tor;
}

double bonding(int x, int y, const double* covrad) {
    return covrad[common_arrays_C::nat[x]] + covrad[common_arrays_C::nat[y]];
}

bool connected_hb(int atom_i, int atom_j, double criterion) {
    using namespace common_arrays_C;
    double& Rab = molkst_C::Rab;
    if ((int)Vab.size() < 4) Vab.resize(4, 0.0);
    if ((int)cell_ijk.size() < 4) cell_ijk.resize(4, 0);
    if (molkst_C::id == 0) {
        Vab[1] = coord[0][atom_i] - coord[0][atom_j];
        Vab[2] = coord[1][atom_i] - coord[1][atom_j];
        Vab[3] = coord[2][atom_i] - coord[2][atom_j];
        Rab = Vab[1]*Vab[1] + Vab[2]*Vab[2] + Vab[3]*Vab[3];
    } else {
        Rab = 1.0e8;
        for (int ii = -molkst_C::l11; ii <= molkst_C::l11; ++ii)
            for (int jj = -molkst_C::l21; jj <= molkst_C::l21; ++jj)
                for (int kk = -molkst_C::l31; kk <= molkst_C::l31; ++kk) {
                    double v1 = coord[0][atom_i] - coord[0][atom_j] + tvec[1][1]*ii + tvec[1][2]*jj + tvec[1][3]*kk;
                    double v2 = coord[1][atom_i] - coord[1][atom_j] + tvec[2][1]*ii + tvec[2][2]*jj + tvec[2][3]*kk;
                    double v3 = coord[2][atom_i] - coord[2][atom_j] + tvec[3][1]*ii + tvec[3][2]*jj + tvec[3][3]*kk;
                    double r = v1*v1 + v2*v2 + v3*v3;
                    if (r < Rab) {
                        Rab = r;
                        Vab[1]=v1; Vab[2]=v2; Vab[3]=v3;
                        cell_ijk[1]=ii; cell_ijk[2]=jj; cell_ijk[3]=kk;
                    }
                }
    }
    if (Rab < criterion) { Rab = std::sqrt(Rab); return true; }
    return false;
}

void find_XH_bonds(std::vector<int>& acc, int& nacc,
                   std::vector<int>& h_b, int& nhb) {
    using namespace common_arrays_C;
    double RAH; int is;
    if (molkst_C::method_pm6_dh_plus) { RAH = 1.4; is = 8; }
    else if (molkst_C::method_pm7)    { RAH = 1.4; is = 8; }
    else                              { RAH = 1.15; is = 16; }
    std::vector<char> used(molkst_C::numat + 1, 0);
    nacc = 0; nhb = 0;
    for (int i = 1; i <= molkst_C::numat; ++i) {
        int ni = nat[i];
        if (ni == 7 || ni == 8 || ni == is) {
            ++nacc;
            acc[nacc] = i;
            for (int j = 1; j <= molkst_C::numat; ++j) {
                if (nat[j] == 1 && !used[j]) {
                    if (connected_hb(i, j, RAH * RAH)) {
                        ++nhb;
                        h_b[nhb] = j;
                        used[j] = 1;
                    }
                }
            }
        }
    }
}

// External: parse float from text (full version lives in reada.F90).
double reada(const std::string& s, int pos);
// External stubs.
void web_message(int, const char*) {}

// NOTE: to_screen() lives in to_screen.cpp (real implementation). Do not
// define a duplicate here.

void find_H__Y_bonds(const std::vector<int>& acc_a, int nacc_a,
                     const std::vector<int>& acc_b, int nacc_b,
                     const std::vector<int>& bonding_a_h, int nb_a_h,
                     std::vector<int>& hblist1, std::vector<int>& hblist2,
                     std::vector<int>& hblist3, int max_h_bonds, int& nrpairs) {
    double RAH, cutoff;
    if (molkst_C::keywrd.find("PM6-DH+") != std::string::npos) {
        RAH = 1.4; cutoff = 10.0;
    } else if (molkst_C::method_pm7) {
        RAH = 1.4; cutoff = 7.0;
    } else {
        RAH = 1.15; cutoff = 7.0;
    }
    for (int ii = 1; ii <= nacc_a; ++ii) {
        int i = acc_a[ii];
        for (int jj = 1; jj <= nb_a_h; ++jj) {
            int j = bonding_a_h[jj];
            if (connected_hb(i, j, RAH * RAH)) {
                for (int kk = 1; kk <= nacc_b; ++kk) {
                    int k = acc_b[kk];
                    if (k != i) {
                        if (connected_hb(k, j, cutoff * cutoff)) {
                            if (angle_hb(k, j, i) > 0.5 * 3.14159265358979323846) {
                                int i1;
                                for (i1 = 1; i1 <= nrpairs; ++i1) {
                                    if (hblist2[i1] != j) continue;
                                    if (hblist1[i1] != k) continue;
                                    if (hblist3[i1] != i) continue;
                                    break;
                                }
                                if (i1 != nrpairs + 1) break;
                                for (i1 = 1; i1 <= nrpairs; ++i1) {
                                    if (hblist2[i1] != j) continue;
                                    if (hblist1[i1] != i) continue;
                                    if (hblist3[i1] != k) continue;
                                    break;
                                }
                                if (i1 != nrpairs + 1) break;
                                ++nrpairs;
                                if (nrpairs > max_h_bonds) {
                                    web_message(0, "PM6_DH_plus.html");
                                    mopend("The default array size for hydrogen bonds is too small");
                                    --nrpairs;
                                    return;
                                }
                                hblist3[nrpairs] = k;
                                hblist2[nrpairs] = j;
                                hblist1[nrpairs] = i;
                            }
                        }
                    }
                }
            }
        }
    }
}

void all_h_bonds(std::vector<int>& hblist1, std::vector<int>& hblist2,
                 std::vector<int>& hblist3, int max_h_bonds, int& nrpairs) {
    using namespace common_arrays_C;
    int numat = molkst_C::numat;
    if ((int)acceptor_a.size() < numat * 2 + 1)
        acceptor_a.resize(numat * 2 + 1, 0);
    if ((int)acceptor_b.size() < numat * 2 + 1)
        acceptor_b.resize(numat * 2 + 1, 0);
    if ((int)bonding_a_h.size() < numat * 2 + 1)
        bonding_a_h.resize(numat * 2 + 1, 0);
    if ((int)bonding_b_h.size() < numat * 2 + 1)
        bonding_b_h.resize(numat * 2 + 1, 0);

    int nacceptor_a, nbonding_a_h;
    find_XH_bonds(acceptor_a, nacceptor_a, bonding_a_h, nbonding_a_h);
    if (molkst_C::moperr) return;
    find_H__Y_bonds(acceptor_a, nacceptor_a, acceptor_a, nacceptor_a,
                    bonding_a_h, nbonding_a_h, hblist1, hblist2, hblist3,
                    max_h_bonds, nrpairs);
}

// ------------------------------------------------------------------ prt_hbonds
void prt_hbonds(int D, int H, int A, double energy) {
    using namespace common_arrays_C;
    using namespace elemts_C;
    using molkst_C::numcal;
    using molkst_C::numat;
    using molkst_C::P_Hbonds;
    using molkst_C::maxtxt;
    static int icalcn = -1;
    static bool prt_first = true;
    static double cutoff;
    if (icalcn != numcal) {
        icalcn = numcal;
        H_txt.assign(numat, "");
        H_energy.assign(numat, 0.0);
        prt_first = true;
        P_Hbonds = 0;
    }
    if (prt_first) {
        prt_first = false;
        std::string::size_type pos = molkst_C::keywrd.find(" DISP(");
        cutoff = (pos == std::string::npos) ? 0.0 : reada(molkst_C::keywrd, (int)pos + 1);
        cutoff = -std::fabs(cutoff);
        P_Hbonds = 0;
    }
    (void)cutoff;
    if (energy > -0.5) return;
    double sum1 = distance_hb(D, H);
    P_Hbonds = std::min(numat, P_Hbonds + 1);
    if (txtatm[D] == " ")
        txtatm[D] = "Atom No.: " + std::to_string(D) + "   " + elemnt[nat[D]];
    if (txtatm[H] == " ")
        txtatm[H] = "Atom No.: " + std::to_string(H) + "   " + elemnt[nat[H]];
    if (txtatm[A] == " ")
        txtatm[A] = "Atom No.: " + std::to_string(A) + "   " + elemnt[nat[A]];
    if (maxtxt == 0) {
        H_txt[P_Hbonds] = txtatm[D] + "  " + std::to_string(sum1)
                        + "    " + txtatm[H] + "    " + txtatm[A] + "  "
                        + std::to_string(energy) + " Kcal/mol";
    } else {
        H_txt[P_Hbonds] = "\"" + txtatm[D] + "\"  " + std::to_string(sum1)
                        + "  \"" + txtatm[H] + "\"  \"" + txtatm[A] + "\"  "
                        + std::to_string(energy) + " Kcal/mol";
    }
    H_energy[P_Hbonds] = energy;
}