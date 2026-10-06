// H_bond_correction_PM6_DH_type.cpp — C++ translation.
// PM6-DH+/PM7 hydrogen-bond energy correction driver.

#include "H_bond_correction_PM6_DH_type.h"
#include "H_bond_correction_bits.h"
#include "H_bond_correction_EH_plus.h"
#include "H_bond_correction_EC_plus_ER.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

// External.
void chrge(const std::vector<double>& p, std::vector<double>& vector_out);
void prt_hbonds(int D, int H, int A, double energy);
void mopend(const char*);
void to_screen(const std::string&);

// Covalent radii (1..94), scaled by 4/3 once.
static const double covrad_raw[95] = {
    0,
    0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54,
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39,
    1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26,
    1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57,
    1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,
    1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58,
    1.52, 1.53, 1.54, 1.55
};
static double covrad[95];
static bool covrad_scaled = false;

double PM6_DH_H_bond_corrections(bool l_grad, bool prt) {
    using namespace common_arrays_C;
    using namespace molkst_C;
    if (!covrad_scaled) {
        for (int j = 1; j <= 94; ++j) covrad[j] = 4.0 / 3.0 * covrad_raw[j];
        covrad_scaled = true;
    }

    static int icalcn = -1;
    static int max_h_bonds = 0;
    if (icalcn != numcal) {
        max_h_bonds = 0;
        for (int j = 1; j <= numat; ++j)
            if (nat[j] == 8 || nat[j] == 7) ++max_h_bonds;
        max_h_bonds = (id == 0) ? max_h_bonds * 110 : max_h_bonds * 250;
        icalcn = numcal;
    }
    if (max_h_bonds == 0) {
        E_hb = 0.0;
        N_Hbonds = 0;
        return 0.0;
    }

    hblist.assign(max_h_bonds * 10 + 1, 0);
    std::vector<int> nrbondsa(max_h_bonds + 1, 0);
    std::vector<int> nrbondsb(max_h_bonds + 1, 0);
    int nrpairs = 0;

    // all_h_bonds expects 1-based vectors sized (max_h_bonds).
    std::vector<int> h1(max_h_bonds + 1, 0), h2(max_h_bonds + 1, 0), h3(max_h_bonds + 1, 0);
    if (method_pm6_dh_plus || method_pm7) {
        all_h_bonds(h1, h2, h3, max_h_bonds, nrpairs);
        // Copy into global hblist columns 1,9,5.
        for (int i = 1; i <= nrpairs; ++i) {
            hblist[(1 - 1) * max_h_bonds + i] = h1[i];
            hblist[(9 - 1) * max_h_bonds + i] = h2[i];
            hblist[(5 - 1) * max_h_bonds + i] = h3[i];
        }
        bool l_h = true;
        setup_DH_Plus(nrpairs, nrbondsa.data(), nrbondsb.data(), l_h, covrad);

    } else {
        all_h_bonds(h1, h2, h3, max_h_bonds, nrpairs);
        for (int i = 1; i <= nrpairs; ++i) {
            hblist[(1 - 1) * max_h_bonds + i] = h1[i];
            hblist[(2 - 1) * max_h_bonds + i] = h2[i];
            hblist[(3 - 1) * max_h_bonds + i] = h3[i];
        }
    }

    if (method_pm6_dh2 || method_pm6_dh2x) {
        std::vector<double> vec(numat + 1, 0.0);
        chrge(p, vec);
        for (int i = 1; i <= numat; ++i) {
            int j = nat[i];
            q[i] = parameters_C::tore[j] - vec[i];
        }
    }

    N_Hbonds = 0;
    double E_hb_acc = 0.0;
    std::vector<int> d_list(10, 0), d_l(10, 0);
    double delta = 1.0e-5;

    for (int ii = 1; ii <= nrpairs; ++ii) {
        int D, H, A;
        double sum, EC = 0, ER = 0;
        int nd_list = 0;
        if (method_pm6_dh_plus || method_pm7) {
            H = hblist[(9 - 1) * max_h_bonds + ii];
            A = hblist[(1 - 1) * max_h_bonds + ii];
            D = hblist[(5 - 1) * max_h_bonds + ii];
            sum = EH_plus(ii, hblist.data(), max_h_bonds, nrbondsa.data(), nrbondsb.data());
            E_hb_acc += sum;
            for (int i = 1; i <= 9; ++i) {
                int j = hblist[(i - 1) * max_h_bonds + ii];
                if (j > 0) {
                    int l;
                    for (l = 1; l <= nd_list; ++l)
                        if (d_list[l] == j) break;
                    if (l > nd_list) {
                        ++nd_list;
                        d_list[nd_list] = j;
                    }
                }
            }
        } else {
            H = hblist[(2 - 1) * max_h_bonds + ii];
            A = hblist[(3 - 1) * max_h_bonds + ii];
            D = hblist[(1 - 1) * max_h_bonds + ii];
            sum = EC_plus_ER(D, H, A, q[H], q[A], EC, ER, d_list.data(), nd_list);
            E_hb_acc += EC + ER;
        }
        if (sum < -1.0) ++N_Hbonds;
        if (prt) prt_hbonds(D, H, A, sum);

        if (l_grad && sum < -0.01) {
            for (int j = 1; j <= nd_list; ++j) {
                int k = d_list[j];
                int iii = l123 * (k - 1);
                if (connected_hb(H, k, 64.0)) {
                    int kkkk = (l3u - cell_ijk[3]) +
                               (2 * l3u + 1) * (l2u - cell_ijk[2] +
                               (2 * l2u + 1) * (l1u - cell_ijk[1]));
                    int i_cell = iii + kkkk;
                    for (int i = 1; i <= 3; ++i) {
                        coord[i-1][k] += delta;
                        double sum1;
                        if (method_pm6_dh_plus || method_pm7) {
                            sum1 = EH_plus(ii, hblist.data(), max_h_bonds,
                                           nrbondsa.data(), nrbondsb.data());
                        } else {
                            int ltmp = 0;
                            sum1 = EC_plus_ER(D, H, A, q[H], q[A], EC, ER, d_l.data(), ltmp);
                        }
                        sum1 = (sum1 - sum) / delta;
                        if (std::abs(sum1) < 50.0) {
                            dxyz[i_cell * 3 + i] += sum1;
                        }
                        coord[i-1][k] -= delta;
                    }
                }
            }
        }
    }
    E_hb = E_hb_acc;
    return E_hb_acc;
}

void setup_DH_Plus(int nrpairs, int* nrbondsa, int* nrbondsb,
                   bool& l_h_bonds, double* covrad_in) {
    using namespace common_arrays_C;
    using namespace molkst_C;
    l_h_bonds = true;
    bool hbs1_ok = true, hbs2_ok = true;

    for (int i = 1; i <= nrpairs; ++i) {
        nrbondsa[i] = 0;
        nrbondsb[i] = 0;
        int bondlist[6][4] = {};
        // max_h_bonds mirrors the driver's allocation; recompute from numat.
        int mb = 0;
        for (int j = 1; j <= numat; ++j)
            if (nat[j] == 8 || nat[j] == 7) ++mb;
        mb = (id == 0) ? mb * 110 : mb * 250;
        auto HB = [&](int r, int col) -> int& {
            return hblist[(col - 1) * mb + r];
        };

        for (int j = 1; j <= numat; ++j) {
            double xa_dist = distance_hb(j, HB(i, 1));
            double xb_dist = distance_hb(j, HB(i, 5));
            if (xa_dist < bonding(j, HB(i, 1), covrad_in) && HB(i, 1) != j) {
                ++nrbondsa[i];
                bondlist[nrbondsa[i]][1] = j;
                if (nrbondsa[i] == 5) {
                    int i1 = HB(i, 1);
                    double mx = 0; int l = 1;
                    for (int k = 1; k <= 5; ++k) {
                        double s1 = distance_hb(i1, bondlist[k][1]);
                        if (mx < s1) { mx = s1; l = k; }
                    }
                    for (int k = l; k <= 4; ++k) bondlist[k][1] = bondlist[k + 1][1];
                    nrbondsa[i] = 4;
                }
            }
            if (xb_dist < bonding(j, HB(i, 5), covrad_in) && HB(i, 5) != j) {
                ++nrbondsb[i];
                bondlist[nrbondsb[i]][2] = j;
                if (nrbondsb[i] == 5) {
                    int i1 = HB(i, 5);
                    double mx = 0; int l = 1;
                    for (int k = 1; k <= 5; ++k) {
                        double s1 = distance_hb(i1, bondlist[k][2]);
                        if (mx < s1) { mx = s1; l = k; }
                    }
                    for (int k = l; k <= 4; ++k) bondlist[k][2] = bondlist[k + 1][2];
                    nrbondsb[i] = 4;
                }
            }
        }
        // 1-3 / 1-4 check.
        for (int j = 1; j <= nrbondsa[i]; ++j) {
            if (bondlist[j][1] == HB(i, 5)) HB(i, 10) = -666;
            for (int k = 1; k <= nrbondsb[i]; ++k) {
                if (bondlist[j][1] == bondlist[k][2] && bondlist[k][2] != HB(i, 9))
                    HB(i, 10) = -666;
            }
        }
        if (HB(i, 10) == -666) continue;

        // hbs1 (donor side).
        if (nrbondsa[i] == 3 || nrbondsa[i] == 4) {
            double old_dist = -1;
            for (int k = 1; k <= nrbondsa[i]; ++k) {
                double xh = distance_hb(bondlist[k][1], HB(i, 9));
                if (xh > old_dist) { old_dist = xh; HB(i, 2) = bondlist[k][1]; }
            }
            old_dist = -1;
            for (int k = 1; k <= nrbondsa[i]; ++k) {
                double xh = distance_hb(bondlist[k][1], HB(i, 9));
                if (xh > old_dist && bondlist[k][1] != HB(i, 2)) {
                    old_dist = xh; HB(i, 3) = bondlist[k][1];
                }
            }
            old_dist = -1;
            for (int k = 1; k <= nrbondsa[i]; ++k) {
                double xh = distance_hb(bondlist[k][1], HB(i, 9));
                if (xh > old_dist && bondlist[k][1] != HB(i, 2) && bondlist[k][1] != HB(i, 3)) {
                    old_dist = xh; HB(i, 4) = bondlist[k][1];
                }
            }
        } else if (nrbondsa[i] == 2) {
            double old_dist = -1;
            for (int k = 1; k <= nrbondsa[i]; ++k) {
                double xh = distance_hb(bondlist[k][1], HB(i, 9));
                if (xh > old_dist) { old_dist = xh; HB(i, 2) = bondlist[k][1]; }
            }
            for (int k = 1; k <= nrbondsa[i]; ++k)
                if (bondlist[k][1] != HB(i, 2)) HB(i, 3) = bondlist[k][1];
            HB(i, 4) = (distance_hb(HB(i, 1), HB(i, 9)) < bonding(HB(i,1),HB(i,9),covrad_in))
                       ? HB(i, 9) : HB(i, 1);
        } else if (nrbondsa[i] == 1) {
            HB(i, 2) = bondlist[1][1];
            int nrbondsc = 0;
            for (int k = 1; k <= numat; ++k) {
                double xc = distance_hb(k, HB(i, 2));
                if (xc < bonding(k, HB(i, 2), covrad_in) && HB(i, 2) != k) {
                    ++nrbondsc;
                    bondlist[nrbondsc][3] = k;
                }
            }
            double old_dist = -1;
            for (int k = 1; k <= nrbondsc; ++k) {
                double xh = distance_hb(bondlist[k][3], HB(i, 9));
                if (xh > old_dist) { old_dist = xh; HB(i, 3) = bondlist[k][3]; }
            }
            HB(i, 4) = (distance_hb(HB(i, 1), HB(i, 9)) < bonding(HB(i,1),HB(i,9),covrad_in))
                       ? HB(i, 9) : HB(i, 1);
        } else if (nrbondsa[i] == 0) {
            if (distance_hb(HB(i, 1), HB(i, 9)) < bonding(HB(i,1),HB(i,9),covrad_in)) {
                HB(i, 2) = HB(i, 9); HB(i, 3) = HB(i, 9); HB(i, 4) = HB(i, 9);
            } else {
                HB(i, 2) = HB(i, 1); HB(i, 3) = HB(i, 1); HB(i, 4) = HB(i, 1);
            }
        } else {
            hbs1_ok = false;
        }

        // hbs2 (acceptor side).
        if (nrbondsb[i] == 3 || nrbondsb[i] == 4) {
            double old_dist = -1;
            for (int k = 1; k <= nrbondsb[i]; ++k) {
                double xh = distance_hb(bondlist[k][2], HB(i, 9));
                if (xh > old_dist) { old_dist = xh; HB(i, 6) = bondlist[k][2]; }
            }
            old_dist = -1;
            for (int k = 1; k <= nrbondsb[i]; ++k) {
                double xh = distance_hb(bondlist[k][2], HB(i, 9));
                if (xh > old_dist && bondlist[k][2] != HB(i, 6)) {
                    old_dist = xh; HB(i, 7) = bondlist[k][2];
                }
            }
            old_dist = -1;
            for (int k = 1; k <= nrbondsb[i]; ++k) {
                double xh = distance_hb(bondlist[k][2], HB(i, 9));
                if (xh > old_dist && bondlist[k][2] != HB(i, 6) && bondlist[k][2] != HB(i, 7)) {
                    old_dist = xh; HB(i, 8) = bondlist[k][2];
                }
            }
        } else if (nrbondsb[i] == 2) {
            double old_dist = -1;
            for (int k = 1; k <= nrbondsb[i]; ++k) {
                double xh = distance_hb(bondlist[k][2], HB(i, 9));
                if (xh > old_dist) { old_dist = xh; HB(i, 6) = bondlist[k][2]; }
            }
            for (int k = 1; k <= nrbondsb[i]; ++k)
                if (bondlist[k][2] != HB(i, 6)) HB(i, 7) = bondlist[k][2];
            HB(i, 8) = (distance_hb(HB(i, 5), HB(i, 9)) < bonding(HB(i,5),HB(i,9),covrad_in))
                       ? HB(i, 9) : HB(i, 5);
        } else if (nrbondsb[i] == 1) {
            HB(i, 6) = bondlist[1][2];
            int nrbondsc = 0;
            for (int k = 1; k <= numat; ++k) {
                double xc = distance_hb(k, HB(i, 6));
                if (xc < bonding(k, HB(i, 6), covrad_in) && HB(i, 6) != k) {
                    ++nrbondsc;
                    bondlist[nrbondsc][3] = k;
                }
            }
            double old_dist = -1;
            for (int k = 1; k <= nrbondsc; ++k) {
                double xh = distance_hb(bondlist[k][3], HB(i, 9));
                if (xh > old_dist) { old_dist = xh; HB(i, 7) = bondlist[k][3]; }
            }
            HB(i, 8) = (distance_hb(HB(i, 5), HB(i, 9)) < bonding(HB(i,5),HB(i,9),covrad_in))
                       ? HB(i, 9) : HB(i, 5);
        } else if (nrbondsb[i] == 0) {
            if (distance_hb(HB(i, 5), HB(i, 9)) < bonding(HB(i,5),HB(i,9),covrad_in)) {
                HB(i, 6) = HB(i, 9); HB(i, 7) = HB(i, 9); HB(i, 8) = HB(i, 9);
            } else {
                HB(i, 6) = HB(i, 5); HB(i, 7) = HB(i, 5); HB(i, 8) = HB(i, 5);
            }
        } else {
            hbs2_ok = false;
        }
    }
    if (!hbs1_ok || !hbs2_ok) l_h_bonds = false;
}
