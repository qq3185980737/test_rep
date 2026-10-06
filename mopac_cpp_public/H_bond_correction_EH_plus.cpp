// H_bond_correction_EH_plus.cpp — C++ translation of "EH_plus.F90".
// Korth 2010 third-generation H-bond energy correction (PM7 / PM6-DH+).

#include "H_bond_correction_EH_plus.h"
#include "H_bond_correction_bits.h"

#include <cmath>

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"

// Fortran hblist(max_h_bonds, 10), column-major flat.
// Fortran hblist(i,col) is 1-based in both dimensions.
static inline int HB(const int* hblist, int max_h_bonds, int i, int col) {
    return hblist[(col - 1) * max_h_bonds + i];
}

double EH_plus(int i, const int* hblist, int max_h_bonds,
               const int* nrbondsa, const int* nrbondsb) {
    using namespace common_arrays_C;
    const double pi = funcon_C::pi;
    double angle_cos, angle, torsion_check, angle2_shift, angle2_shift_2,
           torsion_shift, angle2, torsion_check_bac, angle2_cos,
           angle2_cos_2, torsion_correct, torsion_value, torsion_cos,
           torsion_value_2, torsion_cos_2, angle2_cos_new,
           angle2_cos_2_new, torsion_cos_new, torsion_cos_2_new,
           scale_a, scale_nsp3, scale_nsp2, scale_osp3, scale_osp2, scale_b,
           scale_c, hb_dist, xc_dist, ha_dist, damping, hartree2kcal,
           XY_dist, short_o;
    double shortcut = 2.4, longcut = 7.0, covcut = 1.2;
    (void)angle; // declared for Fortran parity, not referenced
    bool torsion_check_set = false, torsion_check_set2 = false;

    if (molkst_C::method_pm7) {
        scale_nsp3 = -0.171271 * funcon_C::a0 * funcon_C::a0;
        scale_osp3 = -0.098822 * funcon_C::a0 * funcon_C::a0;
        scale_nsp2 = -0.171271 * funcon_C::a0 * funcon_C::a0;
    } else {
        scale_nsp3 = -0.16 * funcon_C::a0 * funcon_C::a0;
        scale_osp3 = -0.12 * funcon_C::a0 * funcon_C::a0;
        scale_nsp2 = scale_nsp3;
    }
    scale_osp2 = scale_osp3;
    hartree2kcal = funcon_C::ev * funcon_C::fpc_9;

    double E = 0.0;
    if (HB(hblist, max_h_bonds, i, 10) != -666) {
        // First angle: cos(pi - angle(D-H...A)).
        angle_cos = -std::cos(angle_hb(HB(hblist,max_h_bonds,i,1),
                                       HB(hblist,max_h_bonds,i,9),
                                       HB(hblist,max_h_bonds,i,5)));
        if (angle_cos <= 0) return 0.0;

        torsion_check = 0.0;
        torsion_check_set = false;
        torsion_check_set2 = false;

        // Target angles for donor (hblist(i,1)).
        int n1 = nat[HB(hblist,max_h_bonds,i,1)];
        if (n1 == 8) {  // oxygen
            if (nrbondsa[i] == 1) {
                angle2_shift = pi;
                angle2_shift_2 = pi / (180.0 / 120.0);
                torsion_shift = 0.0;
                torsion_check_set2 = true;
            } else {
                angle2_shift = pi / (180.0 / 109.48);
                angle2_shift_2 = angle2_shift;
                torsion_shift = pi / (180.0 / 54.74);
            }
        } else if (n1 == 7) {  // nitrogen
            if (nrbondsa[i] == 2) {
                angle2_shift = pi / (180.0 / 120.0);
                angle2_shift_2 = angle2_shift;
                torsion_shift = 0.0;
            } else {
                angle2_shift = pi / (180.0 / 109.48);
                angle2_shift_2 = angle2_shift;
                torsion_shift = pi / (180.0 / 54.74);
                torsion_check_set = true;
            }
        }

        if (torsion_check_set) {
            torsion_check = torsion_hb(HB(hblist,max_h_bonds,i,3),
                                       HB(hblist,max_h_bonds,i,2),
                                       HB(hblist,max_h_bonds,i,1),
                                       HB(hblist,max_h_bonds,i,4));
            if (torsion_check <= -pi) torsion_check += 2.0 * pi;
            if (torsion_check > pi) torsion_check -= 2.0 * pi;
            if (torsion_check < 0) torsion_check = -pi - torsion_check;
            else torsion_check = pi - torsion_check;
            torsion_check_bac = torsion_check;
            if (torsion_check < 0) torsion_check = -torsion_check;
            torsion_check = torsion_check * 180.0 / pi;
            torsion_shift += pi / (180.0 / ((54.74 - torsion_check) / 54.74 * 35.26));
            angle2_shift -= pi / (180.0 / ((54.74 - torsion_check) / 54.74 * 19.48));
            angle2_shift_2 = angle2_shift;
            torsion_check = torsion_check_bac;
        }

        angle2 = angle_hb(HB(hblist,max_h_bonds,i,2),
                          HB(hblist,max_h_bonds,i,1),
                          HB(hblist,max_h_bonds,i,9));
        angle2_cos = std::cos(angle2_shift - angle2);
        angle2_cos_2 = std::cos(angle2_shift_2 - angle2);
        if (angle2_cos_2 > angle2_cos) angle2_cos = angle2_cos_2;
        if (angle2_cos <= 0) return 0.0;

        torsion_correct = torsion_hb(HB(hblist,max_h_bonds,i,3),
                                      HB(hblist,max_h_bonds,i,2),
                                      HB(hblist,max_h_bonds,i,1),
                                      HB(hblist,max_h_bonds,i,9));
        if (torsion_correct <= -pi) torsion_correct += 2.0 * pi;
        if (torsion_correct > pi) torsion_correct -= 2.0 * pi;
        if (!torsion_check_set2 || std::abs(torsion_correct * 180.0 / pi) > 90.0) {
            if (torsion_correct < 0) torsion_correct = -pi - torsion_correct;
            else torsion_correct = pi - torsion_correct;
        }

        torsion_cos = 0.0;
        torsion_value = 0.0; torsion_value_2 = 0.0; torsion_cos_2 = 0.0;
        if (torsion_check < 0) {
            torsion_value = torsion_shift - torsion_correct;
            if (torsion_value <= -pi) torsion_value += 2.0 * pi;
            if (torsion_value > pi) torsion_value -= 2.0 * pi;
            torsion_cos = std::cos(torsion_value);
        } else if (torsion_check > 0) {
            torsion_value = -torsion_shift - torsion_correct;
            if (torsion_value <= -pi) torsion_value += 2.0 * pi;
            if (torsion_value > pi) torsion_value -= 2.0 * pi;
            torsion_cos = std::cos(torsion_value);
        } else {
            torsion_value = torsion_shift - torsion_correct;
            torsion_value_2 = -torsion_shift - torsion_correct;
            if (torsion_value <= -pi) torsion_value += 2.0 * pi;
            if (torsion_value > pi) torsion_value -= 2.0 * pi;
            if (torsion_value_2 <= -pi) torsion_value_2 += 2.0 * pi;
            if (torsion_value_2 > pi) torsion_value_2 -= 2.0 * pi;
            torsion_cos = std::cos(torsion_value);
            torsion_cos_2 = std::cos(torsion_value_2);
            if (torsion_cos_2 > torsion_cos) torsion_cos = torsion_cos_2;
        }
        if (distance_hb(HB(hblist,max_h_bonds,i,9), HB(hblist,max_h_bonds,i,1)) >
            distance_hb(HB(hblist,max_h_bonds,i,9), HB(hblist,max_h_bonds,i,2)) &&
            torsion_check_set2) torsion_cos = 0.0;
        if (HB(hblist,max_h_bonds,i,3) == HB(hblist,max_h_bonds,i,4) ||
            HB(hblist,max_h_bonds,i,9) == HB(hblist,max_h_bonds,i,4)) torsion_cos = 1.0;
        if (torsion_cos < 0) return 0.0;

        // Target angles for acceptor (hblist(i,5)).
        torsion_check = 0.0;
        torsion_check_set = false;
        torsion_check_set2 = false;
        int n5 = nat[HB(hblist,max_h_bonds,i,5)];
        if (n5 == 8) {
            if (nrbondsb[i] == 1) {
                angle2_shift = pi;
                angle2_shift_2 = pi / (180.0 / 120.0);
                torsion_shift = 0.0;
                torsion_check_set2 = true;
            } else {
                angle2_shift = pi / (180.0 / 109.48);
                angle2_shift_2 = angle2_shift;
                torsion_shift = pi / (180.0 / 54.74);
            }
        } else if (n5 == 7) {
            if (nrbondsb[i] == 2) {
                angle2_shift = pi / (180.0 / 120.0);
                angle2_shift_2 = angle2_shift;
                torsion_shift = 0.0;
            } else {
                angle2_shift = pi / (180.0 / 109.48);
                angle2_shift_2 = angle2_shift;
                torsion_shift = pi / (180.0 / 54.74);
                torsion_check_set = true;
            }
        }

        if (torsion_check_set) {
            torsion_check = torsion_hb(HB(hblist,max_h_bonds,i,7),
                                       HB(hblist,max_h_bonds,i,6),
                                       HB(hblist,max_h_bonds,i,5),
                                       HB(hblist,max_h_bonds,i,8));
            if (torsion_check <= -pi) torsion_check += 2.0 * pi;
            if (torsion_check > pi) torsion_check -= 2.0 * pi;
            if (torsion_check < 0) torsion_check = -pi - torsion_check;
            else torsion_check = pi - torsion_check;
            torsion_check_bac = torsion_check;
            if (torsion_check < 0) torsion_check = -torsion_check;
            torsion_check = torsion_check * 180.0 / pi;
            torsion_shift += pi / (180.0 / ((54.74 - torsion_check) / 54.74 * 35.26));
            angle2_shift -= pi / (180.0 / ((54.74 - torsion_check) / 54.74 * 19.48));
            angle2_shift_2 = angle2_shift;
            torsion_check = torsion_check_bac;
        }

        angle2 = angle_hb(HB(hblist,max_h_bonds,i,6),
                          HB(hblist,max_h_bonds,i,5),
                          HB(hblist,max_h_bonds,i,9));
        angle2_cos_new = std::cos(angle2_shift - angle2);
        angle2_cos_2_new = std::cos(angle2_shift_2 - angle2);
        if (angle2_cos_2_new > angle2_cos_new) angle2_cos_new = angle2_cos_2_new;
        if (angle2_cos_new <= 0) return 0.0;

        torsion_correct = torsion_hb(HB(hblist,max_h_bonds,i,7),
                                     HB(hblist,max_h_bonds,i,6),
                                     HB(hblist,max_h_bonds,i,5),
                                     HB(hblist,max_h_bonds,i,9));
        if (torsion_correct <= -pi) torsion_correct += 2.0 * pi;
        if (torsion_correct > pi) torsion_correct -= 2.0 * pi;
        if (!torsion_check_set2 || std::abs(torsion_correct * 180.0 / pi) > 90.0) {
            if (torsion_correct < 0) torsion_correct = -pi - torsion_correct;
            else torsion_correct = pi - torsion_correct;
        }

        torsion_cos_new = 0.0;
        torsion_value = 0.0; torsion_value_2 = 0.0; torsion_cos_2_new = 0.0;
        if (torsion_check < 0) {
            torsion_value = torsion_shift - torsion_correct;
            if (torsion_value <= -pi) torsion_value += 2.0 * pi;
            if (torsion_value > pi) torsion_value -= 2.0 * pi;
            torsion_cos_new = std::cos(torsion_value);
        } else if (torsion_check > 0) {
            torsion_value = -torsion_shift - torsion_correct;
            if (torsion_value <= -pi) torsion_value += 2.0 * pi;
            if (torsion_value > pi) torsion_value -= 2.0 * pi;
            torsion_cos_new = std::cos(torsion_value);
        } else {
            torsion_value = torsion_shift - torsion_correct;
            torsion_value_2 = -torsion_shift - torsion_correct;
            if (torsion_value <= -pi) torsion_value += 2.0 * pi;
            if (torsion_value > pi) torsion_value -= 2.0 * pi;
            if (torsion_value_2 <= -pi) torsion_value_2 += 2.0 * pi;
            if (torsion_value_2 > pi) torsion_value_2 -= 2.0 * pi;
            torsion_cos_new = std::cos(torsion_value);
            torsion_cos_2_new = std::cos(torsion_value_2);
            if (torsion_cos_2_new > torsion_cos_new) torsion_cos_new = torsion_cos_2_new;
        }
        if (distance_hb(HB(hblist,max_h_bonds,i,9), HB(hblist,max_h_bonds,i,5)) >
            distance_hb(HB(hblist,max_h_bonds,i,9), HB(hblist,max_h_bonds,i,6)) &&
            torsion_check_set2) torsion_cos_new = 0.0;
        if (HB(hblist,max_h_bonds,i,7) == HB(hblist,max_h_bonds,i,8) ||
            HB(hblist,max_h_bonds,i,9) == HB(hblist,max_h_bonds,i,8)) torsion_cos_new = 1.0;
        torsion_cos_new = std::abs(torsion_cos_new);

        // Scale factors.
        if (n1 == 7) {
            scale_a = (nrbondsa[i] >= 3) ? scale_nsp3 : scale_nsp2;
        } else {
            scale_a = (nrbondsa[i] >= 2) ? scale_osp3 : scale_osp2;
        }
        if (n5 == 7) {
            scale_b = (nrbondsb[i] >= 3) ? scale_nsp3 : scale_nsp2;
        } else {
            scale_b = (nrbondsb[i] >= 2) ? scale_osp3 : scale_osp2;
        }
        scale_c = (scale_a + scale_b) / 2.0;

        ha_dist = distance_hb(HB(hblist,max_h_bonds,i,9), HB(hblist,max_h_bonds,i,1));
        hb_dist = distance_hb(HB(hblist,max_h_bonds,i,9), HB(hblist,max_h_bonds,i,5));
        xc_dist = std::min(ha_dist, hb_dist);

        if (molkst_C::method_pm7) {
            XY_dist = std::max(ha_dist, hb_dist) - xc_dist;
            if (XY_dist > 0.5) {
                damping = 1.0 - 1.0 / (1.0 + std::exp(-60.0 * (xc_dist / covcut - 1.0)));
            } else {
                damping = 1.0;
            }
            xc_dist = distance_hb(HB(hblist,max_h_bonds,i,1), HB(hblist,max_h_bonds,i,5));
            damping = damping / (1.0 + std::exp(-100.0 * (xc_dist / shortcut - 1.0)));
            damping = damping * (1.0 - 1.0 / (1.0 + std::exp(-10.0 * (xc_dist / longcut - 1.0))));
            XY_dist = distance_hb(HB(hblist,max_h_bonds,i,1), HB(hblist,max_h_bonds,i,5));

            E = scale_c / (XY_dist * XY_dist) * angle_cos * angle_cos *
                (1.0 - (1.0 - angle2_cos * torsion_cos * angle2_cos_new * torsion_cos_new) *
                        (1.0 - angle2_cos * torsion_cos * angle2_cos_new * torsion_cos_new)) *
                hartree2kcal * damping;
            if (n1 == 8 && n5 == 8) {
                short_o = -2.5 * std::exp(-80.0 * std::pow(std::max(XY_dist - 2.67, 0.0), 2)) *
                        std::pow(angle_cos, 4);
                E += short_o;
            }
        } else {
            damping = 1.0 - 1.0 / (1.0 + std::exp(-60.0 * (xc_dist / covcut - 1.0)));
            xc_dist = distance_hb(HB(hblist,max_h_bonds,i,1), HB(hblist,max_h_bonds,i,5));
            damping = damping / (1.0 + std::exp(-100.0 * (xc_dist / shortcut - 1.0)));
            damping = damping * (1.0 - 1.0 / (1.0 + std::exp(-10.0 * (xc_dist / longcut - 1.0))));
            double XY = distance_hb(HB(hblist,max_h_bonds,i,1), HB(hblist,max_h_bonds,i,5));
            E = scale_c / (XY * XY) * angle_cos * angle_cos *
                angle2_cos * angle2_cos * torsion_cos * torsion_cos *
                angle2_cos_new * angle2_cos_new * torsion_cos_new * torsion_cos_new *
                hartree2kcal * damping;
        }
    }
    return E;
}
