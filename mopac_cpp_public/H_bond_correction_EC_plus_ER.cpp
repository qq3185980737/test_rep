// H_bond_correction_EC_plus_ER.cpp — full C++ translation of
// EC_plus_ER (Korth/Pitonak/Rezac/Hobza transferable H-bond correction,
// attraction + repulsion terms).
#include "H_bond_correction_EC_plus_ER.h"

#include <cmath>

#include "bangle.h"
#include "common_arrays_C.h"
#include "dihed.h"
#include "funcon_C.h"
#include "H_bond_correction_bits.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace parameters_C;

double EC_plus_ER(int D, int H, int A, double q1, double q2,
                  double& EC, double& ER, int* d_list, int& nd_list) {
    static bool first = true;
    if (first) {
        first = false;
        dh2_a_parameters[1] = 1.48;   // Nitrogen
        dh2_a_parameters[2] = 1.56;   // Oxygen, generic
        dh2_a_parameters[3] = 1.55;   // Oxygen, acid
        dh2_a_parameters[4] = 0.96;   // Oxygen, peptide
        dh2_a_parameters[5] = 0.76;   // Oxygen, water
        dh2_a_parameters[6] = 0.85;   // Sulfur
    }
    EC = 0.0;
    ER = 0.0;
    double multiplier_a = 0.0;

    // First angle D-H-A
double angle = 0.0;
    bangle(coord, D, H, A, angle);
d_list[1] = D;
    d_list[2] = H;
    d_list[3] = A;
    nd_list = 3;
    double angle_cos = -std::cos(angle);
    if (angle_cos < 0.0) return 0.0;

    // target angles
    double torsion_check = 0.0;
    bool torsion_check_set = false;
    bool torsion_check_set2 = false;
    double angle2_shift = 0.0, angle2_shift_2 = 0.0, torsion_shift = 0.0;
    double torsion2 = 0.0;
    bool peptide = false;
    int C_of_CO = 0, N_of_HNCO = 0;

    if (nat[A] == 8 || nat[A] == 16) {
        if (nbonds[A] == 1) {
            angle2_shift = pi;
            angle2_shift_2 = pi / (180.0 / 120.0);
            torsion_shift = 0.0;
            torsion_check_set2 = true;
        } else {
            angle2_shift = pi / (180.0 / 109.48);
            angle2_shift_2 = angle2_shift;
            torsion_shift = pi / (180.0 / 54.74);
        }
        if (nat[A] == 8) {
            multiplier_a = 0.0;
            if (nbonds[A] == 2) {
                if (nat[ibonds[1][A]] == 1 && nat[ibonds[2][A]] == 1)
                    multiplier_a = dh2_a_parameters[5];   // Water
            } else if (nat[ibonds[1][A]] == 6) {
                C_of_CO = ibonds[1][A];
                if (nbonds[C_of_CO] == 3) {
                    int nH = 0, nC = 0, nO = 0, nN = 0;
                    for (int i = 1; i <= 3; ++i) {
                        if (nat[ibonds[i][C_of_CO]] == 1) ++nH;
                        if (nat[ibonds[i][C_of_CO]] == 6) ++nC;
                        if (nat[ibonds[i][C_of_CO]] == 7) {
                            ++nN;
                            N_of_HNCO = ibonds[i][C_of_CO];
                        }
                        if (nat[ibonds[i][C_of_CO]] == 8) ++nO;
                    }
                    if (nC + nH == 1 && nN == 1 && nO == 1) {
                        peptide = false;
                        // Check that H-N-C-O exists and is trans
                        nH = 0;
                        for (int i = 1; i <= nbonds[N_of_HNCO]; ++i) {
                            if (nat[ibonds[i][N_of_HNCO]] == 1) {
                                ++nH;
                                double sum = 0.0;
                                dihed(coord, A, C_of_CO, N_of_HNCO,
                                      ibonds[i][N_of_HNCO], sum);
                                sum = std::min(sum, 2.0 * pi - sum);
                                if (sum > 0.5 * pi) peptide = true;
                            }
                        }
                        if (peptide && nH > 0)
                            multiplier_a = dh2_a_parameters[4];  // peptide oxygen
                    } else if (nO == 2) {
                        multiplier_a = dh2_a_parameters[3];      // acid oxygen
                    }
                }
            }
            if (multiplier_a < 1e-20) multiplier_a = dh2_a_parameters[2];  // generic
        } else {
            multiplier_a = dh2_a_parameters[6];                 // sulfur
        }
    } else if (nat[A] == 7) {
        multiplier_a = dh2_a_parameters[1];                     // nitrogen
        if (nbonds[A] == 2) {
            angle2_shift = pi / (180.0 / 120.0);
            angle2_shift_2 = angle2_shift;
            torsion_shift = 0.0;
        } else {
            angle2_shift = pi / (180.0 / 109.48);
            angle2_shift_2 = angle2_shift;
            torsion_shift = pi / (180.0 / 54.74);
            torsion_check_set = true;                            // NR3 group
        }
    }

    // extrapolation between tetrahedral and planar NR3 group
    int R1 = 0, R2 = 0, R3 = 0;
    double sum_max = 0.0;
    if (nbonds[A] == 1) {
        R1 = ibonds[1][A];
        for (int ii = 1; ii <= nbonds[R1]; ++ii) {
            int i = ibonds[ii][R1];
            if (!connected_hb(i, H, 1000.0 * 1000.0)) return 0.0;
            if (sum_max < molkst_C::Rab) {
                R2 = i;
                sum_max = molkst_C::Rab;
            }
        }
        R3 = H;
    } else if (nbonds[A] == 2) {
        for (int ii = 1; ii <= nbonds[A]; ++ii) {
            int i = ibonds[ii][A];
            if (!connected_hb(i, H, 1000.0 * 1000.0)) return 0.0;
            if (sum_max < molkst_C::Rab) {
                sum_max = molkst_C::Rab;
                R1 = i;
            }
        }
        for (int ii = 1; ii <= nbonds[A]; ++ii) {
            int i = ibonds[ii][A];
            if (i == R1) continue;
            R2 = i;
        }
        R3 = H;
    } else if (nbonds[A] == 3) {
        for (int ii = 1; ii <= nbonds[A]; ++ii) {
            int i = ibonds[ii][A];
            if (!connected_hb(i, H, 1000.0 * 1000.0)) return 0.0;
            if (sum_max < molkst_C::Rab) {
                sum_max = molkst_C::Rab;
                R1 = i;
            }
        }
        sum_max = 0.0;
        for (int ii = 1; ii <= nbonds[A]; ++ii) {
            int i = ibonds[ii][A];
            if (i == R1) continue;
            if (!connected_hb(i, H, 1000.0 * 1000.0)) return 0.0;
            if (sum_max < molkst_C::Rab) {
                sum_max = molkst_C::Rab;
                R2 = i;
            }
        }
        for (int ii = 1; ii <= nbonds[A]; ++ii) {
            int i = ibonds[ii][A];
            if (i == R1 || i == R2) continue;
            R3 = i;
        }
        if (torsion_check_set) dihed(coord, R2, R1, A, R3, torsion2);
        torsion_check = torsion2;
        if (torsion_check < -pi) torsion_check += 2.0 * pi;
        if (torsion_check > pi)  torsion_check -= 2.0 * pi;
        if (torsion_check < 0.0) torsion_check = -pi - torsion_check;
        else                     torsion_check =  pi - torsion_check;
        double sum = torsion_check;
        if (sum < 0.0) sum = -sum;
        sum = 180.0 / pi * sum;
        torsion_shift = torsion_shift
            + pi / (180.0 / ((54.74 - sum) / 54.74 * 35.26));
        angle2_shift = angle2_shift
            - pi / (180.0 / ((54.74 - sum) / 54.74 * 19.48));
        angle2_shift_2 = angle2_shift;
    }
    if (R1 == 0) {
return 0.0; }
// second angle R1-A-H
    double angle2 = 0.0;
    bangle(coord, R1, A, H, angle2);
++nd_list;
    d_list[nd_list] = R1;
    double angle2_cos = std::cos(angle2_shift - angle2);
    double angle2_cos_2 = std::cos(angle2_shift_2 - angle2);
    if (angle2_cos_2 > angle2_cos) angle2_cos = angle2_cos_2;
    if (angle2_cos <= 0.0) return 0.0;

    // torsion angle, correction of NR3 torsion for through-bond case
double torsion_ref = 0.0;
    dihed(coord, R2, R1, A, H, torsion_ref);
++nd_list;
    d_list[nd_list] = R2;
    double torsion_correct = torsion_ref;
    if (torsion_correct < -pi) torsion_correct += 2.0 * pi;
    if (torsion_correct >  pi) torsion_correct -= 2.0 * pi;
    if (!torsion_check_set2 || std::fabs(torsion_correct) > 0.5 * pi) {
        if (torsion_correct < 0.0) torsion_correct = -pi - torsion_correct;
        else                       torsion_correct =  pi - torsion_correct;
    }

    double torsion_cos = 0.0;
    if (torsion_check < 0.0) {             // negative torsion occupied by -NR3
        double torsion = torsion_shift - torsion_correct;
        if (torsion < -pi) torsion += 2.0 * pi;
        if (torsion >  pi) torsion -= 2.0 * pi;
        torsion_cos = std::cos(torsion);
    } else if (torsion_check > 0.0) {      // positive torsion occupied by -NR3
        double torsion = -torsion_shift - torsion_correct;
        if (torsion < -pi) torsion += 2.0 * pi;
        if (torsion >  pi) torsion -= 2.0 * pi;
        torsion_cos = std::cos(torsion);
    } else {                               // planar -NR3 or general case
        double torsion = torsion_shift - torsion_correct;
        double torsion_2 = -torsion_shift - torsion_correct;
        if (torsion < -pi) torsion += 2.0 * pi;
        if (torsion >  pi) torsion -= 2.0 * pi;
        if (torsion_2 < -pi) torsion_2 += 2.0 * pi;
        if (torsion_2 >  pi) torsion_2 -= 2.0 * pi;
        torsion_cos = std::cos(torsion);
        double torsion_cos_2 = std::cos(torsion_2);
        if (torsion_cos_2 > torsion_cos) torsion_cos = torsion_cos_2;
    }
    if (!connected_hb(A, H, 1000.0 * 1000.0)) {
return 0.0; }
    double r = molkst_C::Rab;
    if (!connected_hb(R1, H, 1000.0 * 1000.0)) {
return 0.0; }
    if (torsion_check_set2 && r > molkst_C::Rab) torsion_cos = -1.0;
    if (torsion_cos <= 0.0) {
return 0.0; }
// distance cutoff (ad-hoc anti-discontinuity): smooth via truncation
    r = truncation(r, 1.80, 0.05);
    r = r / a0;
    const double expo = 3.0;       // "b" in equation 2
    const double rep_pre = 0.65;   // "c" in equation 2
    const double rep_exp = 5.0;    // "d" in equation 2
    double unit_part = fpc_9 * ev;
    double attraction = multiplier_a * q1 * q2 / std::pow(r, expo) * unit_part;
    double repulsion = rep_pre * std::pow(rep_exp, -r) * unit_part;
    double c = (attraction + repulsion) * angle_cos * angle2_cos * torsion_cos;
    EC = attraction * angle_cos * angle2_cos * torsion_cos;
    ER = repulsion * angle_cos * angle2_cos * torsion_cos;
    return c;
}
