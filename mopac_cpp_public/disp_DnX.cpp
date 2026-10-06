// disp_DnX.cpp — C++ translation of MOPAC 2016 "disp_DnX.F90".
//
// disp_DnX: Rezac & Hobza halogen-bonding correction for PM6-DH2X /
// PM6-DH+.  print_post_scf_corrections: prints the D3/H-bond correction
// summary at the end of a job.
//
// distance / connected live in H_bond_correction_bits.F90 and are
// translated there as distance_hb / connected_hb; reada lives in reada.F90
// (M02).  We bind to those real translations.

#include "disp_DnX.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "H_bond_correction_bits.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "reada.h"

using namespace common_arrays_C;
using namespace molkst_C;

double disp_DnX(bool l_grad) {
    static double a_X[4][4] = {}, b_X[4][4] = {};
    static int iX[101] = {}, jOorN[101] = {};
    static bool first = true;
    if (first) {
        first = false;
        if (method_PM6_DH2X) {
            a_X[1][1]=1.0489e12; a_X[2][1]=1.0226e5; a_X[3][1]=1.2751e12;
            a_X[1][2]=4.6783e8; a_X[2][2]=9.6021e3; a_X[3][2]=6.0912e5;
            b_X[1][1]=-9.946; b_X[2][1]=-3.236; b_X[3][1]=-9.534;
            b_X[1][2]=-6.867; b_X[2][2]=-2.9; b_X[3][2]=-4.154;
        } else {
            a_X[1][1]=1.049e12; a_X[2][1]=5.560e4; a_X[3][1]=5.237e8;
            a_X[1][2]=1.871e9; a_X[2][2]=2.160e4; a_X[3][2]=2.436e6; a_X[3][3]=1.051e6;
            b_X[1][1]=-9.95; b_X[2][1]=-3.04; b_X[3][1]=-6.77;
            b_X[1][2]=-7.44; b_X[2][2]=-3.30; b_X[3][2]=-4.71; b_X[3][3]=-3.82;
        }
        iX[17]=1; iX[35]=2; iX[53]=3;
        jOorN[7]=1; jOorN[8]=2; jOorN[16]=3;
    }
    double sum = 0.0;
    for (int i = 1; i <= numat; ++i) {
        if (nat[i]!=17 && nat[i]!=35 && nat[i]!=53) continue;  // crude, but fast
        int k = nat[i];
        for (int j = 1; j <= numat; ++j) {
            if (nat[j]!=7 && nat[j]!=8 && nat[j]!=16) continue;
            if (k != 53 && nat[j] == 16) continue;             // If sulfur, only iodine
            int l = nat[j];
            double Rab = distance_hb(i, j);
            double sum2 = a_X[iX[k]][jOorN[l]];
            double sum3 = b_X[iX[k]][jOorN[l]];
            sum += sum2 * std::exp(sum3 * Rab);
            if (l_grad) {
                if (connected_hb(i, j, 8.0 * 8.0)) {
                    // kkkk is the cell that atom j is in, relative to atom i
                    int iii = l123 * (i - 1);
                    int jjj = l123 * (j - 1);
                    int kkkk = (l3u - cell_ijk[3]) +
                        (2 * l3u + 1) * (l2u - cell_ijk[2] +
                        (2 * l2u + 1) * (l1u - cell_ijk[1]));
                    int i_cell = iii + kkkk;
                    int j_cell = jjj - kkkk;
                    double fact = sum2 * sum3 * std::exp(sum3 * Rab) / Rab;
                    for (int m = 1; m <= 3; ++m) {
                        dxyz[i_cell * 3 + m] += Vab[m] * fact;
                        dxyz[j_cell * 3 + m] -= Vab[m] * fact;
                    }
                }
            }
        }
    }
    return sum;
}

void print_post_scf_corrections() {
    double sum1 = 0.0, sum = 0.0;
    int k = 0, j = 0;
    if (keywrd.find(" DISP(") != std::string::npos) {
        std::printf("\n%47s\n", " List of hydrogen bonds found");
        std::printf("%3s%12s%16s%11s%23s%17s\n", "No.", "Donor",
                    "R(D-H)", "Hydrogen", "Acceptor", "H-bond energy");
        size_t p = keywrd.find(" DISP(");
        sum1 = -std::fabs(reada(keywrd, (int)p + 5));
        for (;;) {
            sum = 0.0;
            for (int i = 1; i <= P_Hbonds; ++i)
                if (sum > H_energy[i]) {
                    sum = H_energy[i];
                    j = i;
                }
            if (sum > sum1) break;
            ++k;
            std::printf("%5d %s\n", k, H_txt[j].c_str());
            H_energy[j] = 10.0;
        }
    }
    if (keywrd.find("0SCF") != std::string::npos) {
        std::printf("\n%10s%17.5f%s\n", "DISPERSION ENERGY       =", E_disp, " KCAL/MOL");
        std::printf("%10s%17.5f%s\n\n", "H-BOND ENERGY           =", E_hb, " KCAL/MOL");
    }
}
