// H_bond_correction_PM6_DH_Dispersion.cpp
#include "H_bond_correction_PM6_DH_Dispersion.h"
#include <cmath>
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

extern "C" bool connected(int, int, double);

static double C_[87] = {
    0.0, 0.16, 0.084, 0.0, 0.0, 5.79, 1.65, 1.11, 0.70,
    0.57, 0.45, 0.0, 0.0, 0.0, 0.0, 3.25, 5.79,
    5.97, 3.71, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 11.60, 4.47, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 25.80, 16.50, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};
static double R_[87] = {
    0.0, 156, 140, 0, 0, 180, 170, 155, 152,
    147, 154, 0, 0, 0, 0, 180, 180,
    175, 188, 0, 0, 0, 0, 0, 0,
    0, 0, 140, 0, 0, 0, 0, 0,
    0, 0, 185, 202, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 198, 216, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
};
static double N_[87] = {
    0.0, 0.80, 1.42, 0.0, 0.0, 2.16, 2.50, 2.82, 3.15,
    3.48, 3.81, 0.0, 0.0, 0.0, 0.0, 4.50, 4.80,
    5.10, 5.40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 2.90, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.00, 6.30, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 6.95, 7.25, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

double PM6_DH_Disp(int set_a, int nsa) {
    double alpha = 20.0, s = 1.04, cscale = 0.89;
    if (method_pm7) { alpha = 15.450118; s = 1.226593; cscale = 2.286419; }
    double E_disp = 0.0;
    for (int ii = 1; ii <= nsa; ++ii) {
        int i = (nsa == numat) ? ii : set_a;
        int ni = nat[i];
        if (ni > 86) continue;
        if (R_[ni] == 0.0 || C_[ni] == 0.0 || N_[ni] == 0.0) continue;
        int j_start = (nsa == numat) ? ii + 1 : 1;
        for (int j = j_start; j <= numat; ++j) {
            if (i == j) continue;
            int nj = nat[j];
            if (nj > 86) continue;
            double C6i = (ni == 6) ? (nbonds[i] == 4 ? 0.95 : 1.65) : C_[ni];
            double C6j = (nj == 6) ? (nbonds[j] == 4 ? 0.95 : 1.65) : C_[nj];
            if (R_[nj] == 0.0 || C_[nj] == 0.0 || N_[nj] == 0.0) continue;
            double C6 = 2.0 * std::pow(C6i*C6i*C6j*C6j*N_[ni]*N_[nj], 1.0/3.0) /
                        (std::pow(C6i*N_[nj]*N_[nj], 1.0/3.0) + std::pow(C6j*N_[ni]*N_[ni], 1.0/3.0));
            double R0 = (R_[ni]*R_[ni]*R_[ni] + R_[nj]*R_[nj]*R_[nj]) /
                        (R_[ni]*R_[ni] + R_[nj]*R_[nj]) / 1000.0 * 2.0;
            if (connected(i, j, 100.0*100.0)) {
                if (id == 0) {
                    double Rij = Rab * 0.1;
                    double damp = 1.0 / (1.0 + std::exp(-alpha * (Rij/(s*R0) - 1.0)));
                    double E_disp_tmp = C6 / std::pow(Rij, 6) * damp / (1000.0*4.184);
                    E_disp -= E_disp_tmp;
                }
            }
        }
    }
    return E_disp * cscale;
}

double PM6_DH_Dispersion(bool l_grad, std::vector<double>& dxyz) {
    double delta = 1e-5;
    double e = PM6_DH_Disp(0, numat);
    E_disp = e;
    if (l_grad) {
        for (int k = 1; k <= numat; ++k) {
            double sum2 = PM6_DH_Disp(k, 1);
            for (int i = 1; i <= 3; ++i) {
                coord[i - 1][k] += delta;
                double sum = PM6_DH_Disp(k, 1);
                double sum1 = (sum2 - sum) / delta;
                if (std::abs(sum1) < 50.0) dxyz[(k-1)*3 + i] -= sum1;
                coord[i - 1][k] -= delta;
            }
        }
    }
    return e;
}

// connected: true if atom_i and atom_j are within sqrt(criterion) of each
// other (translation of H_bond_correction_bits.F90 lines 1-49). Sets Rab
// (Bohr), Vab and cell_ijk as side effects.
extern "C" bool connected(int atom_i, int atom_j, double criterion) {
  using namespace molkst_C;
  using namespace common_arrays_C;
  if (Vab.size() < 4) Vab.resize(4);
  if (cell_ijk.size() < 4) cell_ijk.resize(4);
  if (id == 0) {
    Vab[1] = coord[0][atom_i] - coord[0][atom_j];
    Vab[2] = coord[1][atom_i] - coord[1][atom_j];
    Vab[3] = coord[2][atom_i] - coord[2][atom_j];
    Rab = Vab[1] * Vab[1] + Vab[2] * Vab[2] + Vab[3] * Vab[3];
  } else {
    Rab = 1.e8;
    for (int ii = -l11; ii <= l11; ++ii)
      for (int jj = -l21; jj <= l21; ++jj)
        for (int kk = -l31; kk <= l31; ++kk) {
          double v1 = coord[0][atom_i] - coord[0][atom_j] + tvec[1][1] * ii + tvec[1][2] * jj + tvec[1][3] * kk;
          double v2 = coord[1][atom_i] - coord[1][atom_j] + tvec[2][1] * ii + tvec[2][2] * jj + tvec[2][3] * kk;
          double v3 = coord[2][atom_i] - coord[2][atom_j] + tvec[3][1] * ii + tvec[3][2] * jj + tvec[3][3] * kk;
          double r2 = v1 * v1 + v2 * v2 + v3 * v3;
          if (r2 < Rab) {
            Rab = r2;
            Vab[1] = v1; Vab[2] = v2; Vab[3] = v3;
            cell_ijk[1] = ii; cell_ijk[2] = jj; cell_ijk[3] = kk;
          }
        }
  }
  bool conn = (Rab < criterion);
  if (conn) Rab = std::sqrt(Rab);
  return conn;
}
