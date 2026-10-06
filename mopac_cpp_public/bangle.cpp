// bangle.cpp — C++ translation of MOPAC 2016 "bangle.F90".

#include "bangle.h"

#include <algorithm>
#include <cmath>

#include "common_arrays_C.h"
#include "molkst_C.h"

void bangle(const std::vector<std::vector<double>>& xyz, int i, int j, int k,
            double& angle) {
    double d2ij, d2jk, d2ik;
    if (molkst_C::id == 0) {
        d2ij = (xyz[0][i] - xyz[0][j]) * (xyz[0][i] - xyz[0][j]) +
               (xyz[1][i] - xyz[1][j]) * (xyz[1][i] - xyz[1][j]) +
               (xyz[2][i] - xyz[2][j]) * (xyz[2][i] - xyz[2][j]);
        d2jk = (xyz[0][j] - xyz[0][k]) * (xyz[0][j] - xyz[0][k]) +
               (xyz[1][j] - xyz[1][k]) * (xyz[1][j] - xyz[1][k]) +
               (xyz[2][j] - xyz[2][k]) * (xyz[2][j] - xyz[2][k]);
        d2ik = (xyz[0][i] - xyz[0][k]) * (xyz[0][i] - xyz[0][k]) +
               (xyz[1][i] - xyz[1][k]) * (xyz[1][i] - xyz[1][k]) +
               (xyz[2][i] - xyz[2][k]) * (xyz[2][i] - xyz[2][k]);
    } else {
        d2ij = 1.0e8;
        d2jk = 1.0e8;
        d2ik = 1.0e8;
        for (int ii = -molkst_C::l11; ii <= molkst_C::l11; ++ii) {
            for (int jj = -molkst_C::l21; jj <= molkst_C::l21; ++jj) {
                for (int kk = -molkst_C::l31; kk <= molkst_C::l31; ++kk) {
                    double Vab1 = xyz[0][i] - xyz[0][j] +
                                  common_arrays_C::tvec[1][1] * ii +
                                  common_arrays_C::tvec[1][2] * jj +
                                  common_arrays_C::tvec[1][3] * kk;
                    double Vab2 = xyz[1][i] - xyz[1][j] +
                                  common_arrays_C::tvec[2][1] * ii +
                                  common_arrays_C::tvec[2][2] * jj +
                                  common_arrays_C::tvec[2][3] * kk;
                    double Vab3 = xyz[2][i] - xyz[2][j] +
                                  common_arrays_C::tvec[3][1] * ii +
                                  common_arrays_C::tvec[3][2] * jj +
                                  common_arrays_C::tvec[3][3] * kk;
                    double Rab = Vab1 * Vab1 + Vab2 * Vab2 + Vab3 * Vab3;
                    if (Rab < d2ij) d2ij = Rab;

                    Vab1 = xyz[0][k] - xyz[0][j] +
                           common_arrays_C::tvec[1][1] * ii +
                           common_arrays_C::tvec[1][2] * jj +
                           common_arrays_C::tvec[1][3] * kk;
                    Vab2 = xyz[1][k] - xyz[1][j] +
                           common_arrays_C::tvec[2][1] * ii +
                           common_arrays_C::tvec[2][2] * jj +
                           common_arrays_C::tvec[2][3] * kk;
                    Vab3 = xyz[2][k] - xyz[2][j] +
                           common_arrays_C::tvec[3][1] * ii +
                           common_arrays_C::tvec[3][2] * jj +
                           common_arrays_C::tvec[3][3] * kk;
                    Rab = Vab1 * Vab1 + Vab2 * Vab2 + Vab3 * Vab3;
                    if (Rab < d2jk) d2jk = Rab;

                    Vab1 = xyz[0][i] - xyz[0][k] +
                           common_arrays_C::tvec[1][1] * ii +
                           common_arrays_C::tvec[1][2] * jj +
                           common_arrays_C::tvec[1][3] * kk;
                    Vab2 = xyz[1][i] - xyz[1][k] +
                           common_arrays_C::tvec[2][1] * ii +
                           common_arrays_C::tvec[2][2] * jj +
                           common_arrays_C::tvec[2][3] * kk;
                    Vab3 = xyz[2][i] - xyz[2][k] +
                           common_arrays_C::tvec[3][1] * ii +
                           common_arrays_C::tvec[3][2] * jj +
                           common_arrays_C::tvec[3][3] * kk;
                    Rab = Vab1 * Vab1 + Vab2 * Vab2 + Vab3 * Vab3;
                    if (Rab < d2ik) d2ik = Rab;
                }
            }
        }
    }
    double xy = std::sqrt(d2ij * d2jk);
    if (xy < 1.0e-20) {
        angle = 0.0;
        return;
    }
    double temp = 0.5 * (d2ij + d2jk - d2ik) / xy;
    temp = std::min(1.0, temp);
    temp = std::max(-1.0, temp);
    angle = std::acos(temp);
}
