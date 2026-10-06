// ijbo.cpp — C++ translation of ijbo.F90 (MOPAC 2016).
#include "ijbo.h"

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "overlaps_C.h"

// Given atom numbers ii and jj:
//   IJBO = -1 if the atoms are separated by more than SQRT(CUTOF1), else
//   IJBO = -2 if the atoms are separated by more than SQRT(CUTOF2), else
//   IJBO = N, where "N+1" is the address of the starting location in the
//           H, P, and F arrays.
int ijbo(int ii, int jj) {
    if (MOZYME_C::lijbo) {
        return MOZYME_C::nijbo[ii][jj];
    }
    double r = (common_arrays_C::coord[0][ii] - common_arrays_C::coord[0][jj]) *
                   (common_arrays_C::coord[0][ii] - common_arrays_C::coord[0][jj]) +
               (common_arrays_C::coord[1][ii] - common_arrays_C::coord[1][jj]) *
                   (common_arrays_C::coord[1][ii] - common_arrays_C::coord[1][jj]) +
               (common_arrays_C::coord[2][ii] - common_arrays_C::coord[2][jj]) *
                   (common_arrays_C::coord[2][ii] - common_arrays_C::coord[2][jj]);
    if (r > overlaps_C::cutof1) return -1;
    if (r > overlaps_C::cutof2) return -2;
    //  use lookup array to find ijbo
    int ind_i, ind_j;
    if (ii > jj) {
        ind_i = ii;
        ind_j = jj;
    } else {
        ind_i = jj;
        ind_j = ii;
    }
    int il = MOZYME_C::iij[ind_i];
    int iu = MOZYME_C::numij[ind_i];
    int j = (il + iu + 1) / 2;
    int jm1 = 0, jm2 = 0;
    for (;;) {
        int ik = MOZYME_C::ijall[j];
        if (ik < ind_j) {
            il = j;
            j = (j + iu + 1) / 2;
        } else {
            if (ik == ind_j) break;
            iu = j;
            j = (j + il) / 2;
            if (j == jm2) return -2;
            jm1 = j;
            jm2 = jm1;
        }
    }
    return MOZYME_C::iijj[j];
}
