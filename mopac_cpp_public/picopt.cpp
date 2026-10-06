// picopt.cpp
#include "picopt.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
using namespace MOZYME_C;
using common_arrays_C::loc;
using common_arrays_C::nbonds;
using common_arrays_C::ibonds;
using common_arrays_C::labels;
using molkst_C::natoms;
using molkst_C::numat;
using molkst_C::nvar;
using molkst_C::ndep;
using molkst_C::numcal;
using symmetry_C::locdep;

static int imol = 0;

void picopt(int loop) {
    std::vector<int> iopt_loc(natoms + 1, 0);
    if (loop == -1) {
        for (int i = 1; i <= natoms; ++i) iopt_loc[i] = 1;
    } else {
        for (int i = 1; i <= nvar; ++i) iopt_loc[loc[1][i]] = 2;
        for (int i = 1; i <= ndep; ++i) iopt_loc[locdep[i]] = 2;
        int j = 0;
        for (int i = 1; i <= natoms; ++i) {
            if (labels[i] != 99) {
                ++j;
                iopt_loc[j] = iopt_loc[i];
            }
        }
        if (imol == numcal) {
            for (int i = 1; i <= numat; ++i) {
                if (iopt_loc[i] == 2) {
                    for (int jj = 1; jj <= nbonds[i]; ++jj) {
                        if (iopt_loc[ibonds[jj][i]] == 0)
                            iopt_loc[ibonds[jj][i]] = 1;
                    }
                }
            }
        }
    }
    imol = numcal;
    numred = 0;
    for (int i = 1; i <= numat; ++i) {
        if (iopt_loc[i] != 0) {
            ++numred;
            jopt[numred] = i;
        }
    }
}
