// convert_storage.cpp — C++ translation of MOPAC 2016 "convert_storage.F90".
// ijbo / memory_error / mopend are external stubs.

#include "convert_storage.h"

#include <algorithm>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace MOZYME_C;

namespace {
int ijbo(int, int) { return -1; }
void memory_error(const char*) {}
void mopend(const char*) {}
}

void convert_mat_packed_to_triangle(const std::vector<double>& matrix_packed,
                                    std::vector<double>& matrix_triangle) {
    int linear = norbs * (norbs + 1) / 2;
    matrix_triangle.assign(linear + 1, 0.0);
    for (int i = 1; i <= numat; ++i) {
        for (int j = 1; j <= i; ++j) {
            if (ijbo(i, j) >= 0) {
                int ij = ijbo(i, j);
                int il = nfirst[i], iu = nlast[i];
                int jl = nfirst[j], ju = nlast[j];
                for (int ii = il; ii <= iu; ++ii)
                    for (int jj = jl; jj <= std::min(ju, ii); ++jj) {
                        ij++;
                        matrix_triangle[ii * (ii - 1) / 2 + jj] = matrix_packed[ij];
                    }
            }
        }
    }
}

void convert_lmo_packed_to_square(std::vector<std::vector<double>>& c_square) {
    int nocc = nelecs / 2;
    int nvir = norbs - nocc;
    if (isort.empty()) {
        isort.assign(norbs + 1, 0);
        for (int i = 1; i <= nocc; ++i) isort[i] = i;
        for (int i = 1; i <= nvir; ++i) isort[nocc + i] = i;
    }
    for (int iunsrt = 1; iunsrt <= nocc; ++iunsrt) {
        int i = isort[iunsrt];
        for (int j = 1; j <= norbs; ++j) c_square[j][iunsrt] = 0.0;
        int ka = ncocc[i];
        for (int jj = nncf[i] + 1; jj <= nncf[i] + ncf[i]; ++jj) {
            int j = icocc[jj];
            for (int k = nfirst[j]; k <= nlast[j]; ++k) {
                ka++;
                c_square[k][iunsrt] = cocc[ka];
            }
        }
    }
    for (int iunsrt = 1; iunsrt <= nvir; ++iunsrt) {
        int i = isort[iunsrt + nocc];
        int ii = iunsrt + nocc;
        for (int j = 1; j <= norbs; ++j) c_square[j][ii] = 0.0;
        int ka = ncvir[i];
        for (int jj = nnce[i] + 1; jj <= nnce[i] + nce[i]; ++jj) {
            int j = icvir[jj];
            for (int k = nfirst[j]; k <= nlast[j]; ++k) {
                ka++;
                c_square[k][ii] = cvir[ka];
            }
        }
    }
}
