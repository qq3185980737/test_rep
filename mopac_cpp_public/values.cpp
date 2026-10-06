// values.cpp — C++ translation of MOPAC 2016 "values.F90".
// Calculates the energy levels of the localised MOs and sorts the LMOs into
// increasing energy order (valuen helper included).
#include "values.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include <cmath>
#include <string>
#include <vector>
using namespace MOZYME_C;

namespace {
using molkst_C::norbs;
using molkst_C::mpack;
using chanel_C::iw;
using common_arrays_C::f;
using common_arrays_C::eigs;
}  // namespace

// ijbo: packed index for atom-atom block (external, Fortran ijbo semantics).
extern int ijbo(int n1, int n2);

// valuen(fao, nfao, nocc, nncf, ncf, nnn, icocc, nico, ncocc, iorbs,
//        cocc, nco, fdiat, eigf, isort)
static void valuen(const double* fao, int nfao, int nocc, const std::vector<int>& nncf,
                   const std::vector<int>& ncf, int nnn, const std::vector<int>& icocc,
                   int nico, const std::vector<int>& ncocc, const std::vector<int>& iorbs,
                   const std::vector<double>& cocc, int nco, std::vector<double>& fdiat,
                   std::vector<double>& eigf, std::vector<int>& isort) {
    for (int i = 1; i <= nocc; ++i) {
        int loopi = ncocc[i];
        int l = 0;
        double sum;
        for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) {
            int j1 = icocc[j];
            sum = 0.0;
            for (int k = l + 1; k <= l + iorbs[j1]; ++k)
                sum = sum + cocc[k + loopi] * cocc[k + loopi];
            l = l + iorbs[j1];
        }
            sum = 0.0;
            int jl = loopi;
        for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) {
            int j1 = icocc[j];
            for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                ++jl;
                int kl = loopi;
                for (int k = nncf[i] + 1; k <= nncf[i] + ncf[i]; ++k) {
                    int k1 = icocc[k];
                    if (ijbo(k1, j1) >= 0) {
                        int ii;
                        if (k1 > j1) {
                            // LOWER TRIANGLE
                            ii = ijbo(k1, j1) + jx - iorbs[j1];
                            for (int i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                ii = ii + iorbs[j1];
                                fdiat[i4] = fao[ii];
                            }
                        } else if (k1 < j1) {
                            // UPPER TRIANGLE
                            ii = ijbo(k1, j1) + iorbs[k1] * (jx - 1);
                            for (int i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                ++ii;
                                fdiat[i4] = fao[ii];
                            }
                        } else {
                            // DIAGONAL TERM
                            ii = ijbo(k1, j1) + (jx * (jx + 1)) / 2;
                            for (int i4 = jx + 1; i4 <= iorbs[k1]; ++i4) {
                                ii = ii + jx;
                                fdiat[i4] = fao[ii];
                                ii = ii + i4 - jx;
                            }
                            ii = ijbo(k1, j1) + (jx * (jx - 1)) / 2;
                            for (int j4 = 1; j4 <= jx; ++j4) {
                                ++ii;
                                fdiat[j4] = fao[ii];
                            }
                        }
                        double sum1 = 0.0;
                        for (int i4 = 1; i4 <= iorbs[k1]; ++i4) {
                            ++kl;
                            sum1 = sum1 + fdiat[i4] * cocc[kl];
                        }
                        sum = sum + cocc[jl] * sum1;
                    }
                }
            }
        }
        eigf[i] = sum;
    }
    // Sort eigenvalues.
    for (int i = 1; i <= nocc; ++i) fdiat[i] = eigf[i];
    for (int i = 1; i <= nocc; ++i) {
        double sum = 1.0e9;
        int k = 0;
        for (int j = 1; j <= nocc; ++j) {
            if (fdiat[j] < sum) {
                k = j;
                sum = fdiat[j];
            }
        }
        fdiat[k] = 1.0e10;
        isort[i] = k;
    }
}

void values(const std::string& type) {
    // MOZYME module data (global linkage, declared in MOZYME_C.h).
    using MOZYME_C::nvirtual;
using MOZYME_C::noccupied;
    using MOZYME_C::iorbs;

    if ((int)isort.size() < norbs + 1) isort.resize(norbs + 1, 0);
    std::vector<double> fdiat(norbs + 1, 0.0);
    if (type == "OCCUPIED") {
        valuen(f.data(), mpack, noccupied, nncf, ncf, norbs, icocc, icocc_dim, ncocc,
               iorbs, cocc, cocc_dim, fdiat, eigs, isort);
    } else if (type == "VIRTUAL") {
        std::vector<double> eigf(nvirtual + 1, 0.0);
        std::vector<int> isort_v(nvirtual + 1, 0);
        for (int i = 1; i <= nvirtual; ++i) {
            eigf[i] = eigs[noccupied + i];
            isort_v[i] = 0;
        }
        valuen(f.data(), mpack, nvirtual, nnce, nce, norbs, icvir, icvir_dim, ncvir,
               iorbs, cvir, cvir_dim, fdiat, eigf, isort_v);
        for (int i = 1; i <= nvirtual; ++i) {
            eigs[noccupied + i] = eigf[i];
            isort[noccupied + i] = isort_v[i];
        }
    } else {
        std::fprintf(stderr, " Error\n");
        extern void mopend(const std::string&);
        mopend("Error");
        return;
    }
}
