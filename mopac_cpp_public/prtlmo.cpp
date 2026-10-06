// prtlmo.cpp
#include "prtlmo.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
using namespace MOZYME_C;
using common_arrays_C::nat;
using common_arrays_C::eigs;
using molkst_C::norbs;
using MOZYME_C::iorbs;
using MOZYME_C::noccupied;
using MOZYME_C::nvirtual;
extern const char* elemnt(int);
extern int iw;

static void prtlmn(const std::vector<int>& nncx, const std::vector<int>& icxxx,
                   const std::vector<int>& ncxxx, const std::vector<double>& cxxx,
                   const std::vector<int>& ncx, const std::vector<double>& eigv,
                   int mmos, int i_offset) {
    std::vector<double> eigs_temp(mmos + 1);
    for (int i = 1; i <= mmos; ++i) eigs_temp[i] = eigv[i];
    for (int i = 1; i <= mmos; ++i) {
        double eig_min = 1e7; int ii = 1;
        for (int j = 1; j <= mmos; ++j)
            if (eigs_temp[j] < eig_min) { eig_min = eigs_temp[j]; ii = j; }
        isort[i] = ii;
        eigs_temp[ii] = 1e8;
    }
    int cap = std::max(1, noccupied * 20);
    std::vector<double> w(cap + 1);
    std::vector<int> iscrch(cap + 1);
    std::vector<int> jat(norbs + 1, 0);
    double cnst = 2e-4;
    int j = 0;
    for (int iunsrt = 1; iunsrt <= mmos; ++iunsrt) {
        int i = isort[iunsrt];
        int n = ncxxx[i];
        int k = 0;
        for (int ii = nncx[i] + 1; ii <= nncx[i] + ncx[i]; ++ii) {
            int l = icxxx[ii];
            int nj = iorbs[l];
            double sum = 0;
            for (int m = 1; m <= nj; ++m) sum += cxxx[m + n] * cxxx[m + n];
            if (sum > cnst) {
                ++j; ++k;
                iscrch[j] = l;
                w[j] = sum * 2.0;
            }
            n += nj;
        }
        jat[i] = k;
    }
    int ju = 0;
    std::vector<double> xbig(101, 0);
    std::vector<int> ibig(101, 0), jbig(101, 0);
    for (int iunsrt = 1; iunsrt <= mmos; ++iunsrt) {
        int i = isort[iunsrt];
        int jl = ju + 1;
        ju = jl + jat[i] - 1;
        int k = 0;
        for (int jj = jl; jj <= ju; ++jj) {
            ++k;
            xbig[k] = w[jj];
            ibig[k] = iscrch[jj];
        }
        int kk = std::min(20, k);
        int loop;
        for (loop = 1; loop <= kk; ++loop) {
            double sum = -1; int l = 0, lj = 0;
            for (int jj = loop; jj <= k; ++jj)
                if (xbig[jj] > sum) { sum = xbig[jj]; l = ibig[jj]; lj = jj; }
            if (sum < cnst) break;
            xbig[lj] = xbig[loop];
            ibig[lj] = ibig[loop];
            jbig[loop] = (int)(sum / cnst);
            ibig[loop] = l;
        }
        std::printf("  LMO %d E=%.4f", iunsrt + i_offset, eigv[i]);
        for (int jj = 1; jj < loop; ++jj)
            std::printf("  %s:%d", elemnt(nat[ibig[jj]]), jbig[jj]);
        std::printf("\n");
    }
}

void prtlmo() {
    if (isort.empty()) isort.assign(norbs + 1, 0);
    std::printf(" LOCALIZED MOLECULAR ORBITALS\n");
    prtlmn(nncf, icocc, ncocc, cocc, ncf, eigs, noccupied, 0);
    std::vector<double> ev(nvirtual + 1, 0);
    for (int i = 1; i <= nvirtual; ++i) ev[i] = eigs[noccupied + i];
    prtlmn(nnce, icvir, ncvir, cvir, nce, ev, nvirtual, noccupied);
}
