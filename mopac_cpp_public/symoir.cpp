// symoir.cpp — C++ translation of MOPAC 2016 "symoir.F90".
// Assign irreducible representations to eigenvectors by matching characters.
// charmo is fully translated; charmvi/charmst have no 2016 source (permanent
// gap) and are supplied by the caller as stubs.
#include "symoir.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "charmo.h"
#include "chanel_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"

using namespace meci_C;
using namespace molkst_C;
using namespace symmetry_C;

// charmvi/charmst: no 2016 MOPAC source exists (legacy calls kept for shape);
// permanent placeholder returning 0.0 until M05 batch confirms semantics.
double charmvi(const std::vector<std::vector<double>>&, int, int,
               const std::vector<std::vector<double>>&, int) { return 0.0; }
double charmst(const std::vector<std::vector<double>>&, const std::vector<int>&,
               int, int, const std::vector<std::vector<double>>&, int, bool&) { return 0.0; }

extern double charmvi(const std::vector<std::vector<double>>& vects, int i, int j,
                      const std::vector<std::vector<double>>& r, int nvecs);
extern double charmst(const std::vector<std::vector<double>>& vects,
                      const std::vector<int>& ntype, int i, int j,
                      const std::vector<std::vector<double>>& r, int nvecs, bool& first);

void symoir(int itype, double* vects, double* eigs, int nvecs, double* r, int imat) {
    (void)imat;
    const double toler = 0.1;
    const bool r3 = (name == "R3");
    int korb = 0;
    std::vector<int> ntype(nvecs + 1, 0);
    for (int i = 1; i <= numat; ++i) {
        const int jj = jndex[i];
        for (int j = 1; j <= jj; ++j) { ++korb; ntype[korb] = 100 * i + 9 + j; }
    }
    const int nfind = (itype != 3) ? nvecs : lab;
    std::vector<int> icount(nirred + 1, 0);
    std::string names = "????";
    if (nclass == 1) names = jx[1];
    for (int i = 1; i <= nfind; ++i) { jndex[i] = i; namo[i] = names; }
    if (nclass == 1) return;

    // Rebuild vects as (nvecs+1) x (nvecs+1) direct-index rows for charmo.
    std::vector<std::vector<double>> vv(nvecs + 1, std::vector<double>(nvecs + 1, 0.0));
    for (int row = 1; row <= nvecs; ++row)
        for (int col = 1; col <= nvecs; ++col)
            vv[row][col] = vects[(col - 1) * nvecs + (row - 1)];  // F90 column-major
    std::vector<std::vector<double>> rr(4, std::vector<double>(4, 0.0));
    for (int a = 1; a <= 3; ++a)
        for (int b = 1; b <= 3; ++b) rr[a][b] = r[(b - 1) * 3 + (a - 1)];

    std::vector<std::vector<double>> carmat(nfind + 1, std::vector<double>(nclass + 1, 0.0));
    bool first = true;
    for (int j = 1; j <= nclass; ++j) {
        for (int i = 1; i <= nfind; ++i) {
            if (itype == 1) carmat[i][j] = charmo(vv, ntype, i, j, rr, nvecs, first);
            else if (itype == 2) carmat[i][j] = charmvi(vv, i, j, rr, nvecs);
            else carmat[i][j] = charmst(vv, ntype, i, j, rr, nvecs, first);
        }
    }
    if (itype == 3) (void)charmst(vv, ntype, -1, 0, rr, nvecs, first);

    const bool debug = (keywrd.find("SYMOIR") != std::string::npos);
    if (debug) {
        std::printf(" Characters of Transform\n");
        for (int i = 1; i <= nfind; ++i) {
            std::printf("%5d", i);
            for (int j = 1; j <= nclass; ++j) std::printf("%12.6f", carmat[i][j]);
            std::printf("\n");
        }
    }
    int i = 0;
l70:
    {
        const int ik = i + 1;
        std::vector<double> tchar(nclass + 1, 0.0);
l90:
        ++i;
        if (i > nfind) {
            if (itype == 1) {
                for (int ii = 1; ii <= nfind; ++ii) {
                    char c0 = namo[ii][0];
                    if (c0 >= 'A' && c0 <= 'Z') namo[ii][0] = (char)(c0 - 'A' + 'a');
                }
            }
            if (debug) {
                std::printf("\n Number of Irreducible Representations of each Class\n");
                for (int q = 1; q <= nirred; ++q)
                    std::printf("%5d %s", icount[q], jx[q].c_str());
            }
            return;
        }
        for (int j = 1; j <= nclass; ++j) tchar[j] += carmat[i][j];
        if (tchar[1] > 5.1 && !r3) goto l70;
        for (int k = 1; k <= nirred; ++k) {
            bool ok = true;
            for (int j = 1; j <= nclass; ++j) {
                if (std::fabs(tchar[j] - group[(j - 1) * 20 + (k - 1)]) > toler) { ok = false; break; }
            }
            if (!ok) continue;
            ++icount[k];
            if (i - ik + 1 > 0) {
                for (int q = ik; q <= i; ++q) { jndex[q] = icount[k]; namo[q] = jx[k]; }
            }
            goto l70;
        }
        if (i < nfind) if (eigs[i] - eigs[i - 1] > 0.1) goto l70;  // F90 eigs(i+1)-eigs(i), 0-based
        goto l90;
    }
}
