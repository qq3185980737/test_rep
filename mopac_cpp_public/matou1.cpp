// matou1.cpp — C++ translation of "matou1.F90".
// Prints a square matrix of eigenvectors and eigenvalues (row labels
// depend on iflag).  F90 1-based indexing maps directly; the F90
// to_screen_C copy block (allocate->copy->deallocate, a no-op) is dropped.
#include "matou1.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "molkst_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "symmetry_C.h"

using namespace molkst_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace symmetry_C;

namespace {

// trim trailing blanks (Fortran trim semantics for H_txt slices)
std::string ftrim(const std::string& s) {
    size_t e = s.size();
    while (e > 0 && s[e - 1] == ' ') --e;
    return s.substr(0, e);
}

}  // namespace

void matou1(double* a, double* b, int ncx, int& nr, int ndim, int iflag) {
    (void)ndim;  // F90 uses it only for the dropped to_a/to_b copy block
    std::vector<int> natom(nr + 1, 0);
    std::vector<std::string> itext(nr + 1), jtext(nr + 1);
    static const char* xyz[3] = {" x", " y", " z"};
    static const char* atorbs[9] = {"S ", "Px", "Py", "Pz", "x2", "xz",
                                    "z2", "yz", "xy"};
    bool allprt = keywrd.find("ALLVEC") != std::string::npos;
    int i = (int)keywrd.find(" VECTORS") + 1;   // 1-based position
    int nup = 7, ndown = 9;
    if (i > 0) {
        int j = i + 1;
        for (; j <= i + 20 && j <= (int)keywrd.size(); ++j)
            if (keywrd[j - 1] == ' ') break;
        if (keywrd.substr(i - 1, j - i + 1).find(')') != std::string::npos) {
            //  User-specified counts:  VECTORS=(12,34)
            int k = i + 1;
            for (; k <= j; ++k)
                if (keywrd[k - 1] >= '0' && keywrd[k - 1] <= '9') break;
            ndown = keywrd[k - 1] - '0';
            for (++k; k <= j; ++k) {
                if (keywrd[k - 1] < '0' || keywrd[k - 1] > '9') break;
                ndown = ndown * 10 + keywrd[k - 1] - '0';
            }
            for (; k <= j; ++k)
                if (keywrd[k - 1] >= '0' && keywrd[k - 1] <= '9') break;
            nup = keywrd[k - 1] - '0';
            for (++k; k <= j; ++k) {
                if (keywrd[k - 1] < '0' || keywrd[k - 1] > '9') break;
                nup = nup * 10 + keywrd[k - 1] - '0';
            }
        }
        ndown = ndown - 1;
        if (nup == -16) {   // no digits before ',': split remaining half/half
            nup = (ndown + 1) / 2;
            ndown = ndown - nup;
        }
    }
    int nc = ncx;
    int nfix = 0;
    if (iflag > 2 && iflag != 5) goto label_50;
    if (!allprt) {
        int nsave = ncx;
        nfix = std::max(nalpha, nclose);
        if (iflag == 2 && nc > 16) nc = nfix + nup;
        nc = std::min(nsave, nc);
    }
    if (numat == 0) goto label_50;
    if (nlast[numat] != nr) goto label_50;
    for (i = 1; i <= numat; ++i) {
        int jlo = nfirst[i], jhi = nlast[i];
        int l = nat[i];
        if (iflag <= 2) {
            for (int t = jlo; t <= jhi; ++t) itext[t] = atorbs[t - jlo];
            for (int t = jlo; t <= jhi; ++t) jtext[t] = elemnt[l];
            for (int t = jlo; t <= jhi; ++t) natom[t] = i;
        } else {
            jhi = 3 * (i - 1);
            for (int t = 1; t <= 3; ++t) {
                itext[t + jhi] = xyz[t - 1];
                jtext[t + jhi] = elemnt[l];
                natom[t + jhi] = i;
            }
        }
    }
    goto label_70;
label_50:
    nr = std::abs(nr);
    if (iflag == 3) {
        for (i = 1; i <= nr; ++i) {
            itext[i] = "  ";
            jtext[i] = elemnt[nat[i]];
            natom[i] = i;
        }
    } else {
        for (i = 1; i <= nr; ++i) {
            itext[i] = "  ";
            jtext[i] = "  ";
            natom[i] = i;
        }
    }
label_70:
    int ka = 1, kc = 8;
    if (!allprt) {
        if (iflag == 2 && norbs > 16) ka = nfix - ndown;
        ka = std::max(1, ka);
        if (iflag == 2 && norbs > 16) kc = ka + 7;
    }
label_90:
    int kb = std::min(kc, nc);
    // 130 format: (/,/,2x,' Root No.',i8,11i10)
    std::printf("\n\n  Root No.%8d", ka);
    for (int t = ka + 1; t <= kb; ++t) std::printf("%10d", t);
    if (iflag == 2 || iflag == 5) {
        // 180 format: (/,13x,10(i5,1x,a4))
        std::printf("\n             ");
        for (int t = ka; t <= kb; ++t) std::printf("%5d %4s", jndex[t], namo[t].c_str());
    }
    if (b[1] != 0.0) {
        if (iflag == 5) {
            std::printf("\n            ");
            for (int t = ka; t <= kb; ++t) std::printf("%10.1f", b[t]);
        } else {
            std::printf("\n            ");
            for (int t = ka; t <= kb; ++t) std::printf("%10.3f", b[t]);
        }
    }
    std::printf("  \n");   // 160 format: '  '
    int la = 1, lc = 40;
label_100:
    int lb = std::min(lc, nr);
    for (i = la; i <= lb; ++i) {
        if (itext[i] == " S") std::printf("  \n");
        // 170 format: (' ',2(1x,a2),i5,f10.4,10f10.4)
        std::printf(" %2s %2s%5d%10.4f", itext[i].c_str(), jtext[i].c_str(),
                    natom[i], a[(ka - 1) * nr + (i - 1)]);
        for (int t = ka + 1; t <= kb; ++t)
            std::printf("%10.4f", a[(t - 1) * nr + (i - 1)]);
        std::printf("\n");
    }
    if (lb == nr) goto label_120;
    la = lc + 1;
    lc = lc + 40;
    goto label_100;
label_120:
    if (kb == nc) return;
    ka = kc + 1;
    kc = kc + 8;
    goto label_90;
}
