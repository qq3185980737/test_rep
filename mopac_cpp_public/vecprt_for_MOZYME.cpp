// vecprt_for_MOZYME.cpp — C++ translation of "vecprt_for_MOZYME.F90".
// MOZYME variant of vecprt: remaps the packed interaction array aa through
// ijbo before printing, or prints selected atoms when l_atom is set.

#include "vecprt_for_MOZYME.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "ijbo.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

static const char* atorbs[9] = {" S", "PX", "PY", "PZ", "X2",
                                "XZ", "Z2", "YZ", "XY"};
static constexpr int maxarr = 200;

void vecprt_for_MOZYME(double* aa, int numm) {
    int numb = std::abs(numm);
    if (numb > maxarr) {
        std::printf("VECPRT CAN ONLY PRINT ARRAYS OF SIZE .LE. %d\n", maxarr);
        std::printf("AN ATTEMPT WAS MADE TO PRINT AN ARRAY OF SIZE %d\n", numb);
        numb = maxarr;
    }
    int linear = numb * (numb + 1) / 2;
    std::vector<double> a(linear + 1, 0.0);

    std::vector<std::string> itext(maxarr + 1, "  ");
    std::vector<std::string> jtext(maxarr + 1, "  ");
    std::vector<int> natom(maxarr + 1, 0);

    if (numat != 0 && numat == numb) {
        // Option (1): over selected atoms via l_atom.
        int j = 0, l = 0;
        for (int i = 1; i <= numb; ++i) {
            if (l_atom[i]) {
                ++j;
                itext[j] = "  ";
                jtext[j] = elemnt[nat[i]];
                natom[j] = gui ? j : i;
                for (int k = 1; k <= i; ++k)
                    if (l_atom[k]) {
                        ++l;
                        a[l] = aa[(i * (i - 1)) / 2 + k - 1];
                    }
            }
        }
        numb = j;
        linear = l;
    } else if (numat != 0 && nlast[numat] == numm) {
        // Option (2): over atomic orbitals; remap through ijbo.
        for (int i = 1; i <= linear; ++i) a[i] = 0.0;
        int ij = 0;
        int ll = 0;
        for (int i = 1; i <= numat; ++i) {
            for (int j = 1; j <= i - 1; ++j) {
                if (ijbo(i, j) >= 0) {
                    ++ij;
                    int l = ijbo(i, j);
                    for (int ii = nfirst[i]; ii <= nlast[i]; ++ii)
                        for (int jj = nfirst[j]; jj <= nlast[j]; ++jj) {
                            ++l;
                            ll = (ii * (ii - 1)) / 2 + jj;
                            if (ll <= linear) a[ll] = aa[l - 1];
                        }
                }
            }
            ++ij;
            int l = ijbo(i, i);
            for (int ii = nfirst[i]; ii <= nlast[i]; ++ii)
                for (int jj = nfirst[i]; jj <= ii; ++jj) {
                    ++l;
                    ll = (ii * (ii - 1)) / 2 + jj;
                    if (ll <= linear) a[ll] = aa[l - 1];
                }
            if (ll == linear) break;
        }
        bool outer_done = false;
        for (int i = 1; i <= numat && !outer_done; ++i) {
            int jlo = nfirst[i], jhi = nlast[i];
            int l = nat[i];
            int k = 0;
            for (int j = jlo; j <= jhi; ++j) {
                ++k;
                itext[j] = atorbs[k - 1];
                jtext[j] = elemnt[l];
                natom[j] = i;
                if (j == numb) { outer_done = true; break; }
            }
        }
    } else {
        // Option (3): generic array.
        for (int i = 1; i <= numb; ++i) {
            itext[i] = "  ";
            jtext[i] = "  ";
            natom[i] = i;
        }
        for (int i = 1; i <= linear; ++i) a[i] = aa[i - 1];
    }

    // Scale diagonal terms.
    double sumax = 1.0;
    for (int i = 1; i <= numb; ++i)
        sumax = std::max(sumax, std::fabs(a[i * (i + 1) / 2]));
    int ie = static_cast<int>(std::log10(sumax));
    if (ie == 1 || ie == 2) ie = 0;
    double fact = std::pow(10.0, -ie);
    if (std::fabs(fact - 1.0) > 0.001) {
        std::printf("Diagonal Terms should be Multiplied by%16.6f\n", 1.0 / fact);
        for (int i = 1; i <= numb; ++i) a[i * (i + 1) / 2] *= fact;
    }

    std::string sep(126, '-');
    int limit = numb * (numb + 1) / 2;
    int kk = 8, na = 1;
    while (true) {
        int ll = 0;
        int m = std::min(numb + 1 - na, 6);
        int ma = 2 * m + 1;
        int mend = na + m - 1;
        std::printf("\n\n%13s", "");
        for (int i = na; i <= mend; ++i)
            std::printf(" %2s %2s%3d  ", itext[i].c_str(), jtext[i].c_str(),
                        natom[i]);
        std::printf("\n ");
        for (int i = 1; i <= ma; ++i) std::printf("%-6s", "------");
        std::printf("\n");
        bool done = false;
        for (int i = na; i <= numb; ++i) {
            ++ll;
            int k = i * (i - 1) / 2;
            int lend = std::min(k + (mend - na + 1), k + i);
            k += na;
            if (kk + ll > 50) {
                std::printf("\n\n%13s", "");
                for (int n = na; n <= mend; ++n)
                    std::printf(" %2s %2s%3d  ", itext[n].c_str(),
                                jtext[n].c_str(), natom[n]);
                std::printf("\n ");
                for (int n = 1; n <= ma; ++n) std::printf("%-6s", "------");
                std::printf("\n");
                kk = 4;
                ll = 0;
            }
            std::printf(" %2s %2s%5d", itext[i].c_str(), jtext[i].c_str(),
                        natom[i]);
            for (int n = k; n <= lend; ++n) std::printf("%11.6f", a[n]);
            std::printf("\n");
            if (lend >= limit) { done = true; break; }
        }
        if (done) break;
        kk = kk + ll + 4;
        na = mend + 1;
        if (kk + numb + 1 - na > 50) kk = 4;
    }
}
