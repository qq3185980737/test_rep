// vecprt.cpp — C++ translation of "vecprt.F90".
// Prints a lower-half-triangle of a square matrix stored in packed form.
// If MOZYME mode, delegates to vecprt_for_MOZYME. Output unit iw -> stdout.

#include "vecprt.h"
#include "vecprt_for_MOZYME.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

static const char* atorbs[9] = {" S", "PX", "PY", "PZ", "X2",
                                "XZ", "Z2", "YZ", "XY"};

void vecprt(std::vector<double> const& a, int numm) { vecprt(const_cast<double*>(a.data()), numm); }
void vecprt(double* a, int numm) {
    if (mozyme) {
        vecprt_for_MOZYME(a, numm);
        return;
    }
    const int numb = std::abs(numm);
    // Fortran a(*) is 1-based packed; copy into a padded local vector.
    const int limit = numb * (numb + 1) / 2;
    std::vector<double> av(limit + 1, 0.0);
    for (int i = 1; i <= limit; ++i) av[i] = a[i - 1];

    std::vector<int> natom(numb + 1, 0);
    std::vector<std::string> itext(numb + 1, "  ");
    std::vector<std::string> jtext(numb + 1, "  ");

    double sumax = 1.0;
    for (int i = 1; i <= numm; ++i)
        sumax = std::max(sumax, std::fabs(av[i * (i + 1) / 2]));
    int ie = static_cast<int>(std::log10(sumax));
    if (ie == 1 || ie == 2) ie = 0;
    double fact = std::pow(10.0, -ie);
    if (std::fabs(fact - 1.0) > 0.001) {
        std::printf("Diagonal Terms should be Multiplied by%12.1f\n", 1.0 / fact);
        for (int i = 1; i <= numm; ++i) av[i * (i + 1) / 2] *= fact;
    }

    if (numat != 0 && numat == numm) {
        // Option (1): over atoms.
        for (int i = 1; i <= numat; ++i) {
            itext[i] = "  ";
            jtext[i] = elemnt[nat[i]];
            natom[i] = i;
        }
    } else if (numat != 0 && nlast[numat] == numm) {
        // Option (2): over atomic orbitals.
        for (int i = 1; i <= numat; ++i) {
            int jlo = nfirst[i];
            int jhi = nlast[i];
            int l = nat[i];
            for (int k = 0; k <= jhi - jlo; ++k) {
                itext[jlo + k] = atorbs[k];
                jtext[jlo + k] = elemnt[l];
                natom[jlo + k] = i;
            }
        }
    } else {
        for (int i = 1; i <= numb; ++i) {
            itext[i] = "  ";
            jtext[i] = "  ";
            natom[i] = i;
        }
    }

    std::string sep(126, '-');  // 21 * 6 dashes
    int kk = 8, na = 1;
    while (na <= numb) {
        int ll = 0;
        int m = std::min(numb + 1 - na, 6);
        int ma = 2 * m + 1;
        int mend = na + m - 1;
        // header row: 13 spaces, then m groups of " a2 a2 i3  "
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
            int lend = std::min(k + m, k + i);
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
            for (int n = k; n <= lend; ++n) std::printf("%11.6f", av[n]);
            std::printf("\n");
            if (lend >= limit) { done = true; break; }
        }
        if (done) break;
        kk = kk + ll + 4;
        na = mend + 1;
        if (kk + numb + 1 - na <= 50) {
            // continue same page
        } else {
            kk = 4;
        }
    }

    if (std::fabs(fact - 1.0) > 0.001)
        for (int i = 1; i <= numm; ++i) av[i * (i + 1) / 2] /= fact;
}
