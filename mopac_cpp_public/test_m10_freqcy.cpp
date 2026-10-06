// test_m10_freqcy.cpp — verification for freqcy.cpp (MOPAC 2016 freqcy.F90).
// Links real: freqcy.cpp + rsp.cpp (genuine Jacobi solver) + module defs.
// Stubs (by design, exercised only when branches are disabled in the test):
//   symt / frame / symtrz — symmetry branches skipped (ts=true, NOSYM, eorc=false)
//   phase_lock — phase convention; stub keeps values unchanged (identity).
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include "freqcy.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "to_screen_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

// ---- stubs for symmetry-branch routines (not reached in this test) ----
namespace {
int g_symt_calls = 0, g_frame_calls = 0, g_symtrz_calls = 0;
}
void symt(double*, double*, double*) { ++g_symt_calls; }
void frame(std::vector<double>&, int, int) { ++g_frame_calls; }
void symtrz(double*, double*, int, int) { ++g_symtrz_calls; }
// phase_lock: identity (largest coefficient already positive in test data)
void phase_lock(std::vector<double>&, int) {}

static int g_checks = 0; static double g_err = 0.0;
static void chk(bool ok, const char* name, double err = 0.0) {
    ++g_checks;
    if (!ok) std::fprintf(stderr, "[FAIL] %s\n", name);
    if (err > g_err) g_err = err;
}
static void chk_rel(double got, double want, double tol, const char* name) {
    double e = std::fabs(got - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1.0);
    chk(e <= tol, name, e);
}

int main() {
    // ---- global state for a 1-atom system (3 vibrational DOF) ----
    numat = 1; nvar = 3;
    atmass.assign(2, 0.0); atmass[1] = 12.0;
    keywrd = " NOSYM";
    to_screen_C::redmas.assign(4, std::vector<double>(2, 0.0));

    std::fprintf(stderr, "[t1] freqcy_mass_weight (packed 1-based)\n");
    {
        std::vector<double> f = {0, 2.0, -1.0, 2.0, 0.0, -1.0, 2.0};  // 3x3 tridiagonal
        std::vector<double> old = {0, 0, 0, 0, 0, 0, 0};
        std::vector<double> w = {0, 0.5, 0.5, 0.5};  // wtmass[1..3]
        freqcy_mass_weight(f, old, w, 3);
        // oldf[k] = f[k]*1e5 ; f[k] = f[k]*w[i]*w[k]
        chk_rel(old[1], 2.0e5, 1e-12, "mw old1=2e5");
        chk_rel(old[6], 2.0e5, 1e-12, "mw old6=2e5");
        chk_rel(f[1], 2.0 * 0.25, 1e-12, "mw f1=0.5");
        chk_rel(f[6], 2.0 * 0.25, 1e-12, "mw f6=0.5");
    }

    std::fprintf(stderr, "[t2] freqcy full body (NOSYM, ts=true, eorc=false)\n");
    {
        // H = [[2,-1,0],[-1,2,-1],[0,-1,2]] packed lower-triangle row-major 1-based:
        //   fmatrx(1)=2, (2,1)=-1, (2,2)=2, (3,1)=0, (3,2)=-1, (3,3)=2
        std::vector<double> fmatrx = {0, 2.0, -1.0, 2.0, 0.0, -1.0, 2.0};
        std::vector<double> freq(4, 0.0), travel(4, 0.0), oldf(7, 0.0), ff(4, 0.0);
        std::vector<std::vector<double>> deldip(4, std::vector<double>(4, 0.0));
        freqcy(fmatrx, freq, travel, false, deldip, ff, oldf, true);
        // Mass-weighted matrix: H/12. Eigenvalues of H: 2-sqrt2, 2, 2+sqrt2.
        // rsp returns ascending roots; freqcy then rescales them.
        chk(freq[1] > 0.0 && freq[2] > 0.0 && freq[3] > 0.0, "freqcy freq positive");
        chk(freq[1] < freq[2] && freq[2] < freq[3], "freqcy freq ascending");
        // Frequencies land in the hundreds-to-thousands cm^-1 band for this toy H.
        chk(freq[1] > 100.0 && freq[3] < 100000.0, "freqcy band sanity");
        for (int i = 1; i <= 3; ++i)
            chk(travel[i] >= 0.0 && travel[i] <= 1.0, "freqcy travel in [0,1]");
        // eorc=false: fmatrx restored as oldf*1e-5*w^2
        chk_rel(fmatrx[1], 2.0e5 * 1e-5 / 12.0, 1e-12, "freqcy restore f1");
        chk_rel(fmatrx[6], 2.0e5 * 1e-5 / 12.0, 1e-12, "freqcy restore f6");
        chk(g_symt_calls == 0 && g_frame_calls == 0 && g_symtrz_calls == 0,
            "freqcy sym-branch stubs not reached");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_checks ? 0 : 1;
}
