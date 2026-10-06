// test_m09_gdisp.cpp — verification for exchng.cpp + gdisp.cpp (MOPAC 2016).
// exchng: pure value-exchange routine.  gdisp: D3 dispersion gradient driver
// (ncoord / get_dC6_dCNij / gdisp), verified on a 2-atom toy system.
#include <cstdio>
#include <cmath>
#include <vector>
#include "exchng.h"
#include "gdisp.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

static int g_checks = 0; static int g_fail = 0; static double g_err = 0.0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    double e = std::fabs(got - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1.0);
    if (e > g_err) g_err = e;
    if (e > tol) { ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.12g want %.12g\n", name, got, want); }
}

int main() {
    std::fprintf(stderr, "[t1] exchng (value exchange)\n");
    {
        double a = 1.5, c = 2.5, t = 3.5;
        double b = 0.0, d = 0.0, q = 0.0;
        std::vector<double> x = {0, 1.0, 2.0, 3.0}, y(4, 0.0);
        exchng(a, b, c, d, t, q, x, y, 3);
        chk_rel(b, 1.5, 1e-12, "exchng b=1.5");
        chk_rel(d, 2.5, 1e-12, "exchng d=2.5");
        chk_rel(q, 3.5, 1e-12, "exchng q=3.5");
        chk_rel(y[1], 1.0, 1e-12, "exchng y1=1");
        chk_rel(y[3], 3.0, 1e-12, "exchng y3=3");
    }

    std::fprintf(stderr, "[t2] ncoord (coordination number)\n");
    {
        numat = 2;
        nat = std::vector<int>{0, 6, 6};
        std::vector<int> nat_loc = nat;
        std::vector<double> rcov(95, 0.0); rcov[1] = 1.0; rcov[6] = 1.0;
        // atoms at (0,0,0) and (2,0,0) au: rco=2.0, rr=1.0 -> damp=0.5
        std::vector<std::vector<double>> xyz(4, std::vector<double>(3, 0.0));
        xyz[1][1] = 0.0; xyz[1][2] = 2.0;
        std::vector<double> cn(3, 0.0);
        ncoord(2, rcov, nat_loc, xyz, cn);
        chk_rel(cn[1], 0.5, 1e-12, "ncoord cn1=0.5");
        chk_rel(cn[2], 0.5, 1e-12, "ncoord cn2=0.5");
    }

    std::fprintf(stderr, "[t3] get_dC6_dCNij (single reference point)\n");
    {
        const int maxc = 5, max_elem = 94;
        c6ab_t c6ab(max_elem + 1,
            std::vector<std::vector<std::vector<std::vector<double>>>>(max_elem + 1,
                std::vector<std::vector<std::vector<double>>>(maxc + 1,
                    std::vector<std::vector<double>>(maxc + 1, std::vector<double>(4, 0.0)))));
        // one reference: izi=6, izj=6, a=1, b=1 -> c6=10, cri=1, crj=1
        c6ab[6][6][1][1][1] = 10.0;
        c6ab[6][6][1][1][2] = 1.0;
        c6ab[6][6][1][1][3] = 1.0;
        double c6check = 0.0, dc6i = 0.0, dc6j = 0.0;
        get_dC6_dCNij(maxc, max_elem, c6ab, 1, 1, 1.0, 1.0, 6, 6, c6check, dc6i, dc6j);
        chk_rel(c6check, 10.0, 1e-12, "getc6 c6check=10");
        chk_rel(dc6i, 0.0, 1e-12, "getc6 dc6i=0");
        chk_rel(dc6j, 0.0, 1e-12, "getc6 dc6j=0");
    }

    std::fprintf(stderr, "[t4] gdisp end-to-end (2 atoms, zero gradient at optimal)\n");
    {
        numat = 2;
        nat = std::vector<int>{0, 6, 6};
        std::vector<double> rcov(95, 0.0); rcov[1] = 1.0; rcov[6] = 1.0;
        // r0ab[6][6] must be positive for damping; use 3.0
        std::vector<std::vector<double>> r0ab(95, std::vector<double>(95, 0.0));
        r0ab[6][6] = 3.0;
        const int maxc = 5, max_elem = 94;
        c6ab_t c6ab(max_elem + 1,
            std::vector<std::vector<std::vector<std::vector<double>>>>(max_elem + 1,
                std::vector<std::vector<std::vector<double>>>(maxc + 1,
                    std::vector<std::vector<double>>(maxc + 1, std::vector<double>(4, 0.0)))));
        c6ab[6][6][1][1][1] = 20.0;
        c6ab[6][6][1][1][2] = 1.0;
        c6ab[6][6][1][1][3] = 1.0;
        std::vector<int> mxc(95, 0); mxc[6] = 1;
        std::vector<std::vector<double>> xyz(4, std::vector<double>(3, 0.0));
        xyz[1][1] = 0.0; xyz[1][2] = 3.0;  // r = R0 = 3.0 au
        std::vector<std::vector<double>> dxyz_temp(4, std::vector<double>(3, 0.0));
        double s6 = 1.0, rs6 = 1.560, alp6 = 14.0;
        gdisp(xyz, r0ab, rs6, alp6, c6ab, s6, mxc, rcov, dxyz_temp);
        // Force balance: dxyz on atom1 must be exactly opposite atom2, and finite.
        chk(std::isfinite(dxyz_temp[1][1]) && std::isfinite(dxyz_temp[2][1]), "gdisp finite");
        chk_rel(dxyz_temp[1][1] + dxyz_temp[1][2], 0.0, 1e-9, "gdisp force balance");
        chk_rel(dxyz_temp[2][1] + dxyz_temp[2][2], 0.0, 1e-9, "gdisp y balance");
        chk(dxyz_temp[1][1] != 0.0, "gdisp nonzero force at R0");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_fail ? 1 : 0;
}
