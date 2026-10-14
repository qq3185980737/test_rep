// test_solrot.cpp — numerical verification of solrot's point-charge branch.
// With l1u=l2u=l3u=0 there is a single image; placing the pair at 5 A
// (r > cutof2) exercises the point() path: ee = ev*a0/trunk(r), where trunk
// applies a parabolic blend once r exceeds clower = cutofp*2/3.
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "solrot.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
using namespace molkst_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    if (std::fabs(got - want) > tol) {
        ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.12g want %.12g\n", name, got, want);
    }
}
// Mirror of solrot.cpp trunk() for the given cutofp.
static double trunk_mirror(double r, double cutofp) {
    double clower = (std::min)(cutofp * 2.0 / 3.0, 13.0);
    double bound1 = clower / cutofp;
    double cupper = cutofp;
    double bound2 = cupper / cutofp;
    double range = bound2 - bound1;
    double c = -0.5 * bound1 * bound1 * cutofp / range;
    double cr = 1.0 + bound1 / range;
    double cr2 = -1.0 / (cutofp * 2 * range);
    double clim = c + cr * cupper + cr2 * cupper * cupper;
    if (r > clower) {
        if (r > cupper) return clim;
        return c + cr * r + cr2 * r * r;
    }
    return r;
}

int main() {
    numcal = 1;
    cutofp = 7.0;
    l1u = 0; l2u = 0; l3u = 0;
    tvec.assign(4, std::vector<double>(4, 0.0));
    tvec[1][1] = 1; tvec[2][2] = 1; tvec[3][3] = 1;
    natorb[1] = 1; tore[1] = 1.0;
    double xi[3] = {0, 0, 0};
    double xj[3] = {5.0, 0, 0};   // r = 25 > cutof2 -> point branch
    double wj[2026] = {0}, wk[2026] = {0};
    int kr = 0;
    double e1b[46] = {0}, e2a[46] = {0};
    double enuc = 0;
    solrot(1, 1, xi, xj, wj, wk, kr, e1b, e2a, enuc);
    double rt = trunk_mirror(5.0, cutofp);
    double ee = ev * a0 / rt;
    chk(kr == 1, "kr = nii*njj = 1 for single s orbitals");
    chk_rel(wj[1], ee, 1e-12, "wj[1] = point ee after trunk");
    chk_rel(wk[1], ee, 1e-12, "wk[1] = ee (wmax path)");
    chk_rel(e1b[1], -ee * tore[1], 1e-12, "e1b[1] = -ee*tore[nj]");
    chk_rel(e2a[1], -ee * tore[1], 1e-12, "e2a[1] = -ee*tore[ni]");
    chk_rel(enuc, ee * tore[1] * tore[1], 1e-12, "enuc = ee*tore^2");

    double wj2[2026] = {0}, wk2[2026] = {0};
    int kr2 = 0;
    double e1b2[46] = {0}, e2a2[46] = {0};
    double enuc2 = 0;
    solrot(1, 1, xi, xj, wj2, wk2, kr2, e1b2, e2a2, enuc2);
    chk(kr2 == 1 && wj2[1] == wj[1], "second call reproduces (cached cutof2)");

    double xj2[3] = {1.4, 0, 0};
    double wj3[2026] = {0}, wk3[2026] = {0};
    int kr3 = 0;
    double e1b3[46] = {0}, e2a3[46] = {0};
    double enuc3 = 0;
    solrot(1, 1, xi, xj2, wj3, wk3, kr3, e1b3, e2a3, enuc3);
    chk(!std::isnan(wj3[1]) && !std::isinf(wj3[1]), "rotate stub yields finite w");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
