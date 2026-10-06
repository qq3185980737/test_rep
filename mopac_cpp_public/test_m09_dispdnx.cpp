// test_m09_dispdnx.cpp — verification for disp_DnX.cpp (MOPAC 2016).
// t1: halogen-bond energy (Cl...N, 3.0 A, PM6-DH+); t2: gradient branch
// with cell bookkeeping; t3: print_post_scf_corrections threshold logic.
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "disp_DnX.h"
#include "H_bond_correction_bits.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

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
    // PM6-DH+ (not DH2X): Cl(17)-N(7), a=1.049e12 b=-9.95
    method_PM6_DH2X = false;
    id = 0; l123 = 1; l1u = l2u = l3u = 0;
    numat = 2;
    nat = std::vector<int>{0, 17, 7};
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[1][2] = 3.0;   // 3.0 A apart on x
    Vab.assign(4, 0.0); Vab[1] = 1.0; Vab[2] = 0.0; Vab[3] = 0.0;
    cell_ijk = std::vector<int>{0, 0, 0, 0};
    dxyz.assign(10, 0.0);

    std::fprintf(stderr, "[t1] disp_DnX energy (Cl...N 3.0 A)\n");
    {
        double want = 1.049e12 * std::exp(-9.95 * 3.0);
        double got = disp_DnX(false);
        chk_rel(got, want, 1e-9, "disp_DnX ClN sum");
        chk_rel(disp_DnX(false), want, 1e-9, "disp_DnX idempotent");
    }

    std::fprintf(stderr, "[t2] disp_DnX gradient branch\n");
    {
        dxyz.assign(10, 0.0);
        double got = disp_DnX(true);
        double fact = 1.049e12 * (-9.95) * std::exp(-9.95 * 3.0) / 3.0;
        chk_rel(got, 1.049e12 * std::exp(-9.95 * 3.0), 1e-9, "grad sum unchanged");
        // connected_hb overwrites Vab with the (unnormalized) coord diff:
        // Vab[1]=-3.0, Rab set to distance 3.0
        chk_rel(dxyz[1], -3.0 * fact, 1e-9, "grad dxyz[i] = Vab1*fact");
        chk_rel(dxyz[4], 3.0 * fact, 1e-9, "grad dxyz[j] = -Vab1*fact");
        chk(dxyz[2] == 0.0 && dxyz[3] == 0.0, "grad y/z untouched");
    }

    std::fprintf(stderr, "[t3] print_post_scf_corrections threshold loop\n");
    {
        P_Hbonds = 2;
        H_energy = std::vector<double>{0, -5.0, -2.0};
        H_txt = { "", "hb_strong", "hb_weak" };
        keywrd = " DISP(-4.0)";
        print_post_scf_corrections();
        chk_rel(H_energy[1], 10.0, 1e-12, "print strong removed");
        chk_rel(H_energy[2], -2.0, 1e-12, "print weak kept");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_fail ? 1 : 0;
}
