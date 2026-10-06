// test_m10_small.cpp — M10 small-property batch:
// dfield + nuchar + volume + polar_C + prtpka.
// Acceptance: ALL PASS + ASan clean + zero warnings.
#include "dfield.h"
#include "nuchar.h"
#include "volume.h"
#include "polar_C.h"
#include "prtpka.h"

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "cosmo_C.h"
#include "chanel_C.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using common_arrays_C::nat;
using common_arrays_C::p;
using common_arrays_C::dxyz;
using common_arrays_C::coord;
using common_arrays_C::grad;
using common_arrays_C::xparam;
using molkst_C::numat;
using molkst_C::efield;
using molkst_C::escf;
using molkst_C::mozyme;
using molkst_C::emin;
using molkst_C::mpack;
using molkst_C::moperr;
using parameters_C::tore;
using parameters_C::uss;
using parameters_C::betas;
using parameters_C::zs;
using parameters_C::gss;
using cosmo_C::iseps;
using cosmo_C::useps;
using cosmo_C::lpka;
using chanel_C::iw;
using funcon_C::ev;
using funcon_C::a0;
using funcon_C::fpc_9;

// ---- Stubs for heavy external dependencies (host provides the real ones) --
void cosini(bool) {}
void coscav() {}
void mkbmat() {}
void moldat(int) {}
void calpar() {}
void compfg(const std::vector<double>&, bool, double&, bool,
            std::vector<double>&, bool) {}
void switch_method() {}
void chrge(const std::vector<double>&, std::vector<double>& q) {
    for (size_t i = 1; i < q.size(); ++i) q[i] = 0.0;  // zero electron density
}
namespace linear_cosmo {
void ini_linear_cosmo() {}
void coscavz(std::vector<double>&, int) {}
}  // namespace linear_cosmo

static int fails = 0;
#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        if (cond) {                                                         \
            std::printf("PASS: %s\n", msg);                                 \
        } else {                                                            \
            std::printf("FAIL: %s (line %d)\n", msg, __LINE__);             \
            ++fails;                                                        \
        }                                                                   \
    } while (0)
#define CHECK_NEAR(a, b, tol, msg)                                          \
    do {                                                                    \
        double _a = (a), _b = (b);                                          \
        if (std::fabs(_a - _b) <= (tol)) {                                  \
            std::printf("PASS: %s\n", msg);                                 \
        } else {                                                            \
            std::printf("FAIL: %s got %g want %g (line %d)\n", msg, _a, _b, \
                        __LINE__);                                          \
            ++fails;                                                        \
        }                                                                   \
    } while (0)

// 1) dfield: single C atom in an applied field; zero density -> net charge 4.
void t_dfield() {
    numat = 1;
    nat.assign(4, 0);
    nat[1] = 6;
    p.assign(16, 0.0);
    dxyz.assign(16, 0.0);
    efield[1] = 1.0;
    efield[2] = 2.0;
    efield[3] = 3.0;
    dfield();
    double fldcon = ev / a0 * fpc_9;
    double q2 = tore[6];  // net charge = tore (zero electron density)
    CHECK_NEAR(dxyz[1], efield[1] * q2 * fldcon, 1e-6, "dfield: dxyz x");
    CHECK_NEAR(dxyz[2], efield[2] * q2 * fldcon, 1e-6, "dfield: dxyz y");
    CHECK_NEAR(dxyz[3], efield[3] * q2 * fldcon, 1e-6, "dfield: dxyz z");
}

// 2) nuchar: extract numbers incl. Fortran D-exponent.
void t_nuchar() {
    char line[] = "1.5, 2.5e-2 3.5D-3";
    double value[41] = {0.0};
    int nvalue = 0;
    nuchar(line, (int)std::strlen(line), value, nvalue);
    CHECK(nvalue == 3, "nuchar: 3 numbers found");
    CHECK_NEAR(value[1], 1.5, 1e-12, "nuchar: value1");
    CHECK_NEAR(value[2], 0.025, 1e-12, "nuchar: value2");
    CHECK_NEAR(value[3], 0.0035, 1e-12, "nuchar: value3 (D exponent)");
}

// 3) volume: 1D length / 2D area / 3D volume.
void t_volume() {
    double v1[9] = {3.0, 0, 0, 0, 0, 0, 0, 0, 0};
    CHECK_NEAR(volume(v1, 1), 3.0, 1e-12, "volume: ndim=1");
    double v2[9] = {3.0, 0, 0, 0, 4.0, 0, 0, 0, 0};
    CHECK_NEAR(volume(v2, 2), 12.0, 1e-12, "volume: ndim=2 rect");
    double v3[9] = {2.0, 0, 0, 0, 2.0, 0, 0, 0, 2.0};
    CHECK_NEAR(volume(v3, 3), 8.0, 1e-12, "volume: ndim=3 cube");
    double vt[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    CHECK_NEAR(volume(vt, 3), 1.0, 1e-12, "volume: unit cell");
    // tilted cell: a=(1,1,0) b=(1,0,1) c=(0,1,1) -> det = -2 -> |..| = 2
    double vt2[9] = {1, 1, 0, 1, 0, 1, 0, 1, 1};
    CHECK_NEAR(volume(vt2, 3), 2.0, 1e-12, "volume: tilted cell");
}

// 4) polar_C module data.
void t_polar_C() {
    polar_C::omega = 0.5;
    polar_C::alpavg = 12.3;
    polar_C::dummy = -1.0;
    CHECK_NEAR(polar_C::omega, 0.5, 1e-15, "polar_C: omega");
    CHECK_NEAR(polar_C::alpavg, 12.3, 1e-15, "polar_C: alpavg");
    CHECK_NEAR(polar_C::dummy, -1.0, 1e-15, "polar_C: dummy");
}

// 5a) prtpka: no O-H -> nh=0 -> no=0.
void t_prtpka_none() {
    numat = 1;
    nat.assign(4, 0);
    nat[1] = 6;
    coord.assign(4, std::vector<double>(2, 0.0));
    p.assign(16, 0.0);
    grad.assign(16, 0.0);
    xparam.assign(16, 0.0);
    mozyme = false;
    mpack = 15;
    moperr = false;
    escf = 0.0;
    emin = 0.0;
    iseps = useps = lpka = false;
    iw = 0;  // silence output
    int s[5] = {0}, u[5] = {0};
    double sp[5] = {0}, up[5] = {0};
    int no = -1;
    prtpka(s, sp, u, up, no);
    CHECK(no == 0, "prtpka: no O-H -> no=0");
}

// 5b) prtpka: water (3 atoms, two O-H) -> both H ionizable, nh=2, same pKa.
void t_prtpka_water() {
    numat = 3;
    nat.assign(4, 0);
    nat[1] = 8;
    nat[2] = 1;
    nat[3] = 1;
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1] = 0.0;
    coord[2][1] = 0.0;
    coord[3][1] = 0.0;
    coord[1][2] = 0.96;
    coord[2][2] = 0.0;
    coord[3][2] = 0.0;
    double ang = 104.5 * 3.14159265358979323846 / 180.0;
    coord[1][3] = 0.96 * std::cos(ang);
    coord[2][3] = 0.96 * std::sin(ang);
    coord[3][3] = 0.0;
    p.assign(16, 0.0);
    grad.assign(16, 0.0);
    xparam.assign(16, 0.0);
    mozyme = false;
    mpack = 15;
    moperr = false;
    escf = 0.0;
    emin = 0.0;
    iseps = useps = lpka = false;
    iw = 0;
    int s[5] = {0}, u[5] = {0};
    double sp[5] = {0}, up[5] = {0};
    int no = -1;
    prtpka(s, sp, u, up, no);
    CHECK(no == 2, "prtpka: water -> no=2");
    // pKa = q*c1 + dist*c2 + c3 with q(H)=1, dist=0.96
    double expect = 1.0 * (-288.0530769) + 0.96 * 28.6888717 + 89.1172382;
    CHECK_NEAR(sp[1], expect, 1e-6, "prtpka: pKa_sorted[1] formula");
    CHECK_NEAR(sp[2], expect, 1e-6, "prtpka: pKa_sorted[2] formula");
    CHECK(s[1] == 2 && s[2] == 3, "prtpka: ipKa_sorted = {2,3}");
    // after sorting, the selected entries are reset to the converted value
    CHECK_NEAR(up[1], expect, 1e-6, "prtpka: unsorted reset to converted");
}

// 5c) Parameters_for_PKA: constants + parameter-table overrides.
void t_parameters_for_pka() {
    double c1, c2, c3;
    Parameters_for_PKA(c1, c2, c3);
    CHECK_NEAR(c1, -288.0530769, 1e-9, "PKA: c1");
    CHECK_NEAR(c2, 28.6888717, 1e-9, "PKA: c2");
    CHECK_NEAR(c3, 89.1172382, 1e-9, "PKA: c3");
    CHECK_NEAR(uss[1], -9.3555403, 1e-9, "PKA: uss(H)");
    CHECK_NEAR(betas[1], -2.6789813, 1e-9, "PKA: betas(H)");
    CHECK_NEAR(zs[1], 1.2458568, 1e-9, "PKA: zs(H)");
    CHECK_NEAR(gss[1], 14.5964443, 1e-9, "PKA: gss(H)");
    CHECK_NEAR(uss[6], -49.9390029, 1e-9, "PKA: uss(C)");
    CHECK_NEAR(betas[6], -12.6446141, 1e-9, "PKA: betas(C)");
    CHECK_NEAR(zs[6], 1.6412075, 1e-9, "PKA: zs(C)");
    CHECK_NEAR(zs[8], 4.3991035, 1e-9, "PKA: zs(O)");
    CHECK_NEAR(gss[8], 16.3875625, 1e-9, "PKA: gss(O)");
    CHECK_NEAR(zs[17], 4.2327248, 1e-9, "PKA: zs(Cl)");
    CHECK_NEAR(zs[53], 4.9453882, 1e-9, "PKA: zs(I)");
}

int main() {
    std::printf("== M10 small-property batch ==\n");
    std::printf("[t1 dfield]\n");
    t_dfield();
    std::printf("[t2 nuchar]\n");
    t_nuchar();
    std::printf("[t3 volume]\n");
    t_volume();
    std::printf("[t4 polar_C]\n");
    t_polar_C();
    std::printf("[t5 prtpka]\n");
    t_prtpka_none();
    t_prtpka_water();
    t_parameters_for_pka();
    if (fails == 0) {
        std::printf("ALL PASS: M10 small-property group\n");
        return 0;
    }
    std::printf("FAILURES: %d\n", fails);
    return 1;
}
