// test_batchM05b.cpp — M05b subsystem batch: symopr/rotmol/xyzcry/bldsym/dtran2/dtrans/frame/rotate
// ASCII only. MSVC + ASan.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "symopr.h"
#include "rotmol.h"
#include "xyzcry.h"
#include "bldsym.h"
#include "dtran2.h"
#include "dtrans.h"
#include "frame.h"
#include "rotate.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

using namespace symmetry_C;
using namespace molkst_C;
using namespace common_arrays_C;

// ---- stubs for permanent Fortran gaps (NDDO core) and axis (test-controlled) ----
extern "C" void axis_(double&, double&, double&, double* rot) {
    for (int k = 0; k < 16; ++k) rot[k] = 0.0;
    rot[1 * 4 + 1] = 1.0;  // rot[1][1]
    rot[2 * 4 + 2] = 1.0;  // rot[2][2]
    rot[3 * 4 + 3] = 1.0;  // rot[3][3]
}
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj,
                        double* w, int* kr, double* enuc) {
    *kr = 1; *enuc = 0.0; w[0] = 0.0;
}
extern "C" void elenuc_(int*, int*, int*, int*, double* en) {
    for (int i = 0; i < 171; ++i) en[i] = 0.0;
}
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double* rij, int* ni, int* nj) {}

static int nPass = 0, nFail = 0;
static void chk(const char* name, double got, double exp, double tol) {
    if (std::fabs(got - exp) <= tol) { ++nPass; }
    else { ++nFail; std::printf("FAIL %s: got %.12g exp %.12g\n", name, got, exp); }
}

int main() {
    // ---- 1,2: symopr (column-major r) ----
    // r = [[1,2,0],[3,4,0],[0,0,1]]
    double r[9] = {1, 3, 0, 2, 4, 0, 0, 0, 1};
    double loc[3] = {1, 2, 3};
    symopr(1, loc, 1, r);
    chk("symopr+1 x", loc[0], 7.0, 1e-12);
    chk("symopr+1 y", loc[1], 10.0, 1e-12);
    chk("symopr+1 z", loc[2], 3.0, 1e-12);
    loc[0] = 1; loc[1] = 2; loc[2] = 3;
    symopr(1, loc, -1, r);
    chk("symopr-1 x", loc[0], 5.0, 1e-12);
    chk("symopr-1 y", loc[1], 11.0, 1e-12);
    chk("symopr-1 z", loc[2], 3.0, 1e-12);

    // ---- 3: rotmol: 90 deg rotation, coord (1,0,0) -> (0,-1,0) ----
    {
        double rr[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        double c[3] = {1, 0, 0};
        rotmol(1, c, 1.0, 0.0, 1, 2, rr);
        chk("rotmol x", c[0], 0.0, 1e-12);
        chk("rotmol y", c[1], -1.0, 1e-12);
        chk("rotmol z", c[2], 0.0, 1e-12);
    }

    // ---- 4,5: xyzcry ----
    {
        double tv[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        double dx[3] = {1, 2, 3};
        xyzcry(tv, 1, dx, 0);
        chk("xyzcry orth x", dx[0], 1.0, 1e-12);
        chk("xyzcry orth y", dx[1], 2.0, 1e-12);
        chk("xyzcry orth z", dx[2], 3.0, 1e-12);
    }
    {
        double tv[9] = {2, 0, 0, 0, 3, 0, 0, 0, 4};
        double dx[3] = {2, 6, 12};
        xyzcry(tv, 1, dx, 0);
        chk("xyzcry diag x", dx[0], 4.0, 1e-12);
        chk("xyzcry diag y", dx[1], 18.0, 1e-12);
        chk("xyzcry diag z", dx[2], 48.0, 1e-12);
    }

    // ---- 6,7,8: bldsym ----
    {
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 21; ++k) elem[i][j][k] = 0.0;
        bldsym(1, 1);
        chk("bldsym1 11", elem[1][1][1], 1.0, 1e-12);
        chk("bldsym1 22", elem[2][2][1], -1.0, 1e-12);
        chk("bldsym1 33", elem[3][3][1], -1.0, 1e-12);
        chk("bldsym1 12", elem[1][2][1], 0.0, 1e-12);
        chk("bldsym1 21", elem[2][1][1], 0.0, 1e-12);
    }
    {
        bldsym(9, 2);
        chk("bldsym9 11", elem[1][1][2], std::cos(6.2831853071796 / 4.0), 1e-12);
        chk("bldsym9 22", elem[2][2][2], std::cos(6.2831853071796 / 4.0), 1e-12);
        chk("bldsym9 21", elem[2][1][2], 1.0, 1e-12);
        chk("bldsym9 12", elem[1][2][2], -1.0, 1e-12);
    }
    {
        bldsym(20, 3);
        chk("bldsym20 12", elem[1][2][3], 1.0, 1e-12);
        chk("bldsym20 21", elem[2][1][3], 1.0, 1e-12);
        chk("bldsym20 33", elem[3][3][3], -1.0, 1e-12);
        chk("bldsym20 11", elem[1][1][3], 0.0, 1e-12);
    }

    // ---- 9: dtran2 with identity r -> identity t ----
    {
        std::vector<std::vector<double>> r3(3, std::vector<double>(3, 0.0));
        r3[0][0] = 1; r3[1][1] = 1; r3[2][2] = 1;
        std::vector<std::vector<std::vector<double>>> t(12,
            std::vector<std::vector<double>>(5, std::vector<double>(5, 0.0)));
        dtran2(r3, t, 1);
        for (int a = 1; a <= 5; ++a)
            for (int b = 1; b <= 5; ++b)
                chk("dtran2 id", t[1][a - 1][b - 1], (a == b ? 1.0 : 0.0), 1e-12);
    }

    // ---- 10: dtrans with identity -> d unchanged ----
    {
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 21; ++k) elem[i][j][k] = 0.0;
        elem[1][1][1] = 1; elem[2][2][1] = 1; elem[3][3][1] = 1;
        nclass = 1;
        std::vector<double> d(6, 0.0);
        d[1] = 1.0;
        bool first = true;
        std::vector<std::vector<double>> orient(3, std::vector<double>(3, 0.0));
        orient[0][0] = 1; orient[1][1] = 1; orient[2][2] = 1;
        dtrans(d, 1, first, orient);
        chk("dtrans d1", d[1], 1.0, 1e-12);
        chk("dtrans d2", d[2], 0.0, 1e-12);
        chk("dtrans d3", d[3], 0.0, 1e-12);
    }

    // ---- 11: rotate degenerate (rij < 2e-5) -> all zero, kr=1 ----
    {
        double xi[3] = {0, 0, 0}, xj[3] = {0, 0, 0};
        double w[2026] = {}; for (int i = 0; i < 2026; ++i) w[i] = 123.0;
        int kr = 5;
        double e1b[45] = {}, e2a[45] = {};
        double enuc = 5.0;
        rotate(1, 1, xi, xj, w, kr, e1b, e2a, enuc);
        chk("rotate kr", (double)kr, 5.0, 0.0);  // kr unchanged on small-rij return
        chk("rotate enuc", enuc, 0.0, 0.0);
        double mx = 0.0;
        for (int i = 1; i < 2026; ++i) if (std::fabs(w[i]) > mx) mx = std::fabs(w[i]);
        chk("rotate w zero", mx, 0.0, 0.0);
    }

    // ---- 12: frame (axis_ = identity) -> positive diagonal accumulations ----
    {
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 21; ++k) elem[i][j][k] = 0.0;
        numat = 1;
        id = 0;
        itemp_2 = 0;
        coord.resize(4);
        for (auto& row : coord) row.resize(8, 0.0);
        atmass.resize(8, 1.0);
        std::vector<double> fmat(22, 0.0);
        frame(fmat, 1, 0);
        chk("frame f1", fmat[1], 41000.0, 1e-9);
        chk("frame f2", fmat[2], 0.0, 1e-9);
        chk("frame f3", fmat[3], 42000.0, 1e-9);
    }

    std::printf("M05b: %d PASS / %d FAIL\n", nPass, nFail);
    return nFail == 0 ? 0 : 1;
}
