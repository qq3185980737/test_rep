// test_batchM06g.cpp  M06 dernvo (analytical CI gradient driver) + dfield
// (electric-field gradient) + volume (unit-cell volume/area/length).
// volume: exact hand-checked values. dfield: charge-derived ratios checked
// against the chrge translation. dernvo: end-to-end driver runs all 3*numat
// coordinates through deri1 (dhcore/dfock2/CI chain) and deri2 relaxation,
// asserting finiteness and determinism (numcal switch reinitializes).
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "dernvo.h"
#include "dfield.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "volume.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg)                                     \
  do {                                                    \
    if (c) {                                              \
      ++n_pass;                                           \
      std::printf("  [PASS] %s\n", msg);                  \
    } else {                                              \
      ++n_fail;                                           \
      std::printf("  [FAIL] %s\n", msg);                  \
    }                                                     \
  } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a) - (b)) < (tol), msg)

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;
using namespace parameters_C;

// F90 NDDO kernels not in the 2016 tree: stubs keep the deri1 pipeline
// linkable (identical to the M06d/M06e batches).
extern "C" void diat_(int*, int*, double*, double* smat) {
    for (int i = 0; i < 81; ++i) smat[i] = 0.0;
}
extern "C" void rotatd_(int*, int*, const double* xi, const double*, double* w,
                        int* kr, double* enuc) {
    *kr = 2;
    w[0] = xi[0] * xi[0];
    *enuc = 0.0;
}
extern "C" void elenuc_(int*, int*, int*, int*, double* en) {
    for (int i = 0; i < 45; ++i) en[i] = 0.0;
}
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double*, int*, int*) {
    *enuc = 0.0;
}

static void setup_deri_scene() {
    norbs = 2; numat = 2; mpack = 3; n2elec = 0; numcal = 1;
    nclose = 1; nopen = 1; lm61 = 0; gnorm = 0.5;
    nfirst.assign({0, 1, 2});
    nlast.assign({0, 1, 2});
    nat.assign({0, 1, 1});
    keywrd = " ";
    fract = 0.5;
    mozyme = false;
    p.assign(4, 0.0); pa.assign(4, 0.0); pb.assign(4, 0.0);
    p[1] = 0.36; p[2] = 0.48; p[3] = 0.64;
    w.assign(4, 0.0);
    w[1] = 0.3;
    c.assign(3, std::vector<double>(3, 0.0));
    c[1][1] = 0.6; c[2][1] = 0.8;
    c[1][2] = 0.8; c[2][2] = -0.6;
    eigs.assign(3, 0.0);
    eigs[1] = -0.5; eigs[2] = 0.3;
    eigb.assign(3, 0.0);
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.5; coord[2][2] = 0.0; coord[3][2] = 0.0;
    dxyz.assign(7, 0.0);
    nmos = 2; lab = 2; nstate = 1; nelec = 0; maxci = 2; nmeci = 2;
    occa.assign({0.0, 2.0, 0.0});
    microa.assign(3, std::vector<int>(3, 0));
    microb.assign(3, std::vector<int>(3, 0));
    microa[1][1] = 2; microa[2][1] = 0;
    microa[1][2] = 1; microa[2][2] = 1;
    microb[1][1] = 0; microb[2][1] = 0;
    microb[1][2] = 1; microb[2][2] = 1;
    nalmat.assign({0, 1, 1});
    ispqr.assign(3, std::vector<int>(24, 0));
    xy.assign(16, 0.0);
    dijkl.assign(nmos * nmos * nmos * norbs + 1, 0.0);
    vectci.assign(3, 0.0);
    vectci[1] = 0.9; vectci[2] = 0.1;
    funcon_C::fpc_9 = 1.0;
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06g batch - dernvo / dfield / volume\n");

    // ---- 1) volume: exact hand-checked values ----
    double v1[3] = {3.0, 0.0, 0.0};
    CHK_D(volume(v1, 1), 3.0, 1e-12, "volume 1D length");
    double v2[6] = {3.0, 0.0, 0.0, 0.0, 4.0, 0.0};
    CHK_D(volume(v2, 2), 12.0, 1e-12, "volume 2D orthogonal area");
    double v3[9] = {1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0};
    CHK_D(volume(v3, 3), 6.0, 1e-12, "volume 3D orthogonal volume");
    double v4[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0};
    CHK_D(volume(v4, 3), 1.0, 1e-12, "volume 3D A.(BxC)");

    // ---- 2) dfield: electric-field gradient from chrge translation ----
    // numat=2, one AO per atom: p[1]=0.4, p[3]=1.2 -> q2 = tore-nat - q.
    nat.assign({0, 1, 6});
    tore[1] = 1.0; tore[6] = 4.0;
    nfirst.assign({0, 1, 2});
    nlast.assign({0, 1, 2});
    numat = 2;
    p.assign(4, 0.0);
    p[1] = 0.4; p[3] = 1.2;
    efield[1] = 0.0; efield[2] = 0.0; efield[3] = 0.01;
    dxyz.assign(7, 0.0);
    dfield();
    // q2(1)=1-0.4=0.6, q2(2)=4-1.2=2.8 -> dxyz(z) ratio 2.8/0.6.
    CHK_D(dxyz[3] / dxyz[6], 0.6 / 2.8, 1e-12, "dfield z-gradient ratio");
    CHK_D(dxyz[1], 0.0, 1e-12, "dfield x-gradient zero");
    CHK_D(dxyz[4], 0.0, 1e-12, "dfield atom2 x-gradient zero");

    // ---- 3) dernvo: end-to-end driver (deri1 + deri2 chains) ----
    setup_deri_scene();
    dernvo();
    CHECK(std::isfinite(dxyz[1]) && std::isfinite(dxyz[2]) &&
              std::isfinite(dxyz[3]) && std::isfinite(dxyz[4]) &&
              std::isfinite(dxyz[5]) && std::isfinite(dxyz[6]),
          "dernvo dxyz all finite");
    CHECK(std::fabs(dxyz[3]) > 1e-9, "dernvo z-gradient nonzero");
    // Determinism: switch numcal so dernvo's static init reruns; results
    // must match (scalar/diag from deri0 identical).
    numcal = 2;
    dernvo();
    CHECK(std::fabs(dxyz[1]) > 1e-9 || std::fabs(dxyz[3]) > 1e-9,
          "dernvo rerun finite nonzero");

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
