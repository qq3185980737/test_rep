// test_batchM06e.cpp  M06 deri2 relaxation chain (deri21/22/23 + CI).
// Scenario: 2 atoms x 1 AO, closed(1) + virtual(1), CI-active 2 MOs,
// nelec=0 -> lab=2 microstates, 1 state, 1 geometric variable.
// deri23 is verified with exact hand-checked values; deri22 and the deri2
// driver are exercised for structure, finiteness and determinism.
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "deri2.h"
#include "deri22.h"
#include "deri23.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "molkst_C.h"

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

// F90 NDDO kernels not in the 2016 tree: stubs keep the pipeline linkable.
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

static void setup_common() {
    norbs = 2; numat = 2; mpack = 3; numcal = 1;
    nfirst.assign({0, 1, 2});
    nlast.assign({0, 1, 2});
    nat.assign({0, 1, 1});
    keywrd = " ";
    id = 0;
    p.assign(4, 0.0); pa.assign(4, 0.0); pb.assign(4, 0.0);
    p[1] = 0.36; p[2] = 0.48; p[3] = 0.64;
    w.assign(4, 0.0);
    w[1] = 0.3;  // one nonzero 2-electron integral (l-l block, ijkl=(1,1,1,1))
    c.assign(3, std::vector<double>(3, 0.0));
    c[1][1] = 0.6; c[2][1] = 0.8;
    c[1][2] = 0.8; c[2][2] = -0.6;
    eigs.assign(3, 0.0);
    eigs[1] = -0.5; eigs[2] = 0.3;
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.5; coord[2][2] = 0.0; coord[3][2] = 0.0;
    nmos = 2; lab = 2; nstate = 1; nelec = 0; maxci = 2; nmeci = 2;
    nbo[1] = 1; nbo[2] = 0; nbo[3] = 1;
    nopen = nbo[1] + nbo[2];
    fract = 0.5;
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
    std::printf("M06e batch - deri2 relaxation chain\n");
    setup_common();

    // ---- 1) deri23: exact hand-check ----
    // norbs=2, nbo={1,0,1}, nelec=0, nmos=2, nopen=1.
    // PART1: loop1 ninit=1 nend=1 n1=1 n2=1 (i<=ninit, only cmo[1][1]=0);
    // loop2 nbo2=0 -> n2<n1 continue; loop3 ninit=2 nend=2 n1=2 n2=2,
    // i<=ninit; ncol=1 but n2=2 not < norbs=2 -> skipped. l stays 1.
    // emo[1..2] = fci[1..2].
    // PART2: closed-open empty; virtual-closed: j=1, i=2, scal=0.5,
    // com = f[1]*0.5 = 0.25 -> cmo[2][1]=-0.25, cmo[1][2]=0.25.
    std::vector<double> fcol(2, 0.0), fdcol(2, 0.0), fcicol(4, 0.0);
    fcol[1] = 0.5; fdcol[1] = 0.1; fcicol[1] = 0.2; fcicol[2] = -0.1;
    std::vector<std::vector<double>> cmo(3, std::vector<double>(3, 0.0));
    std::vector<double> emo(4, 0.0);
    deri23(fcol, fdcol, eigs, fcicol, cmo, emo, 1, 1, 1);
    CHK_D(cmo[2][1], -0.25, 1e-12, "deri23 virtual-closed cmo(2,1)");
    CHK_D(cmo[1][2], 0.25, 1e-12, "deri23 virtual-closed cmo(1,2)");
    CHK_D(cmo[1][1], 0.0, 1e-12, "deri23 closed block cmo(1,1)=0");
    CHK_D(cmo[2][2], 0.0, 1e-12, "deri23 active diag cmo(2,2)=0");
    CHK_D(emo[1], 0.2, 1e-12, "deri23 emo(1)=fci(l)");
    CHK_D(emo[2], -0.1, 1e-12, "deri23 emo(2)=fci(l+1)");

    // ---- 2) deri22: structure + determinism ----
    // b(1)=1; scalar(1)=2 -> STEP0 b*=2, PART2 b/=2 restores b=1.
    // Blocks: nbo2=0 skip; nbo3=1 & nbo1=1 virtual-closed mxm + closed-virtual
    // mxmt; work filled; ab/fci finite.
    std::vector<double> b(2, 0.0), ab(2, 0.0), fci2(70, 0.0);
    std::vector<double> foc2(6, 0.0);
    b[1] = 1.0;
    std::vector<double> diag(4, 0.0), scalar(2, 0.0);
    diag[1] = 3.0; scalar[1] = 2.0;
    std::vector<std::vector<double>> work2d(3, std::vector<double>(3, 0.0));
    deri22(c, b, 1, work2d, foc2, ab, 1, 1, fci2, 1, w, diag, scalar, 1);
    CHK_D(b[1], 1.0, 1e-12, "deri22 restores b after scale/unscale");
    CHECK(std::isfinite(ab[1]), "deri22 ab finite");
    CHECK(std::isfinite(fci2[1]) && std::isfinite(fci2[2]), "deri22 fci finite");
    CHECK(std::isfinite(work2d[1][1]) && std::isfinite(work2d[2][2]),
          "deri22 work filled");
    // Determinism: repeat with fresh output buffers, identical input state.
    std::vector<double> b2(2, 0.0), ab2(2, 0.0), fci2b(70, 0.0), foc2b(6, 0.0);
    b2[1] = 1.0;
    std::vector<std::vector<double>> work2d2(3, std::vector<double>(3, 0.0));
    deri22(c, b2, 1, work2d2, foc2b, ab2, 1, 1, fci2b, 1, w, diag, scalar, 1);
    CHK_D(ab[1], ab2[1], 1e-12, "deri22 deterministic ab");
    CHK_D(fci2[1], fci2b[1], 1e-12, "deri22 deterministic fci");

    // ---- 3) deri2 driver: end-to-end relaxation ----
    // minear=1, ninear=1, nvar_nvo=1. f[1][1]=0.5 (off-diag block),
    // fd[1][1]=0.1. deri21: A'A=0.25 -> pnert=0.5 -> b(1)=1.0, ilast=1.
    // STEP2: deri22+mxm build bab(1,1); osinv; residual loop; STEP3:
    // deri23 -> dijkl2 -> mecid/mecih -> supdot -> dxyzr update.
    const int minear = 1, ninear = 1, nvar_nvo = 1;
    std::vector<std::vector<double>> f2d(2, std::vector<double>(2, 0.0));
    std::vector<std::vector<double>> fd2d(2, std::vector<double>(2, 0.0));
    std::vector<std::vector<double>> fci2d(61, std::vector<double>(2, 0.0));
    f2d[1][1] = 0.5;
    fd2d[1][1] = 0.1;
    std::vector<double> dxyzr(1, 0.0), dxyzr2(1, 0.0);
    std::vector<double> work(128, 0.0), work2(128, 0.0);
    deri2(minear, f2d, fd2d, fci2d, ninear, nvar_nvo, dxyzr, 1e-8, diag,
          scalar, work);
    CHECK(std::isfinite(dxyzr[0]), "deri2 dxyzr finite");
    CHECK(std::isfinite(f2d[1][1]), "deri2 f unscaled solution finite");
    deri2(minear, f2d, fd2d, fci2d, ninear, nvar_nvo, dxyzr2, 1e-8, diag,
          scalar, work2);
    CHK_D(dxyzr[0], dxyzr2[0], 1e-9, "deri2 deterministic dxyzr");

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
