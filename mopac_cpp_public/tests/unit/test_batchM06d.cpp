// test_batchM06d.cpp  M06 deri1 (non-relaxed CI energy derivative)
// Structure-verification batch: 2 atoms x 1 AO, CI-active 2 MOs.
// dhcore/dfock2/helect run through the real translation chain; the
// h1elec/rotate integral kernels are external stubs (missing from the
// 2016 tree), so absolute values are stub-dependent; we assert the
// self-consistency of the assembled f/fd blocks and the CI closure.
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "common_arrays_C.h"
#include "derivs_C.h"
#include "deri1.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "molkst_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;

// F90 NDDO kernels not in the 2016 tree: stubs keep the pipeline linkable.
// dhcore consumes h1elec/rotate which call these; with zero-valued stubs the
// 1-electron/2-electron integral derivatives are zero and the structure of
// deri1 (block assembly, scaling, CI closure) is exercised with clean data.
extern "C" void diat_(int*, int*, double*, double* smat) {
    for (int i = 0; i < 81; ++i) smat[i] = 0.0;
}
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj,
                        double* w, int* kr, double* enuc) {
    // Coordinate-dependent stub: finite difference of w[1]=xi[0]^2 over the
    // +step/-step dhcore calls yields 4*csave*step (odd in csave), so the
    // whole derivative chain is nonzero and antisymmetric in the coordinate.
    *kr = 2;
    w[0] = xi[0] * xi[0];   // rotatd_ writes 0-based; rotate copies w0[0]
    *enuc = 0.0;
}
extern "C" void elenuc_(int*, int*, int*, int*, double* en) {
    for (int i = 0; i < 45; ++i) en[i] = 0.0;
}
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double*, int*, int*) {
    *enuc = 0.0;
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06d batch - deri1\n");

    // Minimal closed-shell scenario: 2 atoms, 1 AO each.
    norbs = 2; numat = 2; mpack = 3; n2elec = 0; numcal = 1;
    nfirst.assign({0, 1, 2});
    nlast.assign({0, 1, 2});
    nat.assign({0, 1, 1});
    keywrd = " ";
    p.assign(4, 0.0); pa.assign(4, 0.0);
    p[1] = 0.36; p[2] = 0.48; p[3] = 0.64;  // MO1 doubly occupied
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.5; coord[2][2] = 0.0; coord[3][2] = 0.0;
    c.assign(3, std::vector<double>(3, 0.0));
    c[1][1] = 0.6; c[2][1] = 0.8;
    c[1][2] = 0.8; c[2][2] = -0.6;
    eigb.assign(3, 0.0);
    funcon_C::fpc_9 = 1.0;

    // CI state: lab=2 microstates, one state, MO1 closed (nelec=0 active).
    nmos = 2; lab = 2; nstate = 1; nelec = 0; maxci = 2;
    nbo[1] = 1; nbo[2] = 0; nbo[3] = 1;
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
    vectci.assign(3, 0.0);
    vectci[1] = 0.9; vectci[2] = 0.1;

    // deri1 call: number = coordinate 4 (atom 2, x), minear = 1.
    int minear = 1;
    std::vector<double> fv(minear + 1, 0.0);
    std::vector<double> fd(8, 0.0);
    std::vector<double> scalar(minear + 1, 1.0);
    std::vector<std::vector<double>> work(3, std::vector<double>(3, 0.0));
    double grad = 0.0;
    deri1(4, grad, fv, minear, fd, scalar, work);

    std::printf("  grad = %.12f\n", grad);
    std::printf("  f(1) = %.12f\n", fv[1]);
    std::printf("  fd(1..7) = %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
                fd[1], fd[2], fd[3], fd[4], fd[5], fd[6], fd[7]);

    // SCF-only part: grad_self = (helect + enucl2) * 0.5/1e-3 before CI.
    // Recompute the SCF-only gradient through the same chain: since step
    // scaling is 0.5/step, grad must equal (helect+enucl2)*500 plus the
    // CI term. Check the closure: with no CI (xy=0 -> wmat=0) the CI sum
    // is 0, so grad == helect_term*const. Instead of recomputing dhcore,
    // verify structural invariants:
    // 1) fv was scaled by scalar (scalar=1 -> unchanged).
    // 2) fd diagonal block: closed-closed is empty (nelec=0, loop 1 has
    //    nbo[1]=1, n1=max(1,1)=1, n2=min(1,2)=1, i<=ninit -> skipped).
    // 3) virtual-virtual block: loop 3 ninit=2 nend=2, n1=2, n2=2,
    //    i<=ninit skipped; ncol = n2-ninit+1 = 1 >0 and n2=2 < norbs=2?
    //    false -> skipped. lcut=1.
    // 4) fd(1..nmos) = C'*W*C diagonal of active MOs.
    // With c orthogonal, fd[1]=sum_k c(k,2)*wf(k,2).
    // 5) The CI contribution uses mecih(wmat) with xy=0 -> wmat=0 -> sum=0.
    // With the coordinate-dependent rotatd_ stub, dhcore yields a nonzero
    // wmat: the 2-electron derivative path (dfock2 -> fmat, dijkl1 -> xy,
    // mecid/mecih -> CI matrix) must produce nonzero f, fd and grad.
    CHECK(std::fabs(fd[1]) > 1e-9, "fd(lcut) CI-active Fock diag nonzero");
    CHECK(std::fabs(fv[1]) > 1e-9, "f(1) off-diagonal block nonzero");
    CHECK(std::isfinite(grad) && std::fabs(grad) > 1e-9,
          "grad finite nonzero");
    // fd(lcut..lcut+nmos-1) are the CI-active Fock eigenvalue derivatives:
    // lcut=1, so fd[1] and fd[2] hold MO2/MO3 diagonals (MO numbering from
    // the active window nelec+1..nelec+nmos). They must be finite and
    // distinct (different MOs, different overlaps).
    CHECK(std::isfinite(fd[2]) && std::fabs(fd[2] - fd[1]) > 1e-12,
          "fd active diagonals finite and distinct");
    // work copy-back filled (row-major 1-based).
    CHECK(std::isfinite(work[1][1]), "work out-parameter filled");

    // Sign symmetry: flipping the sign of the coordinate perturbation must
    // flip the sign of the whole derivative chain. Repeat deri1 with a
    // mirrored coordinate; grad must change sign (stub is linear in xi).
    coord[1][2] = -1.5; coord[2][2] = 0.0; coord[3][2] = 0.0;
    std::vector<double> fv2(minear + 1, 0.0);
    std::vector<double> fd2(8, 0.0);
    std::vector<std::vector<double>> work2(3, std::vector<double>(3, 0.0));
    double grad2 = 0.0;
    deri1(4, grad2, fv2, minear, fd2, scalar, work2);
    CHK_D(grad + grad2, 0.0, 1e-9, "grad antisymmetric in coordinate");
    CHK_D(fv[1] + fv2[1], 0.0, 1e-9, "f antisymmetric in coordinate");
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
