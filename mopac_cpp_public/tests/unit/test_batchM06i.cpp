// test_batchM06i.cpp  M06 deritr (finite-difference internal-coordinate
// derivatives via full SCF). Scene: 2 atoms, 1 bond variable, closed shell,
// non-PRECISE (backward difference, delta=0.0002, xderiv=5000). The SCF
// path (hcore+iter) is the partial port; assertions keep to structural
// invariants: errfn finite, energy globals restored, geometry restored.
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "derivs_C.h"
#include "deritr.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "molmec_C.h"
#include "MOZYME_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

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
using namespace molkst_C;

// SCF kernels (hcore/iter and their 12-symbol closure) belong to the
// M03 SCF subsystem, not yet closed out; deritr control flow is validated
// here with stubs. The numerical SCF-derivative chain is a known gap.
void hcore() {}
void iter(double&, bool, bool) {}
void hcore_for_MOZYME() {}
void iter_for_MOZYME(double&) {}
void post_scf_corrections(double&, bool) {}

// F90 NDDO kernels not in the 2016 tree: stubs keep the pipeline linkable.
extern "C" void diat_(int*, int*, double*, double* smat) {
    for (int i = 0; i < 81; ++i) smat[i] = 0.0;
}
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj,
                        double* w, int* kr, double* enuc) {
    *kr = 1; w[1] = 0.0; *enuc = 0.0;
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
    std::printf("M06i batch - deritr\n");
    numcal = 1;
    norbs = 2; numat = 2; natoms = 2; mpack = 3; nvar = 1; ndep = 0;
    l1u = 0; l2u = 0; l3u = 0; l123 = 1; id = 0; n2elec = 0;
    nclose = 1; nopen = 1; nelecs = 2; fract = 2.0;
    cutofp = 9.0;
    keywrd = " ";
    use_ref_geo = false;
    method_pm6 = false;
    method_pm7 = false;
    N_3_present = false;
    Si_O_H_present = false;
    density = 0.0;
    mozyme = false;
    funcon_C::fpc_9 = 14.4; funcon_C::ev = 27.21; funcon_C::a0 = 0.529177;
    nfirst.assign({0, 1, 2});
    nlast.assign({0, 1, 2});
    nat.assign({0, 1, 1});
    p.assign(10, 0.0);
    pa.assign(10, 0.0);
    pb.assign(10, 0.0);
    pa[1] = 0.5; pa[2] = 0.25; pa[3] = 0.125;
    tvec.assign(4, std::vector<double>(4, 0.0));
    nbonds.assign(5, 0);
    ibonds.assign(5, std::vector<int>(5, 0));
    geo.assign(4, std::vector<double>(4, 0.0));
    geoa.assign(4, std::vector<double>(4, 0.0));
    coord.assign(4, std::vector<double>(4, 0.0));
    labels.assign(4, 0);
    na.assign(4, 0);
    nb.assign(4, 0);
    nc.assign(4, 0);
    step = 0.0;
    for (int i = 0; i <= 4; ++i) parameters_C::tore[i] = 1.0;
    molmec_C::nnhco = 0;
    cosmo_C::useps = false;
    numat_old = 0;
    loc.assign(3, std::vector<int>(4, 0));
    loc[1][1] = 2; loc[2][1] = 1;  // atom 2, bond length
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    geo[1][2] = 1.5; geo[2][2] = 0.0; geo[3][2] = 0.0;
    // Scratch for the SCF path.
    derivs_C::fmat.assign(10, 0.0);
    derivs_C::hmat.assign(10, 0.0);
    derivs_C::wmat.assign(10, 0.0);
    errfn.assign(4, 0.0);

    escf = 1.234;
    enuclr = 0.567;
    elect = 0.333;

    deritr();
    CHECK(std::isfinite(errfn[1]), "deritr errfn finite");
    CHK_D(escf, 1.234, 1e-12, "deritr escf restored");
    CHK_D(enuclr, 0.567, 1e-12, "deritr enuclr restored");
    CHK_D(elect, 0.333, 1e-12, "deritr elect restored");
    // F90 semantics: the loop ends with geo(l,k)=xstore-delta (last
    // perturbed point); the caller rebuilds the geometry.
    CHK_D(geo[1][2], 1.5 - 0.0002, 1e-9,
          "deritr geo at x-delta (F90 semantic)");
    // Non-PRECISE: errfn=(E(x)-E(x-delta))*const*5000; with the stub NDDO
    // kernels the two SCF energies coincide, so errfn is zero but finite.
    CHECK(std::fabs(errfn[1]) < 1e8, "deritr errfn bounded");

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
