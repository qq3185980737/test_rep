// test_batchM06c2.cpp   M06 dcart    + derp/dihed     
#include <cstdio>
#include <cmath>
#include <vector>
#include "dcart.h"
#include "dihed.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "cosmo_C.h"
#include "molmec_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

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
// diegrd stub: useps=false in tests, so the COSMO gradient path never runs.
void diegrd(double*) {}

extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double*, int*, int*) {
    *enuc = 0.0;
}

static void test_derp() {
    std::printf("Test 1: derp three-branch numerics\n");
    molkst_C::clower = 2.0;
    molkst_C::cupper = 3.0;
    molkst_C::cutofp = 3.5;
    molkst_C::numcal = 1;
    // r < clower: -1/r^2 = -1/2.25
    CHK_D(derp(1.5), -1.0 / 2.25, 1e-12, "derp r<clower");
    // clower < r < cupper: -(cr+2*cr2*r)/(c+cr*r+cr2*r^2)^2
    // bound1=2/3.5, bound2=3/3.5, range=1/3.5, c=-2, cr=3, cr2=-0.5
    // r=2.5: den = -2+7.5-3.125 = 2.375 ; derp = -0.5/2.375^2
    CHK_D(derp(2.5), -0.5 / (2.375 * 2.375), 1e-12, "derp mid-branch");
    CHK_D(derp(4.0), 0.0, 1e-12, "derp r>cupper");
}

static void test_dihed() {
    std::printf("Test 2: dihed angles\n");
    std::vector<std::vector<double>> xyz(4, std::vector<double>(5, 0.0));
    // case A: i=(1,0,0) j=(0,0,0) k=(0,1,0) l=(0,0,1) -> -90 deg -> 3pi/2
    xyz[1][1] = 1; xyz[2][1] = 0; xyz[3][1] = 0;
    xyz[1][2] = 0; xyz[2][2] = 0; xyz[3][2] = 0;
    xyz[1][3] = 0; xyz[2][3] = 1; xyz[3][3] = 0;
    xyz[1][4] = 0; xyz[2][4] = 0; xyz[3][4] = 1;
    double angle = 0.0;
    dihed(xyz, 1, 2, 3, 4, angle);
    CHK_D(angle, acos(-1.0) / 2.0, 1e-12, "dihed A pi/2 (MOPAC convention)");
    // case B: coplanar zigzag -> 0
    xyz[1][1] = 0; xyz[2][1] = 0; xyz[3][1] = 0;
    xyz[1][2] = 1; xyz[2][2] = 0; xyz[3][2] = 0;
    xyz[1][3] = 1; xyz[2][3] = 1; xyz[3][3] = 0;
    xyz[1][4] = 0; xyz[2][4] = 1; xyz[3][4] = 0;
    dihed(xyz, 1, 2, 3, 4, angle);
    CHK_D(angle, 0.0, 1e-12, "dihed B coplanar 0");
    // case C: perpendicular -> pi/2
    xyz[1][4] = 1; xyz[2][4] = 1; xyz[3][4] = 1;
    dihed(xyz, 1, 2, 3, 4, angle);
    CHK_D(angle, 3.0 * acos(-1.0) / 2.0, 1e-12, "dihed C 3pi/2 (MOPAC convention)");
}

static void test_dcart_structural() {
    std::printf("Test 3: dcart 2-atom structural (stub kernels -> zero dxyz)\n");
    molkst_C::numcal = 2;
    molkst_C::norbs = 2;
    molkst_C::numat = 2;
    molkst_C::id = 0;
    molkst_C::l1u = 0; molkst_C::l2u = 0; molkst_C::l3u = 0;
    molkst_C::l123 = 1;
    molkst_C::cutofp = 9.0;
    molkst_C::keywrd = " ";
    molkst_C::use_ref_geo = false;
    molkst_C::method_pm6 = false;
    molkst_C::method_pm7 = false;
    molkst_C::N_3_present = false;
    molkst_C::Si_O_H_present = false;
    molkst_C::density = 0.0;
    molkst_C::mozyme = false;
    mode = 0;
    funcon_C::fpc_9 = 14.4; funcon_C::ev = 27.21; funcon_C::a0 = 0.529177;
    common_arrays_C::nfirst.assign({0, 1, 2});
    common_arrays_C::nlast.assign({0, 1, 2});
    common_arrays_C::nat.assign({0, 1, 1});
    common_arrays_C::p.assign(10, 0.0);
    common_arrays_C::pa.assign(10, 0.0);
    common_arrays_C::pb.assign(10, 0.0);
    common_arrays_C::tvec.assign(4, std::vector<double>(4, 0.0));
    common_arrays_C::nbonds.assign(5, 0);
    common_arrays_C::ibonds.assign(5, std::vector<int>(5, 0));
    common_arrays_C::geo.assign(4, std::vector<double>(4, 0.0));
    common_arrays_C::geoa.assign(4, std::vector<double>(4, 0.0));
    for (int i = 0; i <= 4; ++i) parameters_C::tore[i] = 1.0;
    molmec_C::nnhco = 0;
    cosmo_C::useps = false;
    molkst_C::numat_old = 0;

    std::vector<std::vector<double>> coord(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.5; coord[2][2] = 0.0; coord[3][2] = 0.0;
    std::vector<std::vector<double>> dxyz(4, std::vector<double>(3, 0.0));
    dcart(coord, dxyz);
    for (int k = 1; k <= 3; ++k)
        for (int i = 1; i <= 2; ++i)
            CHK_D(dxyz[k][i], 0.0, 1e-9, "dcart dxyz zero (stub kernels)");
    // coordinates untouched by the main loop path (nnhco=0)
    CHK_D(coord[1][2], 1.5, 1e-12, "dcart coord preserved");
    // repeat call (icalcn save path) must stay stable
    dcart(coord, dxyz);
    for (int k = 1; k <= 3; ++k) CHK_D(dxyz[k][1], 0.0, 1e-9, "dcart 2nd call stable");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06c2 batch — dcart/derp/dihed\n");
    test_derp();
    test_dihed();
    test_dcart_structural();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
