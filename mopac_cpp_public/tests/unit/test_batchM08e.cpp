// test_batchM08e.cpp  M08 mecid + diagi (microstate energies)
#include <cstdio>
#include <cmath>
#include <vector>
#include "mecid.h"
#include "meci_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08e batch - mecid/diagi\n");
    meci_C::nmos = 2;
    meci_C::lab = 2;
    meci_C::occa.assign({0.0, 2.0, 0.0});
    // microstates: #1 closed {2,0}|{0,0}; #2 doubly-occ {1,1}|{1,1}
    meci_C::microa.assign(3, std::vector<int>(3, 0));
    meci_C::microb.assign(3, std::vector<int>(3, 0));
    meci_C::microa[1][1] = 2; meci_C::microa[2][1] = 0;
    meci_C::microa[1][2] = 1; meci_C::microa[2][2] = 1;
    meci_C::microb[1][1] = 0; meci_C::microb[2][1] = 0;
    meci_C::microb[1][2] = 1; meci_C::microb[2][2] = 1;
    // eigs 1-based: eigs[1]=-1.0, eigs[2]=-0.5
    double eigs[4] = {0, -1.0, -0.5, 0};
    // xy all 0.1 (nmos^4=16)
    std::vector<double> xy(16, 0.1);
    double eiga[4] = {0, 0, 0, 0};
    double diag[4] = {0, 0, 0, 0};
    double gse = 0.0;
    mecid(eigs, gse, eiga, diag, xy.data());
    // eiga(1) = -1.0 - (2*0.1-0.1)*2 = -1.2 ; eiga(2) = -0.5 - (2*0.1-0.1)*2 = -0.7
    CHK_D(eiga[1], -1.2, 1e-12, "eiga(1)");
    CHK_D(eiga[2], -0.7, 1e-12, "eiga(2)");
    // gse = eiga1*2*2 + xy1111*4 = -4.8 + 0.4 = -4.4
    CHK_D(gse, -4.4, 1e-12, "gse");
    // diag(1) closed shell = eiga1 - gse = -1.2 + 4.4 = 3.2
    CHK_D(diag[1], 3.2, 1e-12, "diag(1) closed shell");
    // diag(2) doubly-occ: alpha: eiga1 + 2*(0.1) + eiga2 + 2*(0.1) = -1.0-0.7+0.4... 
    //   = -1.2+0.2 -0.7+0.2 = -1.5 ; beta: eiga1 + eiga2 = -1.9 -> x=-3.4
    //   diag = -3.4 + 4.4 = 1.0
    CHK_D(diag[2], 1.0, 1e-12, "diag(2) doubly occupied");
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
