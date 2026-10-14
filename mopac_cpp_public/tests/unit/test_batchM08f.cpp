// test_batchM08f.cpp  M08 mecih (packed CI matrix) + aabbcd/aabacd/aababc/babbcd/babbbc
#include <cstdio>
#include <cmath>
#include <vector>
#include "aabbcd.h"
#include "aababc.h"
#include "aabacd.h"
#include "babbcd.h"
#include "babbbc.h"
#include "meci_C.h"
#include "mecih.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08f batch - mecih + matrix elements\n");
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
    meci_C::nalmat.assign({0, 1, 1});
    meci_C::ispqr.assign(3, std::vector<int>(24, 0));
    meci_C::is = 2;
    std::vector<double> diag(4, 0.0), cimat(8, 0.0);
    diag[1] = 3.2; diag[2] = 1.0;
    std::vector<double> xy(16, 0.1);
    mecih(diag.data(), cimat.data(), meci_C::nmos, &meci_C::lab, xy.data());
    CHK_D(cimat[1], 3.2, 1e-12, "cimat(1,1)=diag1");
    CHK_D(cimat[2], -0.1, 1e-12, "cimat(2,1)=aabbcd");
    CHK_D(cimat[3], 1.0, 1e-12, "cimat(2,2)=diag2");
    CHECK(meci_C::ispqr[1][1] == 1, "ispqr(1,1)");
    CHECK(meci_C::ispqr[2][1] == 1, "ispqr(2,1)");

    // Direct matrix-element checks on legal nmos=4 patterns (xy all 0.1).
    int c1[5] = {0,1,1,0,0}, d1[5] = {0,0,0,0,0};
    int c2[5] = {0,0,0,1,1};
    std::vector<double> xy4(4*4*4*4, 0.1);
    // aabacd (4 alpha diffs): i=3,j=4,k=1,l=2; xy(3,1,4,2)-xy(3,2,4,1)=0; even -> 0
    CHK_D(aabacd(c1, d1, c2, d1, 4, xy4.data()), 0.0, 1e-12, "aabacd direct (even parity)");
    // babbcd (4 beta diffs): beta1={1,1,0,0}, beta2={0,0,1,1}; i=3,j=4,k=1,l=2; 0
    int e1[5] = {0,1,1,0,0}, e2[5] = {0,0,0,1,1};
    CHK_D(babbcd(d1, e1, d1, e2, 4, xy4.data()), 0.0, 1e-12, "babbcd direct (even parity)");
    // aababc (1 alpha diff): a1={1,1,0,0}, a2={1,0,1,0}; i=2,j=3; sum=-0.4 even -> -0.4
    int f1[5] = {0,1,1,0,0}, g1[5] = {0,0,0,0,0}, f2[5] = {0,1,0,1,0};
    meci_C::occa.assign({0.0, 2.0, 0.0, 0.0, 0.0});
    CHK_D(aababc(f1, g1, f2, 4, xy4.data()), -0.2, 1e-12, "aababc direct");
    // babbbc (1 beta diff): b1={1,1,0,0}, b2={1,0,1,0}; i=2,j=3; sum=-0.4 even -> -0.4
    int h1[5] = {0,1,1,0,0}, h2[5] = {0,1,0,1,0};
    CHK_D(babbbc(d1, h1, h2, 4, xy4.data()), -0.2, 1e-12, "babbbc direct");

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
