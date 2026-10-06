// test_batchM03a.cpp — M03 SCF 核心（密度/能量/电荷/收敛判据）语义断言
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "densit.h"
#include "helect.h"
#include "chrge.h"
#include "scfcri.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_densit() {
    std::printf("Test 1: densit density matrix\n");
    // C = 4x4 identity (c[i][m] = delta): norbs=4, nocc=2, occ=1, nfract=0
    int mdim = 5;
    std::vector<std::vector<double>> c(mdim, std::vector<double>(mdim, 0.0));
    for (int i = 1; i <= 4; ++i) c[i][i] = 1.0;
    std::vector<double> p(11, -1.0);
    densit(c, mdim, 4, 2, 1.0, 0, 0.0, p, 1);
    // P(l) = sum_{m<=nocc} c[i][m]*c[j][m]  =>  P11=1, P22=1, rest 0
    CHECK(std::fabs(p[1] - 1.0) < 1e-12 && std::fabs(p[3] - 1.0) < 1e-12,
          "densit occupied diagonals");
    bool rest = true;
    for (int l = 2; l <= 10; ++l) if (l != 3 && std::fabs(p[l]) > 1e-12) rest = false;
    CHECK(rest, "densit off-diagonal/empty zero");
    // fractional: nocc=1 occ=1 nfract=2 fract=0.5, C columns: col1=(1,0,0,0), col2=(0,1,0,0)
    std::vector<std::vector<double>> c2(mdim, std::vector<double>(mdim, 0.0));
    c2[1][1] = 1.0; c2[2][2] = 1.0;
    std::vector<double> p2(11, 0.0);
    densit(c2, mdim, 4, 1, 1.0, 2, 0.5, p2, 1);
    // P11 = occ*1 + fract*0 = 1; P22 = occ*0 + fract*1 = 0.5
    CHK_D(p2[1], 1.0, 1e-12, "densit fractional P11");
    CHK_D(p2[3], 0.5, 1e-12, "densit fractional P22");
}

static void test_helect() {
    std::printf("Test 2: helect electronic energy\n");
    // n=3 packed lower 1-based p/h/f (size 6)
    double p[7] = {0, 1, 0.5, 1, 0.2, 0.3, 1};
    double h[7] = {0, 1, 2, 3, 1, 2, 4};
    double f[7] = {0, 0, 1, 1, 0, 0, 1};
    // F90: i=2: k=1 ed+=p1(h1+f1)=1; ee+=p2(h2+f2)=1.5; k=2
    //      i=3: k=3 ed+=p3(h3+f3)=4; ee+=p4(h4+f4)+p5(h5+f5)=0.8; k=2+3=5
    //      i=4: k=6 ed+=p6(h6+f6)=5 (ed=10)
    //      ee = 1.5+0.8 + 0.5*10 = 7.3
    double got = helect(3, p, h, f);
    CHK_D(got, 7.3, 1e-12, "helect energy 7.3");
}

static void test_chrge() {
    std::printf("Test 3: chrge atomic charges\n");
    using namespace molkst_C;
    using namespace common_arrays_C;
    molkst_C::numat = 2;
    molkst_C::mozyme = false;
    common_arrays_C::nfirst.resize(3); common_arrays_C::nlast.resize(3);
    common_arrays_C::nfirst[1] = 1; common_arrays_C::nlast[1] = 2;
    common_arrays_C::nfirst[2] = 3; common_arrays_C::nlast[2] = 4;
    // packed 4x4: p[l] 1-based
    std::vector<double> p(11, 0.0);
    p[1] = 2.0; p[3] = 1.0; p[6] = 4.0; p[10] = 0.5;
    std::vector<double> q(3, -1.0);
    chrge(p, q);
    // atom1: p[1]+p[3]=3; atom2: p[6]+p[10]=4.5
    CHK_D(q[1], 3.0, 1e-12, "chrge atom1");
    CHK_D(q[2], 4.5, 1e-12, "chrge atom2");
}

static void test_scfcri() {
    std::printf("Test 4: scfcri convergence criterion\n");
    using namespace molkst_C;
    double selcon = 0.0;
    molkst_C::numcal = 1;
    molkst_C::keywrd = " ";
    for (int i = 0; i < 4; ++i) molkst_C::efield[i] = 0.0;
    scfcri(selcon);
    CHK_D(selcon, 1e-2, 1e-15, "scfcri default 1e-2");
    molkst_C::numcal = 2;
    molkst_C::keywrd = "  PRECIS";
    scfcri(selcon);
    CHK_D(selcon, 1e-4, 1e-16, "scfcri PRECIS -> 1e-4");
    molkst_C::numcal = 3;
    molkst_C::keywrd = "  SCFCRT=1.D-6";
    scfcri(selcon);
    CHK_D(selcon, 1e-6, 1e-15, "scfcri SCFCRT=1.D-6");
    molkst_C::numcal = 4;
    molkst_C::keywrd = "  TS";
    scfcri(selcon);
    CHK_D(selcon, 1e-3, 1e-15, "scfcri TS -> 1e-3");
    // same numcal -> else branch: selcon = scfref (set by SCFCRT=1.D-6 -> 1e-6)
    molkst_C::keywrd = " ";
    scfcri(selcon);
    CHK_D(selcon, 1e-6, 1e-15, "scfcri repeated call keeps scfref");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M03a batch — SCF core: density / energy / charge / criterion\n");
    test_densit();
    test_helect();
    test_chrge();
    test_scfcri();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
