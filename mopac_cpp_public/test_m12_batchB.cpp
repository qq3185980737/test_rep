// test_m12_batchB.cpp — M12 batch B: sort / coe / formxy / genun / collid / outer1(mode0+mode1) / nxtmer
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
#include "sort.h"
#include "coe.h"
#include "formxy.h"
#include "genun.h"
#include "outer1.h"
#include "nxtmer.h"
#include "molkst_C.h"
#include "funcon_C.h"
#include "parameters_C.h"
#include "common_arrays_C.h"

using namespace molkst_C;
using namespace funcon_C;
using namespace parameters_C;
using namespace common_arrays_C;

static int passed = 0, failed = 0;
static void check(const char* name, bool ok) {
    if (ok) { ++passed; std::printf("  [PASS] %s\n", name); }
    else { ++failed; std::printf("  [FAIL] %s\n", name); }
}
static bool close(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) < tol;
}

int main() {
    std::printf("M12 batch B tests\n");

    // 1) sort: ascending selection sort, complex columns move with values
    {
        float val[4] = {0.0f, 3.0f, 1.0f, 2.0f};  // 1-based payload
        std::complex<float> vec[9];
        // 3x3 column-major: col1={1,0,0} col2={0,1,0} col3={0,0,1}
        for (int i = 0; i < 9; ++i) vec[i] = std::complex<float>(0, 0);
        vec[0] = 1; vec[4] = 1; vec[8] = 1;
        sort(val, vec, 3);
        bool ok = close(val[1], 1.0f) && close(val[2], 2.0f) && close(val[3], 3.0f);
        // col order follows: col1={0,1,0} col2={0,0,1} col3={1,0,0}
        ok = ok && close(vec[0].real(), 0.0) && close(vec[1].real(), 1.0) && close(vec[2].real(), 0.0);
        ok = ok && close(vec[3].real(), 0.0) && close(vec[4].real(), 0.0) && close(vec[5].real(), 1.0);
        ok = ok && close(vec[6].real(), 1.0) && close(vec[7].real(), 0.0) && close(vec[8].real(), 0.0);
        check("sort ascending + column swap", ok);
    }

    // 2) coe: x2=1,y2=0,z2=0 -> ca=1,cb=0,sa=0,sb=1 (hand-checked 0-based)
    {
        double c[75];
        double r;
        coe(1.0, 0.0, 0.0, 5, 5, c, r);
        bool ok = close(r, 1.0);
        ok = ok && close(c[36], 1.0);   // c(37)
        ok = ok && close(c[40], 1.0);   // c(41)=ca*sb
        ok = ok && close(c[52], -1.0);  // c(53)=-sb
        ok = ok && close(c[19], 1.0);   // c(20)=ca
        ok = ok && close(c[74], 0.5);   // c(75)=0.5*c2a*sb^2 (c2a=1)
        ok = ok && close(c[44], 0.86602540378444);  // rt34
        ok = ok && close(c[56], -1.0);  // c(57)=ca*c2b (c2b=-1)
        ok = ok && close(c[68], 0.86602540378444);  // rt13*1.5
        ok = ok && close(c[38], -0.5);  // c(39)=cb^2-0.5*sb^2
        ok = ok && close(c[5], -1.0);   // c(6)=-ca*sb
        ok = ok && close(c[17], 1.0);   // c(18)=c2a*sb
        check("coe rotation matrix", ok);
    }

    // 3) formxy: na=nb=1, w={pad,2.5}, ca=cb={pad,1.0}
    {
        std::vector<double> w(2, 0.0), wca(46, 0.0), wcb(46, 0.0), ca(46, 0.0), cb(46, 0.0);
        w[1] = 2.5; ca[1] = 1.0; cb[1] = 1.0;
        int kr = 0;
        formxy(w, kr, wca, wcb, ca, cb, 1, 1);
        // sum = cb[1]*w[1]*0.5 = 1.25 ; wca[1] += 1.25*0.5 = 0.625 (same for wcb)
        bool ok = close(wca[1], 0.625) && close(wcb[1], 0.625) && (kr == 1);
        check("formxy 1x1 electrostatics", ok);
    }

    // 4) genun: sphere points; n=10 -> nequat=5,nvert=2, nu ends at 7
    {
        std::vector<std::vector<double>> u(4, std::vector<double>(64, 0.0));
        int n = 10;
        genun(u, n);
        bool ok = (n == 7);
        ok = ok && close(u[1][1], 0.0) && close(u[3][1], 1.0);
        ok = ok && close(u[1][2], 1.0) && close(u[2][2], 0.0) && close(u[3][2], 0.0);
        ok = ok && close(u[3][7], -1.0);
        check("genun unit vectors", ok);
    }

    // 5) collid: probe rw=1.0
    {
        std::vector<double> cw(4, 0.0), rnbr(201, 0.0);
        std::vector<std::vector<double>> cnbr(4, std::vector<double>(201, 0.0));
        cw[1] = 0; cw[2] = 0; cw[3] = 0;
        rnbr[1] = 1.0;
        cnbr[1][1] = 3.0;   // far -> no collision
        bool ok = !collid(1.0, cw, cnbr, rnbr, 1, 0);
        cnbr[1][1] = 1.5;   // 1.5<2, dd2=2.25<4 -> collision
        ok = ok && collid(1.0, cw, cnbr, rnbr, 1, 0);
        ok = ok && !collid(1.0, cw, cnbr, rnbr, 1, 3);  // ishape=3 skips
        check("collid collision check", ok);
    }

    // 6) outer1 mode=0 (molecular): rij=a0 -> r=1, method_pm7 off
    {
        am[1] = 100.0; tore[1] = 1.0; natorb[1] = 1;
        method_pm7 = false;
        double w[2026] = {0.0}, e1b[45] = {0.0}, e2a[45] = {0.0};
        double c1[3] = {0, 0, 0}, c2[3] = {0, 0, a0};
        int kr = 0;
        double enuc = 0;
        outer1(1, 1, c1, c2, w, kr, e1b, e2a, enuc, 0, false);
        // aee = (0.005+0.005)^2 = 1e-4 ; ww = ev/sqrt(1.0001)
        double ww = ev / std::sqrt(1.0 + 1e-4);
        bool ok = close(w[0], ww, 1e-9);
        ok = ok && close(e1b[0], -ww) && close(e1b[2], -ww) && close(e1b[9], -ww) && close(e1b[44], -ww);
        ok = ok && close(e2a[0], -ww) && close(enuc, ww);
        ok = ok && (kr == 1);
        check("outer1 mode=0 molecular", ok);
    }

    // 7) outer1 mode=1 (solid-state): single cell, r=1A below clower -> trunk=r
    {
        l1u = l2u = l3u = 0;
        cutofp = 10.0; numcal = 1; clower = cupper = 0;
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) tvec[i].resize(4, 0.0);
        // tvec unit cell vectors
        tvec[1][1] = 5.0; tvec[2][2] = 5.0; tvec[3][3] = 5.0;
        am[1] = 100.0; tore[1] = 1.0; natorb[1] = 1;
        double w[2026] = {0.0}, e1b[45] = {0.0}, e2a[45] = {0.0};
        double c1[3] = {0, 0, 0}, c2[3] = {1.0, 0, 0};
        int kr = 0;
        double enuc = 0;
        outer1(1, 1, c1, c2, w, kr, e1b, e2a, enuc, 1, false);
        double ww = a0 * ev / 1.0;   // trunk(1)=1
        bool ok = close(w[0], ww, 1e-9) && close(enuc, ww);
        ok = ok && close(e1b[0], -ww) && close(e1b[25], -ww);
        ok = ok && (kr == 1);
        // far: r=20 > cupper=10 -> trunk=clim
        double c2b[3] = {20.0, 0, 0};
        double w2[2026] = {0.0}, e1b2[45] = {0.0}, e2a2[45] = {0.0};
        kr = 0; enuc = 0;
        outer1(1, 1, c1, c2b, w2, kr, e1b2, e2a2, enuc, 1, false);
        // clim = c+cr*cupper+cr2*cupper^2 (same numcal -> constants cached)
        double bound1 = clower / cutofp, range = 1.0 - bound1;
        double cT = -0.5 * bound1 * bound1 * cutofp / range;
        double crT = 1.0 + bound1 / range;
        double cr2T = -1.0 / (cutofp * 2.0 * range);
        double clim = cT + crT * cupper + cr2T * cupper * cupper;
        ok = ok && close(w2[0], a0 * ev / clim, 1e-9);
        check("outer1 mode=1 solid-state trunk", ok);
    }

    // 8) nxtmer: backbone detection N-C-C(=O)-OH + N-H
    {
        // atoms: 1=N, 2=C, 3=C, 4=O(terminal), 5=N, 6=H
        nat.resize(7); nbonds.resize(7);
        ibonds.assign(16, std::vector<int>(7, 0));
        nat[1] = 7; nat[2] = 6; nat[3] = 6; nat[4] = 8; nat[5] = 7; nat[6] = 1;
        nbonds[1] = 1; nbonds[2] = 2; nbonds[3] = 3; nbonds[4] = 1; nbonds[5] = 2; nbonds[6] = 1;
        ibonds[1][1] = 2;
        ibonds[1][2] = 1; ibonds[2][2] = 3;
        ibonds[1][3] = 2; ibonds[2][3] = 4; ibonds[3][3] = 5;
        ibonds[1][4] = 3;
        ibonds[1][5] = 3; ibonds[2][5] = 6;
        ibonds[1][6] = 5;
        int nbackb[4] = {0, 0, 0, 0};
        nxtmer(1, nbackb);
        bool ok = (nbackb[0] == 2 && nbackb[1] == 3 && nbackb[2] == 4 && nbackb[3] == 5);
        check("nxtmer backbone N-C-C(=O)-NH", ok);
    }

    std::printf("M12 batch B: %d passed, %d failed\n", passed, failed);
    return failed ? 1 : 0;
}
