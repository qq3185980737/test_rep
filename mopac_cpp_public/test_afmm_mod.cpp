// test_afmm_mod.cpp — numeric checks for get_legendre (afmm_ini runs).

#include <cmath>
#include <cstdio>

#include "afmm_mod.h"
#include "chanel_C.h"

using namespace afmm_mod;

static bool approx(double a, double b, double tol = 1e-12) {
    return std::fabs(a - b) < tol;
}

int main() {
    afmm_ini();  // must run cleanly

    RArr pmn;

    // x=0: P0=1, P1=0, P2=-1/2.
    double x = 0.0;
    get_legendre(P, x, pmn);
    bool a2 = approx(pmn(0, 0), 1.0) && approx(pmn(0, 1), 0.0) && approx(pmn(0, 2), -0.5);
    std::printf("legendre x=0: P0=%g P1=%g P2=%g  %s\n", pmn(0,0), pmn(0,1), pmn(0,2), a2?"PASS":"FAIL");

    // x=1: Pn(1)=1.
    x = 1.0;
    get_legendre(P, x, pmn);
    bool a3 = approx(pmn(0,1),1.0) && approx(pmn(0,2),1.0) && approx(pmn(0,3),1.0,1e-9);
    std::printf("legendre x=1: P1=%g P2=%g P3=%g  %s\n", pmn(0,1), pmn(0,2), pmn(0,3), a3?"PASS":"FAIL");

    // x=-1: Pn(-1)=(-1)^n.
    x = -1.0;
    get_legendre(P, x, pmn);
    bool a4 = approx(pmn(0,1),-1.0) && approx(pmn(0,2),1.0) && approx(pmn(0,3),-1.0);
    std::printf("legendre x=-1: P1=%g P2=%g P3=%g  %s\n", pmn(0,1), pmn(0,2), pmn(0,3), a4?"PASS":"FAIL");

    // x=0.5: P2=0.5*(3*0.25-1)=-0.125.
    x = 0.5;
    get_legendre(P, x, pmn);
    bool a5 = approx(pmn(0,2), -0.125, 1e-12);
    std::printf("legendre x=0.5: P2=%g (expect -0.125)  %s\n", pmn(0,2), a5?"PASS":"FAIL");

    return (a2 && a3 && a4 && a5) ? 0 : 1;
}
