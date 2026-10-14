// sym_test.cpp — check diat_ symmetry: diat(1,6,xjuc) vs diat(6,1,xjuc)
// same xjuc must give transposed smat with equal (0,0) element
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include "molkst_C.h"
#include "parameters_C.h"
#include "overlaps_C.h"
#include "MOZYME_C.h"
extern "C" void diat_(int* ni, int* nj, double* xjuc, double* smat);
extern void switch_method();
extern "C" void vdExp(const int* n, const double* a, double* y) {
    for (int i = 0; i < *n; ++i) y[i] = std::exp(a[i]);
}

static void hexdump(const char* tag, const double* smat, int n) {
    printf("%s:", tag);
    for (int i = 0; i < n; ++i) {
        unsigned long long b; std::memcpy(&b, &smat[i], 8);
        printf(" %016llx", (unsigned long long)b);
    }
    printf("\n");
}

int main() {
    using namespace parameters_C;
    // PM7-like Slater exponents: H (1s), C (2s/2p)
    natorb[1] = 1; zs[1] = 1.2458568; zp[1] = 0.0; zd[1] = 0.0;
    natorb[6] = 4; zs[6] = 1.6412075; zp[6] = 1.8;  zd[6] = 0.0;
    npq[1][1] = 1; npq[6][1] = 2;
    overlaps_C::cutof1 = 1e4;  // disable distance cutoff
    double xjuc[3] = {1.1, 0.0, 0.0};  // H-C distance ~1.1 A
    double s1[81] = {}, s2[81] = {};
    int ni = 1, nj = 6;
    diat_(&ni, &nj, xjuc, s1);   // H s x C
    ni = 6; nj = 1;
    diat_(&ni, &nj, xjuc, s2);   // C s x H
    printf("s1[0]  = %.17g\n", s1[0]);
    printf("s2[0]  = %.17g\n", s2[0]);
    printf("equal (0,0): %s\n", s1[0] == s2[0] ? "YES" : "NO");
    hexdump("s1", s1, 12);
    hexdump("s2", s2, 12);

    // --- dhc call 46 reproduction: H-H pair, ni=nj=1, PM7 params ---
    // xi from dhc probe call 46 (FTN & CXX identical):
    //   A = (5e-5, 5e-5, 1.00005), B = (1.5401, 0, 0.9999999999999999)
    double xa[3] = {5e-5, 5e-5, 1.00005};
    double xb[3] = {1.5401, 0.0, 0.9999999999999999};
    double sm[81] = {};
    double s_pre[81] = {};
    // PM7 Slater exponents: zs7(1)=1.260237 (overwrite the MNDO-style 1.2458568 above)
    zs[1] = 1.260237; zp[1] = 0.0; zd[1] = 0.0;
    // raw diat (scaled later by h1elec; here print pre-scaling by calling diat_)
    int ni1 = 1, nj1 = 1;
    double xjuc2[3] = {xb[0]-xa[0], xb[1]-xa[1], xb[2]-xa[2]};
    diat_(&ni1, &nj1, xjuc2, s_pre);
    printf("\n--- dhc call 46 repro (H-H, PM7 zs=1.260237) ---\n");
    printf("pre-scale s_pre[0] = %.17g\n", s_pre[0]);
    hexdump("s_pre", s_pre, 3);
    // now call the actual h1elec (scaled)
    extern void h1elec(int, int, const double*, const double*, double*);
    MOZYME_C::cutofs = 1.0e4;
    // test betas[1] candidates: PM7 -8.389745
    for (double btest : {-8.389745}) {
        betas[1] = btest; betap[1] = 0.0; betad[1] = 0.0;
        h1elec(1, 1, xa, xb, sm);
        unsigned long long b; std::memcpy(&b, &sm[0], 8);
        printf("betas[1]=%.17g -> sm[0]=%.17g hex %016llx  %s\n", btest, sm[0], b,
               b == 0xBFFF5EF3AC24E509ULL ? "MATCH FTN" : (b == 0xBFFF5EF3AC24E508ULL ? "1ulp below (CXX main)" : "other"));
    }
    return 0;
}
