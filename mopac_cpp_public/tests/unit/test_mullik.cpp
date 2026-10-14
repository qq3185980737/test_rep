// test_mullik.cpp — end-to-end smoke: H2 (2 s-orbitals), S=[[1,.5],[.5,1]],
// P=I. Mulliken charges must be 0/0; per-orbital populations ~1.
#include "mullik.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "chanel_C.h"
#include "maps_C.h"
#include "symmetry_C.h"
#include "mod_vars_cuda.h"
#include <cstdio>
#include <cmath>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace chanel_C;
using namespace maps_C;
using namespace symmetry_C;

static bool near(double a, double b, double tol = 1e-8) {
    return std::fabs(a - b) <= tol;
}

int main() {
    bool ok = true;

    numat = 2; norbs = 2; numcal = 1; nelecs = 2; nclose = 1; nopen = 0;
    fract = 1.0; uhf = false; nalpha = 1; nbeta = 0; mozyme = false;
    mod_vars_cuda::lgpu = false;
    keywrd = " AM1 GRAPHF MULLIK";   // GRAPH formatted path + full population analysis
    verson = "2016";
    gpt_fn = "test_mullik.mgf";
    escf = -1.0; rxn_coord = 1.0e9;

    nat.assign(3, 0); nat[1] = 1; nat[2] = 1;
    nfirst.assign(3, 0); nlast.assign(3, 0);
    nfirst[1] = 1; nlast[1] = 1;
    nfirst[2] = 2; nlast[2] = 2;
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 0.7; coord[2][2] = 0.0; coord[3][2] = 0.0;
    tore[1] = 1.0; zs[1] = 1.3; zp[1] = 0.0; zd[1] = 0.0;
    betas[1] = -6.1737870; betap[1] = 0.0;

    // Packed h for S = [[1, 0.5],[0.5, 1]] (lower triangle, 1-based):
    // h(1,1)=1, h(2,1)=0.5, h(2,2)=1  -> h[1],h[2],h[3]
    h.assign(4, 0.0); h[1] = 1.0; h[2] = 0.5; h[3] = 1.0;
    // Density P = I (packed): p(1,1)=1, p(2,1)=0, p(2,2)=1
    p.assign(4, 0.0); p[1] = 1.0; p[2] = 0.0; p[3] = 1.0;
    q.assign(3, 0.0); pb.assign(4, 0.0);
    c.assign(3, std::vector<double>(3, 0.0));   // c[i][j], 1-based
    // Physical H2 MOs (AO coefficients): bonding c1=[.577,.577] (c1^T S c1 = 1
    // for S=[[1,.5],[.5,1]]), antibonding c2=[.707,-.707] (orthogonal to c1).
    c[1][1] = 0.5773502692; c[1][2] = 0.7071067812;
    c[2][1] = 0.5773502692; c[2][2] = -0.7071067812;
    cb.assign(3, std::vector<double>(3, 0.0));
    eigs.assign(3, 0.0); eigb.assign(3, 0.0);
    ifact.assign(4, 0);
    namo.assign(3, "H");
    jndex.assign(3, 0);

    mullik();

    if (!(near(q[1], 0.0) && near(q[2], 0.0))) {
        std::printf("FAIL charges: q = %g %g\n", q[1], q[2]);
        ok = false;
    }
    // Per-orbital Mulliken populations: the summation loop stores each row
    // sum into pb(ifact(i+1)) (the diagonal slot). Reference values from
    // numpy with the beta-rescaled S: row1=0.68406, row2=0.55093.
    double p11 = pb[ifact[2]];   // row 1 population
    double p22 = pb[ifact[3]];   // row 2 population
    if (!(std::fabs(p11 - 0.68406) < 0.05 && std::fabs(p22 - 0.55093) < 0.05)) {
        std::printf("FAIL populations: row1=%g row2=%g\n", p11, p22);
        ok = false;
    }
    // h must be restored to the input overlap matrix.
    if (!(near(h[1], 1.0) && near(h[2], 0.5) && near(h[3], 1.0))) {
        std::printf("FAIL h restore: %g %g %g\n", h[1], h[2], h[3]);
        ok = false;
    }
    // GRAPH file must exist.
    FILE* f = std::fopen("test_mullik.mgf", "r");
    if (!f) { std::printf("FAIL mgf missing\n"); ok = false; }
    else std::fclose(f);

    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
