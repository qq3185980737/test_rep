// test_m10_mullik.cpp — M10 mullik (Mulliken population analysis).
// Acceptance: ALL PASS + ASan clean + zero warnings.
#define _CRT_SECURE_NO_WARNINGS
#include "mullik.h"

#include "maps_C.h"
#include "mod_vars_cuda.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "symmetry_C.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace parameters_C;
using namespace common_arrays_C;
using namespace chanel_C;
using namespace symmetry_C;
using namespace maps_C;
using namespace mod_vars_cuda;

static int fails = 0;
#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        if (cond) {                                                         \
            std::fprintf(stderr, "PASS: %s\n", msg);                        \
        } else {                                                            \
            std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__);    \
            ++fails;                                                        \
        }                                                                   \
    } while (0)
#define CHECK_NEAR(a, b, tol, msg)                                          \
    do {                                                                    \
        double _a = (a), _b = (b);                                          \
        if (std::fabs(_a - _b) <= (tol)) {                                  \
            std::fprintf(stderr, "PASS: %s\n", msg);                        \
        } else {                                                            \
            std::fprintf(stderr, "FAIL: %s got %g want %g (line %d)\n", msg, \
                        _a, _b, __LINE__);                                  \
            ++fails;                                                        \
        }                                                                   \
    } while (0)

// t1: H2 (2 orbitals) — non-graph path.
void t_mullik_h2() {
    numat = 2;
    norbs = 2;
    nelecs = 2;
    nclose = 1;
    nopen = 0;
    fract = 0.0;
    uhf = false;
    numcal = 1;
    keywrd = " MULLIK";
    mod_vars_cuda::lgpu = false;
    rxn_coord = 1.0e9;  // > 1e8 -> no reaction-coordinate write
    escf = -10.0;
    method_pm6 = false;
    verson = "19.0";
    id = 0;
    line = "";
    nalpha = 1;
    nbeta = 1;

    nat.assign(4, 0);
    nat[1] = 1;
    nat[2] = 1;
    nfirst.assign(4, 0);
    nlast.assign(4, 0);
    nfirst[1] = 1;
    nlast[1] = 1;
    nfirst[2] = 2;
    nlast[2] = 2;
    tore[1] = 1.0;
    betas[1] = -2.678;
    betap[1] = 0.0;
    zs[1] = 1.245;
    zp[1] = 0.0;
    zd[1] = 0.0;

    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1] = 0.0;
    coord[1][2] = 0.741;
    tvec.assign(4, std::vector<double>(4, 0.0));

    // packed h (overlap): h11=1, h21=0.8, h22=1
    h.assign(8, 0.0);
    h[1] = 1.0;
    h[2] = 0.8;
    h[3] = 1.0;
    std::vector<double> h_orig = h;

    // packed density: H2 bonding density
    p.assign(8, 0.0);
    p[1] = 0.7;
    p[2] = 0.3;
    p[3] = 0.7;

    // S-orthonormal eigenvectors: a = 1/sqrt(3.6), b = 1/sqrt(0.4)
    c.assign(4, std::vector<double>(4, 0.0));
    double a_ = 1.0 / std::sqrt(3.6);
    double b_ = 1.0 / std::sqrt(0.4);
    c[1][1] = a_;
    c[2][1] = a_;
    c[1][2] = b_;
    c[2][2] = -b_;
    q.assign(4, 0.0);
    eigs.assign(4, 0.0);
    eigb.assign(4, 0.0);
    namo.clear();
    jndex.assign(4, 0);
    jndex[1] = 1;
    jndex[2] = 2;

    ifact.clear();
    pb.clear();
    mullik();

    // q = tore - electron density (H2: each H ~ 1 e)
    std::fprintf(stderr, "t1 q1=%g q2=%g pb1=%g pb3=%g\\n", q[1], q[2], pb[ifact[2]], pb[ifact[3]]);
    CHECK(std::fabs(q[1] - 0.3) < 0.1, "mullik: q(H1) = 1 - 0.7 = 0.3");
    // h must be restored to the input overlap matrix
    CHECK_NEAR(h[1], h_orig[1], 1e-9, "mullik: h restored (diag)");
    CHECK_NEAR(h[2], h_orig[2], 1e-9, "mullik: h restored (off-diag)");
    // per-orbital populations on the diagonal of pb: ~1 e each
    CHECK(std::fabs(pb[ifact[2]] - pb[ifact[3]]) < 1e-6, "mullik: symmetric populations");
    CHECK(pb[ifact[2]] > 0.4, "mullik: population positive");
}

// t2: GRAPHF formatted output file.
void t_mullik_graphf() {
    numat = 1;
    norbs = 1;
    nelecs = 1;
    nclose = 1;
    nopen = 0;
    fract = 0.0;
    uhf = false;
    numcal = 2;
    keywrd = " GRAPHF MULLIK";
    mod_vars_cuda::lgpu = false;
    rxn_coord = 1.0e9;
    escf = -5.0;
    method_pm6 = false;
    verson = "19.0";
    id = 0;
    line = "";
    nalpha = 1;
    nbeta = 1;
    gpt_fn = "test_mullik_out.mgf";

    nat.assign(4, 0);
    nat[1] = 1;
    nfirst.assign(4, 0);
    nlast.assign(4, 0);
    nfirst[1] = 1;
    nlast[1] = 1;
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1] = 0.0;
    coord[2][1] = 0.0;
    coord[3][1] = 0.0;
    tvec.assign(4, std::vector<double>(4, 0.0));

    h.assign(4, 0.0);
    h[1] = 1.0;
    p.assign(4, 0.0);
    p[1] = 1.0;
    c.assign(4, std::vector<double>(4, 0.0));
    c[1][1] = 1.0;
    cb = c;
    q.assign(4, 0.0);
    eigs.assign(4, 0.0);
    eigb.assign(4, 0.0);
    namo.clear();
    jndex.assign(4, 0);
    jndex[1] = 1;
    ifact.clear();
    pb.clear();
    mullik();

    FILE* f = std::fopen("test_mullik_out.mgf", "r");
    bool ok = (f != nullptr);
    std::string head;
    if (ok) {
        char buf[512];
        if (std::fgets(buf, sizeof(buf), f)) head = buf;
        std::fclose(f);
    }
    CHECK(ok, "mullik: GRAPHF file created");
    CHECK(head.find("MOPAC-Graphical data Version 2012") != std::string::npos,
          "mullik: GRAPHF header line");
    std::remove("test_mullik_out.mgf");
}

int main() {
    std::fprintf(stderr, "== M10 mullik ==\n");
    std::fprintf(stderr, "[t1 H2]\n");
    t_mullik_h2();
    std::fprintf(stderr, "[t2 GRAPHF]\n");
    t_mullik_graphf();
    if (fails == 0) {
        std::fprintf(stderr, "ALL PASS: M10 mullik\n");
        return 0;
    }
    std::fprintf(stderr, "FAILURES: %d\n", fails);
    return 1;
}
