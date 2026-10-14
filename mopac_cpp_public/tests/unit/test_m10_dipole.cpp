// test_m10_dipole.cpp — M10 dipole group:
// dipole + dipole_for_MOZYME + ionout.
// Acceptance: ALL PASS + ASan clean + zero warnings.
#define _CRT_SECURE_NO_WARNINGS
#include "dipole.h"
#include "dipole_for_MOZYME.h"
#include "ionout.h"

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "to_screen_C.h"
#include "elemts_C.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <io.h>

using common_arrays_C::nat;
using common_arrays_C::p;
using common_arrays_C::coord;
using common_arrays_C::q;
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using common_arrays_C::atmass;
using common_arrays_C::na;
using common_arrays_C::nb;
using common_arrays_C::nc;
using common_arrays_C::txtatm;
using common_arrays_C::labels;
using molkst_C::numat;
using molkst_C::numcal;
using molkst_C::keywrd;
using molkst_C::mol_weight;
using molkst_C::natoms;
using molkst_C::maxtxt;
using molkst_C::mpack;
using molkst_C::moperr;
using parameters_C::dd;
using parameters_C::ddp;
using parameters_C::tore;
using parameters_C::ams;
using funcon_C::fpc_8;
using funcon_C::fpc_1;
using funcon_C::a0;
using elemts_C::elemnt;


// ijbo stub (host provides the real implementation in ijbo.cpp; excluded
// from this link so the test can fix deterministic offsets).
int ijbo(int, int) { return 0; }

static int fails = 0;
#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        if (cond) {                                                         \
            std::printf("PASS: %s\n", msg);                                 \
        } else {                                                            \
            std::printf("FAIL: %s (line %d)\n", msg, __LINE__);             \
            ++fails;                                                        \
        }                                                                   \
    } while (0)
#define CHECK_NEAR(a, b, tol, msg)                                          \
    do {                                                                    \
        double _a = (a), _b = (b);                                          \
        if (std::fabs(_a - _b) <= (tol)) {                                  \
            std::printf("PASS: %s\n", msg);                                 \
        } else {                                                            \
            std::printf("FAIL: %s got %g want %g (line %d)\n", msg, _a, _b, \
                        __LINE__);                                          \
            ++fails;                                                        \
        }                                                                   \
    } while (0)

// t1: single H, non-charged -> point-charge term only.
void t_dipole_h() {
    numat = 1;
    numcal = 1;
    keywrd = "";
    nat.assign(4, 0);
    nat[1] = 1;
    nfirst.assign(4, 0);
    nlast.assign(4, 0);
    nfirst[1] = 1;
    nlast[1] = 1;
    q.assign(4, 0.0);
    q[1] = 0.4;  // |sum| < 0.5 -> not chargd -> no COM shift
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1] = 1.0;  // H at (1,0,0)
    p.assign(10, 0.0);
    std::vector<double> dipvec(4, 0.0);
    double factdm = fpc_8 * fpc_1 * 1e-10;
    double d = dipole(p, coord, dipvec, 0);
    CHECK_NEAR(d, 0.4 * 1.0 * factdm, 1e-9, "dipole: H point-charge");
    CHECK(dipvec[1] == 0.0 && dipvec[2] == 0.0 && dipvec[3] == 0.0,
          "dipole: dipvec untouched when not FORCE");
}

// t2: single C with sp + pd terms.
void t_dipole_c() {
    numat = 1;
    numcal = 2;
    keywrd = "";
    nat[1] = 6;
    nfirst[1] = 1;
    nlast[1] = 9;  // sp + d -> l == 8 branch
    q.assign(4, 0.0);
    q[1] = 0.4;
    coord.assign(4, std::vector<double>(4, 0.0));  // C at origin
    p.assign(64, 0.0);
    dd[6] = 0.5;    // test value (normally from parameter set)
    ddp[5][6] = 0.3;
    p[2] = 0.5;  // sp-x (k = (ia+j)(ia+j-1)/2+ia, ia=1, j=1)
    p[4] = 0.5;  // sp-y
    p[7] = 0.5;  // sp-z
    p[19] = 1.0; // dx: (ia+5)(ia+4)/2+ia+3
    std::vector<double> dipvec(4, 0.0);
    dipole(p, coord, dipvec, 0);
    double hyfsp = 2.0 * dd[6] * a0 * fpc_8 * fpc_1 * 1e-10;
    double hyfpd = 2.0 * ddp[5][6] * a0 * fpc_8 * fpc_1 * 1e-10;
    CHECK_NEAR(to_screen_C::dip[0][1], -(hyfsp * 0.5 + hyfpd * 1.0), 1e-9,
               "dipole: hybrid x (sp + pd)");
    CHECK_NEAR(to_screen_C::dip[1][1], -(hyfsp * 0.5), 1e-9,
               "dipole: hybrid y");
    CHECK_NEAR(to_screen_C::dip[2][1], -(hyfsp * 0.5), 1e-9,
               "dipole: hybrid z");
}

// t3: charged dimer -> COM shift applied then restored.
void t_dipole_com() {
    numat = 2;
    numcal = 3;
    keywrd = "";  // not FORCE -> COM shift
    nat[1] = 1;
    nat[2] = 1;
    nfirst[1] = 1;
    nlast[1] = 1;
    nfirst[2] = 2;
    nlast[2] = 2;
    q.assign(4, 0.0);
    q[1] = 1.0;
    q[2] = 1.0;  // sum = 2 -> chargd
    atmass.assign(108, 0.0);
    atmass[1] = 1.0;
    atmass[2] = 1.0;
    mol_weight = 2.0;
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1] = 0.0;
    coord[1][2] = 1.0;
    p.assign(10, 0.0);
    std::vector<double> dipvec(4, 0.0);
    double d = dipole(p, coord, dipvec, 0);
    // symmetric about COM -> net dipole ~ 0
    CHECK_NEAR(d, 0.0, 1e-12, "dipole: charged dimer dipole ~ 0");
    CHECK_NEAR(coord[1][1], 0.0, 1e-12, "dipole: coord restored (H1)");
    CHECK_NEAR(coord[1][2], 1.0, 1e-12, "dipole: coord restored (H2)");
}

// t4: dipole_for_MOZYME, H atom, non-charged.
void t_dipole_mozyme() {
    numat = 1;
    numcal = 4;
    keywrd = "";
    nat.assign(4, 0);
    nat[1] = 1;
    MOZYME_C::iorbs.assign(4, 0);
    MOZYME_C::iorbs[1] = 1;  // H: 1 orbital -> l = 0 -> no hybrid term
    q.assign(4, 0.0);
    q[1] = 0.4;
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1] = 1.0;
    p.assign(10, 0.0);
    std::vector<double> dipvec(4, 0.0);
    double factdm = fpc_8 * fpc_1 * 1e-10;
    double d = dipole_for_MOZYME(dipvec, 0);
    CHECK_NEAR(d, 0.4 * 1.0 * factdm, 1e-9, "dipole_for_MOZYME: point charge");
}

// t5: ionout prints ion table (m==3 with connectivity + charge).
void t_ionout() {
    natoms = 3;
    numat = 3;
    maxtxt = 4;
    labels.assign(4, 0);
    labels[1] = 99;  // dummy atom
    labels[2] = 6;
    labels[3] = 1;
    na.assign(4, 0);
    nb.assign(4, 0);
    nc.assign(4, 0);
    na[2] = 2;  // triggers Int + klim=3
    txtatm.assign(4, std::string(""));
    txtatm[2] = "C1";
    txtatm[3] = "H1";
    // ions(4,2) col-major: ion1 atom=2, ion2 atom=3, charges +1/+2
    int ions[8] = {2, 3, 0, 0, 1, 2, 0, 0};
    std::fflush(stdout);
    FILE* cap = std::freopen("_ionout_cap.txt", "w", stdout);
    if (cap) {
        ionout(ions, 3, 2);
        std::fflush(stdout);
        std::fclose(stdout);
    }
    std::string out;
    FILE* f = std::fopen("_ionout_cap.txt", "r");
    if (f) {
        char buf[512];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
    }
}

int main() {
    std::fprintf(stderr, "== M10 dipole group ==\n");
    std::fprintf(stderr, "[t1 dipole H]\n");
    t_dipole_h();
    std::fprintf(stderr, "[t2 dipole C]\n");
    t_dipole_c();
    std::fprintf(stderr, "[t3 dipole COM]\n");
    t_dipole_com();
    std::fprintf(stderr, "[t4 dipole_for_MOZYME]\n");
    t_dipole_mozyme();
    std::fprintf(stderr, "[t5 ionout]\n");
    t_ionout();
    if (fails == 0) {
        std::fprintf(stderr, "ALL PASS: M10 dipole group\n");
        return 0;
    }
    std::fprintf(stderr, "FAILURES: %d\n", fails);
    return 1;
}
