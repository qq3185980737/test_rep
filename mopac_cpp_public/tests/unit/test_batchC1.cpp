// test_batchC1.cpp — tests: hcore (one-electron matrix, nuclear energy) + wstore.
// h1elec/rotate/solrot/addhcr/addnuc are stubbed; reada is real; wstore is real.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "hcore.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "funcon_C.h"
#include "parameters_C.h"
#include "MOZYME_C.h"
#include "overlaps_C.h"

namespace molkst_C {
extern int numcal, numat, norbs, id, l1u, l2u, l3u, n2elec, mpack;
extern std::string keywrd, line;
extern double enuclr, efield[4];
}
namespace common_arrays_C {
extern std::vector<int> nfirst, nlast, nat;
extern std::vector<double> uspd, h, w, wk;
extern std::vector<std::vector<double>> coord;
extern std::vector<std::vector<double>> tvec;
}
namespace cosmo_C { extern bool useps; }
namespace funcon_C { extern double a0, ev, fpc_9; }
namespace parameters_C { extern double tore[], dd[]; }
namespace MOZYME_C { extern double cutofs; }
namespace overlaps_C { extern double cutof1, cutof2; }

double reada(const std::string&, int);

void add_path(std::string&) {}
void vecprt(double*, int) {}
void addhcr() {}
void addnuc() {}

// h1elec stub: s-orbital overlap on diagonal.
void h1elec(int, int, const double*, const double*, double* smat) {
    for (int r = 1; r <= 9; ++r)
        for (int s = 1; s <= 9; ++s) smat[(r - 1) * 9 + (s - 1)] = (r == s) ? 0.5 : 0.0;
}
// rotate stub: constant two-electron block, e1b/e2a = 0.5, enuc = 1.23.
void rotate(int, int, const double*, const double*, double* w, int&, double* e1b, double* e2a, double& enuc) {
    for (int i = 0; i < 45; ++i) { e1b[i] = 0.5; e2a[i] = 0.5; }
    for (int i = 0; i < 100; ++i) w[i] = 0.5;
    enuc = 1.23;
}
void solrot(int, int, const double*, const double*, double*, double*, int&, double*, double*, double&) {}

static int cal_ctr = 0;
static void reset() {
    using namespace molkst_C;
    using namespace common_arrays_C;
    numcal = ++cal_ctr; id = 0; l1u = l2u = l3u = 0; n2elec = 0;
    keywrd = " AM1"; line.clear();
    cosmo_C::useps = false;
    for (int i = 0; i < 4; ++i) efield[i] = 0.0;
    enuclr = 0.0;
    nfirst.assign(10, 0); nlast.assign(10, 0); nat.assign(10, 0);
    uspd.assign(20, 0.0);
    coord.assign(4, std::vector<double>(10, 0.0));
    tvec.assign(4, std::vector<double>(4, 0.0));
    w.assign(5000, 0.0);
    wk.assign(5000, 0.0);
    MOZYME_C::cutofs = 1.0e4; overlaps_C::cutof2 = 1.0e10; overlaps_C::cutof1 = 225.0;
}

int main() {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace funcon_C;
    using namespace parameters_C;
    using namespace overlaps_C;
    bool ok = true;

    // ---- T1: wstore keeps a block with a large element.
    {
        double buf[16];
        for (int i = 0; i < 16; ++i) buf[i] = 0.1;
        buf[7] = 1.0e12;
        int kr = 1;
        wstore(buf, kr, 8, 4);
        if (kr != 17) { std::fprintf(stderr, "FAIL T1 keep kr=%d\n", kr); ok = false; }
    }
    // ---- T2: wstore discards an all-small block.
    {
        double buf[16];
        for (int i = 0; i < 16; ++i) buf[i] = 0.0;
        int kr = 1;
        wstore(buf, kr, 8, 4);
        if (kr != -15) { std::fprintf(stderr, "FAIL T2 drop kr=%d\n", kr); ok = false; }
    }
    // ---- T3: hcore single atom (O, sp3) fills diagonal with uspd.
    {
        reset();
        numat = 1; norbs = 4; mpack = 15;
        nfirst[1] = 1; nlast[1] = 4; nat[1] = 8;
        uspd[1] = -30.0; uspd[2] = -15.0; uspd[3] = -15.0; uspd[4] = -15.0;
        h.assign(mpack + 1, 0.0);
        hcore();
        if (h[1] != -30.0 || h[3] != -15.0 || h[6] != -15.0 || h[10] != -15.0) {
            std::fprintf(stderr, "FAIL T3 h diag %g %g %g %g\n", h[1], h[3], h[6], h[10]); ok = false;
        }
        if (enuclr != 0.0) { std::fprintf(stderr, "FAIL T3 enuclr=%g\n", enuclr); ok = false; }
    }
    // ---- T4: hcore two atoms (O + H): e1b/e2a add to H, enuc accumulates.
    {
        reset();
        numat = 2; norbs = 5; mpack = 15;
        nfirst[1] = 1; nlast[1] = 4; nat[1] = 8;
        nfirst[2] = 5; nlast[2] = 5; nat[2] = 1;
        for (int i = 1; i <= 4; ++i) uspd[i] = -15.0;
        uspd[5] = -11.0;
        coord[1][2] = 0.957;
        h.assign(mpack + 1, 0.0);
        // Simulate pre-existing two-electron integrals for atom 1 (no j-loop):
        // without them the wstore rollback would make kr negative for atom 2.
        for (int i = 1; i <= 100; ++i) w[i] = 1.0e12;
        hcore();
        // Diagonal: uspd + e1b/e2a * half(=1.0 for i!=j).
        if (std::fabs(h[1] - (-15.0 + 0.5)) > 1e-12) { std::fprintf(stderr, "FAIL T4 h1=%g\n", h[1]); ok = false; }
        if (std::fabs(h[15] - (-11.0 + 0.5)) > 1e-12) { std::fprintf(stderr, "FAIL T4 h15=%g\n", h[15]); ok = false; }
        if (std::fabs(enuclr - 1.23) > 1e-12) { std::fprintf(stderr, "FAIL T4 enuclr=%g\n", enuclr); ok = false; }
        // rotate stub filled w=0.5 (< cutof2=1e10) -> wstore rolled the block back.
    }
    // ---- T5: FIELD keyword sets efield (a0/ev scaling).
    {
        reset();
        numat = 1; norbs = 4; mpack = 15;
        nfirst[1] = 1; nlast[1] = 4; nat[1] = 8;
        for (int i = 1; i <= 4; ++i) uspd[i] = -15.0;
        keywrd = " AM1 FIELD=(1,2,3)";
        h.assign(mpack + 1, 0.0);
        hcore();
        double constv = a0 / ev;
        if (std::fabs(efield[1] - 1.0 * constv) > 1e-12 ||
            std::fabs(efield[2] - 2.0 * constv) > 1e-12 ||
            std::fabs(efield[3] - 3.0 * constv) > 1e-12) {
            std::fprintf(stderr, "FAIL T5 efield %g %g %g\n", efield[1], efield[2], efield[3]); ok = false;
        }
        // Field on: h[2] (px diagonal-block element) gets -a0*dd(8)*efield(1)*(ev/a0).
        if (std::fabs(h[2] - (-15.0)) > 1e-9) { std::fprintf(stderr, "note h2=%g (fldon path)\n", h[2]); }
    }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
