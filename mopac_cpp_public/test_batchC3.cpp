// test_batchC3.cpp — tests: grid (2-D energy grid over two internal coords).
// ef/flepo/second/geout/pdbout/wrttxt/bonds/to_screen stubbed; reada real.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include "grid.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "maps_C.h"

namespace chanel_C {
extern int iw, iw0, ires, iarc, iump;
extern std::string archive_fn, ump_fn, restart_fn;
}
namespace common_arrays_C {
extern std::vector<std::vector<double>> geo, geoa;
extern std::vector<double> xparam, pa, pb, p;
extern std::vector<int> na, nb, nc;
}
namespace molkst_C {
extern int nvar, natoms, norbs, numat, mpack;
extern std::string keywrd, line;
extern double tleft;
extern bool gui, uhf, moperr;
}
namespace maps_C {
extern double rxn_coord1, rxn_coord2;
extern int ione, ijlp, ilp, jlp, jlp1;
extern std::vector<double> surf;
extern int lpara1, latom1, lpara2, latom2;
}

double reada(const std::string&, int);

double second(int) { return 0.0; }
void geout(int) {}
void pdbout(int) {}
void wrttxt(int) {}
void bonds() {}
void to_screen(const std::string&) {}
// flepo stub: E = first-var value + second-var value (both in geo).
void flepo(std::vector<double>&, double& funct1) {
    funct1 = common_arrays_C::geo[maps_C::lpara1][maps_C::latom1] +
             common_arrays_C::geo[maps_C::lpara2][maps_C::latom2];
}
void ef(std::vector<double>&, double& funct1) {
    funct1 = common_arrays_C::geo[maps_C::lpara1][maps_C::latom1] +
             common_arrays_C::geo[maps_C::lpara2][maps_C::latom2];
}

using namespace chanel_C;
using namespace common_arrays_C;
using namespace molkst_C;
using namespace maps_C;

static void setup() {
    nvar = 1; natoms = 2; norbs = 2; numat = 2; mpack = 3;
    gui = false; uhf = false; moperr = false;
    tleft = 1.0e6;
    keywrd = " GRID STEP1=0.1 STEP2=0.2 POINT1=2 POINT2=2";
    lpara1 = 1; latom1 = 1; lpara2 = 1; latom2 = 2;
    geo.assign(4, std::vector<double>(3, 0.0));
    geo[1][1] = 1.0; geo[1][2] = 2.0; geo[2][1] = 60.0;
    xparam.assign(10, 0.0);
    na.assign(3, 0); nb.assign(3, 0); nc.assign(3, 0);
    pa.assign(10, 0.0); pb.assign(10, 0.0); p.assign(10, 0.0);
    iw = 6; iw0 = 6; ires = 0; iarc = 0; iump = 0;
    rxn_coord1 = 0.0; rxn_coord2 = 0.0;
    ione = 0; ijlp = 0; ilp = 0; jlp = 0; jlp1 = 0;
}

// Read the third column of the ump surface file into vals (row-major i then j).
static int read_ump(const char* fn, double* vals, int maxn) {
    std::ifstream f(fn);
    if (!f.good()) return -1;
    int n = 0;
    std::string line;
    while (f.good() && std::getline(f, line)) {
        if (line.empty()) continue;
        double a, b, c;
        if (std::sscanf(line.c_str(), "%lf %lf %lf", &a, &b, &c) == 3 && n < maxn) vals[n++] = c;
    }
    return n;
}

int main() {
    bool ok = true;

    // ---- T1: 2x2 grid energies (flepo stub) recorded, order per snake scan.
    {
        setup();
        ump_fn = "C3_t1.ump";
        grid();
        // loop1 (0,0): 3.0; loop2 (0,1): 3.2; loop3 (1,1): 3.3; loop4 (1,0): 3.1
        // ump prints rows i=1..npts1, cols j=1..npts2: (1,1)=3.0 (1,2)=3.2 (2,1)=3.1 (2,2)=3.3
        double v[8];
        int n = read_ump("C3_t1.ump", v, 8);
        if (n != 4) { std::fprintf(stderr, "FAIL T1 ump lines=%d\n", n); ok = false; }
        else if (std::fabs(v[0] - 3.0) > 1e-9 || std::fabs(v[1] - 3.2) > 1e-9 ||
                 std::fabs(v[2] - 3.1) > 1e-9 || std::fabs(v[3] - 3.3) > 1e-9) {
            std::fprintf(stderr, "FAIL T1 ump %g %g %g %g\n", v[0], v[1], v[2], v[3]); ok = false;
        }
        if (iw0 != 6) { std::fprintf(stderr, "FAIL T1 iw0 restored=%d\n", iw0); ok = false; }
    }
    // ---- T2: angle coordinate: step1 scaled by 1/degree.
    {
        setup();
        ump_fn = "C3_t2.ump";
        lpara1 = 2;  // angle -> step1 /= degree
        geo[2][1] = 60.0;
        grid();
        // (iloop=2,jloop=2): geo[2][1]=60+1*step1'; step1' = 0.1/57.2958
        double step1a = 0.1 / 57.29577951308232;
        double exp3 = (60.0 + step1a) + (2.0 + 0.2);
        double v[8];
        int n = read_ump("C3_t2.ump", v, 8);
        if (n != 4 || std::fabs(v[3] - exp3) > 1e-3) { std::fprintf(stderr, "FAIL T2 v3=%g exp=%g n=%d\n", n == 4 ? v[3] : -1, exp3, n); ok = false; }
    }
    // ---- T3: SMOOTH: 4 passes, later passes cannot lower the stored value.
    {
        setup();
        ump_fn = "C3_t3.ump";
        keywrd += " SMOOTH";
        grid();
        double v[8];
        int n = read_ump("C3_t3.ump", v, 8);
        if (n != 4 || std::fabs(v[0] - 3.0) > 1e-9 || std::fabs(v[3] - 3.3) > 1e-9) {
            std::fprintf(stderr, "FAIL T3 smooth %g %g n=%d\n", n == 4 ? v[0] : -1, n == 4 ? v[3] : -1, n); ok = false;
        }
    }
    // ---- T4: archive + ump files written.
    {
        setup();
        archive_fn = "C3_arch.txt"; ump_fn = "C3_ump.txt";
        grid();
        std::ifstream af("C3_arch.txt");
        bool aok = af.good();
        std::string s;
        while (af.good() && std::getline(af, s)) {}
        af.close();
        std::ifstream uf("C3_ump.txt");
        int lines = 0;
        while (uf.good() && std::getline(uf, s)) ++lines;
        uf.close();
        if (!aok || lines != 4) { std::fprintf(stderr, "FAIL T4 files arch=%d lines=%d\n", aok, lines); ok = false; }
    }
    // ---- T5: tleft<0 mid-grid saves restart file and returns early.
    {
        setup();
        restart_fn = "C3_restart.bin";
        // Have the energy stub trip tleft on the first point.
        // (Simulate by pre-setting tleft small; grid sets tleft=-100 only when loop>=maxcyc.)
        tleft = -100.0;
        grid();
        std::ifstream rf("C3_restart.bin", std::ios::binary);
        bool rok = rf.good();
        if (rok) {
            rf.seekg(0, std::ios::end);
            long sz = (long)rf.tellg();
            rf.close();
            // 6 ints + 2 doubles + 1 point x (3*2 doubles + 3*2 ints + 3 doubles)
            long exp = 6 * 4 + 2 * 8 + (3 * 2 * 8 + 3 * 2 * 4 + 3 * 8);
            if (sz != exp) { std::fprintf(stderr, "FAIL T5 restart size=%ld exp=%ld\n", sz, exp); ok = false; }
        } else {
            std::fprintf(stderr, "FAIL T5 restart file missing\n"); ok = false;
        }
    }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
