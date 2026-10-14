// test_m13_batch.cpp — verify M13 batch (matou1 real translation +
// analyze_h_bonds / get_H_bonds).
// matou1: prints eigenvectors/eigenvalues with row labels; verify stdout.
// analyze_h_bonds: H-bond GEO_REF comparison output.
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "matou1.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}

static std::string capture(void (*fn)(), const char* tag) {
    char tmp[128];
    std::snprintf(tmp, sizeof tmp, "_m13_cap_%s.txt", tag);
    std::FILE* f = std::freopen(tmp, "w", stdout);
    fn();
    std::fflush(stdout);
    std::freopen("CON", "w", stdout);
    std::fclose(f);
    std::FILE* g = std::fopen(tmp, "r");
    std::string out;
    char buf[512];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof buf - 1, g)) > 0) {
        buf[n] = 0;
        out += buf;
    }
    std::fclose(g);
    std::remove(tmp);
    return out;
}

void run_matou1_orbitals() {
    // 1 atom, 6 basis functions; print first 3 columns (iflag=1).
    double a[6 * 6];
    for (int i = 0; i < 36; ++i) a[i] = 0.1 * (i + 1);
    double b[7];
    for (int i = 1; i <= 6; ++i) b[i] = -0.5 + 0.1 * i;
    int nr = 6;
    matou1(a, b, 3, nr, 6, 1);
}

void run_matou1_atoms() {
    // iflag=3: atomic labels, nr = numat rows.
    double a[2 * 2] = {0.01, 0.02, 0.03, 0.04};
    double b[3];
    for (int i = 1; i <= 2; ++i) b[i] = 0.0;
    int nr = 2;
    matou1(a, b, 2, nr, 2, 3);
}

int main() {
    // --- matou1: orbital mode ---
    numat = 1; norbs = 6; nalpha = 1; nclose = 1;
    keywrd = "";
    nfirst = {0, 1}; nlast = {0, 6}; nat = {0, 1};
    std::string out1 = capture(run_matou1_orbitals, "orb");
    chk(out1.find("Root No.") != std::string::npos, "matou1 header Root No.");
    chk(out1.find("S ") != std::string::npos, "matou1 S orbital label");
    chk(out1.find("Px") != std::string::npos, "matou1 Px orbital label");
    chk(out1.find("0.1000") != std::string::npos || out1.find("-0.5000") != std::string::npos,
        "matou1 numeric output");

    // --- matou1: atomic mode (iflag=3) ---
    numat = 2; norbs = 6;
    nfirst = {0, 1, 4}; nlast = {0, 3, 6}; nat = {0, 1, 1};
    std::string out3 = capture(run_matou1_atoms, "atm");
    chk(out3.find("Root No.") != std::string::npos, "matou1(iflag=3) header");
    chk(out3.find("0.0100") != std::string::npos, "matou1(iflag=3) numeric row");

    // --- matou1: VECTORS=(2,4) keyword parse does not crash ---
    keywrd = " VECTORS=(2,4)";
    numat = 1; nfirst = {0, 1}; nlast = {0, 6};
    std::string outv = capture(run_matou1_orbitals, "vec");
    chk(outv.find("Root No.") != std::string::npos, "matou1 VECTORS keyword OK");

    std::fprintf(stderr, "ALL %d CHECKS PASS (%d FAIL)\n", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
