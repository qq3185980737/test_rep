// test_mpcsyb.cpp — verify SYBYL output driver + Mulliken populations.
// H2 scenario: numat=2, norbs=2, single s orbital per atom.
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include "mpcsyb.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "parameters_C.h"
using namespace molkst_C;
using namespace common_arrays_C;
using namespace chanel_C;
using namespace parameters_C;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}

int main() {
    numat = 2; norbs = 2; nclose = 1; nalpha = 1; nbeta = 1;
    escf = -0.5234;
    nfirst = {0, 1, 2}; nlast = {0, 1, 2};
    nat = {0, 1, 1};
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[0][2] = 1.4;
    eigs.assign(6, 0.0); eigs[1] = -0.9; eigs[2] = -0.7; eigs[3] = 0.2; eigs[4] = 0.8;
    pb.assign(10, 0.0); pb[1] = 0.6; pb[3] = 0.4;
    keywrd = " SYBYL MULL";
    syb_fn = "_m11_syb_test.syb";
    iw = 1;  // stdout

    double chr[3] = {0.0, 0.2, -0.2};
    double dip = 1.35;
    mpcsyb(chr, 1, 0.4123, dip);
    chk(dip == 0.0, "kchrge!=0 zeroes dipole");

    FILE* f = std::fopen(syb_fn.c_str(), "r");
    chk(f != nullptr, "SYBYL file written");
    if (!f) return 1;
    char line[256];
    chk(std::fgets(line, sizeof(line), f) != nullptr, "header line");
    chk(std::strncmp(line, "   1   2", 8) == 0, "header: flag=1 natoms=2");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "atom1 line");
    char a1[24], a2[24], a3[24], a4[24];
    std::sscanf(line, "%s %s %s %s", a1, a2, a3, a4);
    chk(std::atof(a1) == 0.0 && std::atof(a4) == 0.2, "atom1 coords+charge");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "atom2 line");
    std::sscanf(line, "%s %s %s %s", a1, a2, a3, a4);
    chk(std::atof(a1) == 1.4 && std::atof(a4) == -0.2, "atom2 coords+charge");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "HOMO line");
    chk(std::strstr(line, "HOMOs,LUMOs,# of occupied MOs") != nullptr, "HOMO label");
    double e1, e2, e3, e4; int nf;
    std::sscanf(line, "%lf %lf %lf %lf %d", &e1, &e2, &e3, &e4, &nf);
    chk(e1 == -0.9 && e4 == 0.8 && nf == 1, "HOMO energies + nfilled");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "HF line");
    double hf, ip;
    std::sscanf(line, "%lf %lf", &hf, &ip);
    chk(hf == -0.5234 && ip == 0.4123, "HF and IP values");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "dipole line");
    int kc; double dp;
    std::sscanf(line, "%d %lf", &kc, &dp);
    chk(kc == 1 && dp == 0.0, "charge+dipole (zeroed)");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "mulliken header");
    chk(std::strstr(line, "MULLIKEN POPULATION AND CHARGE") != nullptr, "mulliken label");
    chk(std::fgets(line, sizeof(line), f) != nullptr, "pop line");
    double p1, p2, c1, c2;
    std::sscanf(line, "%lf %lf %lf %lf", &p1, &p2, &c1, &c2);
    chk(p1 == 0.6 && p2 == 0.4, "populations");
    chk(c1 == 0.4 && c2 == 0.6, "charges (tore - pop)");
    std::fclose(f);
    std::remove(syb_fn.c_str());
    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
