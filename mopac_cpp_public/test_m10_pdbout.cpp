// test_m10_pdbout.cpp — batch tests for pdbout.F90 translation (M10 batch 10):
// PDB output driver (ATOM records, element symbols, TER, END, HEADER/REMARK).
#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "pdbout.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace molkst_C;

static int g_checks = 0;
static void chk(int cond, const char* msg) {
  ++g_checks;
  if (!cond) { std::fprintf(stderr, "[FAIL] %s\n", msg); std::exit(1); }
}

int main() {
  std::fprintf(stderr, "[t1] pdbout PDB records (H2O)\n");
  {
    numat = 2; natoms = 2;
    nat.assign(3, 0);
    nat[1] = 1; nat[2] = 8;
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 0.96; coord[2][2] = 0.0; coord[3][2] = 0.0;
    ncomments = 0;
    maxtxt = 0;
    input_fn = "water.pdb";
    std::string pad(30, ' ');
    txtatm.assign(3, pad);
    txtatm1.assign(3, pad);
    // PDB-style line: "ATOM  " (0-5), atom name at cols 12-26 (idx 11..), TER at 17..
    txtatm[1] = "ATOM      1  H1  MOL  " + std::string(14, ' ');
    txtatm[2] = "ATOM      2  O   MOL  " + std::string(14, ' ');
    breaks.assign(3, 0);
    nbreaks = 1;
    keywrd = "";
    iw = 6;
    std::fflush(stdout);
    std::string out;
    {
      FILE* cap = std::freopen("_pdbout_cap.txt", "w", stdout);
      if (cap) {
        pdbout(1);
        std::fflush(stdout);
        std::fclose(stdout);
      }
      FILE* f = std::fopen("_pdbout_cap.txt", "r");
      if (f) {
        char buf[2048];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
      }
    }
    chk(out.find("HEADER") != std::string::npos, "pdbout HEADER line");
    chk(out.find("REMARK") != std::string::npos, "pdbout REMARK line");
    chk(out.find("ATOM") != std::string::npos, "pdbout ATOM records");
    chk(out.find(" H") != std::string::npos || out.find(" H1") != std::string::npos,
        "pdbout H element");
    chk(out.find(" O") != std::string::npos || out.find(" O ") != std::string::npos,
        "pdbout O element");
    chk(out.find("END") != std::string::npos, "pdbout END line");
    std::remove("_pdbout_cap.txt");
  }

  std::fprintf(stderr, "[t2] pdbout mode1=-2 (no header, iprt=abs)\n");
  {
    // Reuse t1 state but force a different mode; PDB records must still be
    // emitted (iprt = abs(mode1)); TER appears when i == breaks[1].
    numat = 2; natoms = 2;
    nat[1] = 1; nat[2] = 8;
    breaks.assign(3, 0);
    breaks[1] = 2;   // TER after atom 2
    nbreaks = 1;
    std::fflush(stdout);
    std::string out;
    {
      FILE* cap = std::freopen("_pdbout_cap2.txt", "w", stdout);
      if (cap) {
        pdbout(-2);
        std::fflush(stdout);
        std::fclose(stdout);
      }
      FILE* f = std::fopen("_pdbout_cap2.txt", "r");
      if (f) {
        char buf[2048];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
      }
    }
    chk(out.find("TER") != std::string::npos, "pdbout TER after chain break");
    chk(out.find("END") != std::string::npos, "pdbout END (mode -2)");
    std::remove("_pdbout_cap2.txt");
  }

  std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
  return 0;
}
