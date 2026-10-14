// test_batchC7.cpp — regression tests for modchg (batch C7).
// Each case redirects stdout to a capture file; assertions read it back.
#include "modchg.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace MOZYME_C;
using namespace common_arrays_C;

static int cal_ctr = 0;

// Build a 26-char txtatm entry (PDB-like): atom name at 13-16, residue at
// 18-20, residue number at 23-26 (all 1-based Fortran positions).
static std::string mk(const std::string& atom, const std::string& res,
                      const std::string& num) {
  std::string s = "ATOM";
  s.resize(12, ' ');
  s += atom;
  s.resize(17, ' ');
  s += res;
  s.resize(22, ' ');
  s += num;
  s.resize(26, ' ');
  return s;
}

static void setup(int numat_n, const std::vector<std::string>& labels) {
  numcal = ++cal_ctr;
  keywrd = " ";
  moperr = false;
  mozyme = false;
  numat = numat_n;
  id = 0;
  q.assign(numat_n + 1, 0.0);
  txtatm.assign(numat_n + 1, std::string(26, ' '));
  for (int i = 1; i <= numat_n; ++i) txtatm[i] = labels[i - 1];
  nres = 0;
  at_res.clear();
  res_start.clear();
}

static std::string read_file(const char* path) {
  FILE* f = std::fopen(path, "rb");
  if (!f) return "";
  std::string s;
  char buf[4096];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) s.append(buf, n);
  std::fclose(f);
  return s;
}

static void T1_density_matrix_unavailable() {
  setup(2, {mk(" N ", "ASP", "  1"), mk(" C ", "ASP", "  1")});
  q.clear();  // emulate .not. allocated(q)
  modchg();
}

static void T2_residue_charges() {
  setup(4, {mk(" N ", "ASP", "  1"), mk(" C ", "ASP", "  1"),
            mk(" N ", "GLU", "  2"), mk(" O ", "HOH", "  3")});
  q[1] = 0.7;
  q[2] = 0.2;   // residue ASP total = +0.9  -> CATION
  q[3] = -0.8;  // residue GLU total = -0.8  -> ANION
  q[4] = 0.05;  // residue HOH  = +0.05 -> below threshold, not printed
  modchg();
}

static void T3_build_res_start_etc() {
  setup(4, {mk(" N ", "ASP", "  1"), mk(" C ", "ASP", "  1"),
            mk(" N ", "GLU", "  2"), mk(" O ", "HOH", "  3")});
  build_res_start_etc();
  if (nres != 3) {
    std::fprintf(stderr, "FAIL T3 nres=%d expect 3\n", nres); std::exit(1);
  }
  if (at_res[1] != 1 || at_res[2] != 1 || at_res[3] != 2 || at_res[4] != 3) {
    std::fprintf(stderr, "FAIL T3 at_res = %d %d %d %d expect 1 1 2 3\n",
                 at_res[1], at_res[2], at_res[3], at_res[4]);
    std::exit(1);
  }
  if (res_start[1] != 1 || res_start[2] != 3 || res_start[3] != 4) {
    std::fprintf(stderr, "FAIL T3 res_start = %d %d %d expect 1 3 4\n",
                 res_start[1], res_start[2], res_start[3]);
    std::exit(1);
  }
}

int main() {
  // T1: no density matrix -> early message.
  if (!std::freopen("test_c7_t1.txt", "w", stdout)) {
    std::fprintf(stderr, "FAIL cannot redirect stdout\n"); return 1;
  }
  T1_density_matrix_unavailable();
  std::fflush(stdout);
  std::string o1 = read_file("test_c7_t1.txt");
  if (o1.find("Density matrix is not available") == std::string::npos) {
    std::fprintf(stderr, "FAIL T1 missing density-matrix message\n"); return 1;
  }

  // T2: residue charge report.
  if (!std::freopen("test_c7_out.txt", "w", stdout)) {
    std::fprintf(stderr, "FAIL cannot redirect stdout\n"); return 1;
  }
  T2_residue_charges();
  std::fflush(stdout);
  std::string out = read_file("test_c7_out.txt");
  if (out.find(" NET CHARGE ON RESIDUES") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing residue header\n"); return 1;
  }
  if (out.find(" CATION") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing CATION\n"); return 1;
  }
  if (out.find(" ANION") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing ANION\n"); return 1;
  }
  if (out.find("+0.900") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing ASP charge +0.900\n"); return 1;
  }
  if (out.find("-0.800") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing GLU charge -0.800\n"); return 1;
  }
  if (out.find("HOH") != std::string::npos) {
    std::fprintf(stderr, "FAIL T2 HOH residue should not be printed\n"); return 1;
  }

  // T3: residue bookkeeping.
  T3_build_res_start_etc();

  std::fprintf(stderr, "ALL PASS\n");
  return 0;
}
