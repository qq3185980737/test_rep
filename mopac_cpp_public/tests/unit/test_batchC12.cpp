// test_batchC12.cpp — batch tests for 批C12 (ligand + distance).
// T1 water HOH; T2 phosphate PO4; T3 nheavy; T4 moiety; T5 distance;
// T6 methanol MOH.
#include "ligand.h"
#include "distance.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace MOZYME_C;

static int cal_ctr = 0;
static int failures = 0;

static void check(const char* name, bool ok) {
  if (!ok) {
    std::cerr << "FAIL: " << name << "\n";
    ++failures;
  } else {
    std::cout << "ok: " << name << "\n";
  }
}

static void base_setup(int n, const int* z, const int* lab) {
  numcal = ++cal_ctr;
  moperr = false;
  id = 0;
  keywrd = " ";
  ncomments = 0;
  natoms = n;
  numat = n;
  nat.assign(n + 1, 0);
  labels.assign(n + 1, 0);
  nbonds.assign(n + 1, 0);
  ibonds.assign(16, std::vector<int>(n + 1, 0));
  coord.assign(4, std::vector<double>(n + 1, 0.0));
  txtatm.assign(n + 1, std::string(26, ' '));
  txtatm1.assign(n + 1, std::string(26, ' '));
  all_comments.assign(40, std::string(81, ' '));
  allres.assign(40, std::string(3, ' '));
  for (int i = 1; i <= n; ++i) {
    nat[i] = z[i - 1];
    labels[i] = lab ? lab[i - 1] : z[i - 1];
  }
}

static void bond(int a, int b) {
  nbonds[a] = nbonds[a] + 1;
  nbonds[b] = nbonds[b] + 1;
  ibonds[nbonds[a]][a] = b;
  ibonds[nbonds[b]][b] = a;
}

static void t_water() {
  int z[3] = {8, 1, 1};
  base_setup(3, z, nullptr);
  bond(1, 2);
  bond(1, 3);
  std::vector<int> start_res(4, -200);
  int ires = 0, nfrag = 1;
  ligand(ires, start_res, nfrag);
  check("T1 ires=1", ires == 1);
  check("T1 HOH on O", txtatm[1].substr(15, 3) == "HOH");
  check("T1 HOH on H", txtatm[2].substr(15, 3) == "HOH");
  check("T1 tag HETATM", txtatm[1].substr(0, 6) == "HETATM");
}

static void t_phosphate() {
  int z[5] = {15, 8, 8, 8, 8};
  base_setup(5, z, nullptr);
  bond(1, 2);
  bond(1, 3);
  bond(1, 4);
  bond(1, 5);
  std::vector<int> start_res(4, -200);
  int ires = 0, nfrag = 1;
  ligand(ires, start_res, nfrag);
  check("T2 PO4 element", txtatm[1].substr(12, 3) == " P ");
  check("T2 PO4 residue", txtatm[1].substr(15, 3) == "PO4");
  check("T2 PO4 on O", txtatm[2].substr(12, 3) == " O ");
}

static void t_nheavy() {
  int z[3] = {6, 6, 1};
  base_setup(3, z, nullptr);
  bond(1, 2);
  bond(1, 3);
  check("T3 nheavy(C)=1", nheavy(1) == 1);
  check("T3 nheavy(H)=1 (C neighbor)", nheavy(3) == 1);
}

static void t_moiety() {
  // C1-C2-N3-O4 linear, C1 also attached to H5 (all bonds given manually).
  int z[5] = {6, 6, 7, 8, 1};
  base_setup(5, z, nullptr);
  bond(1, 2);
  bond(2, 3);
  bond(3, 4);
  bond(1, 5);
  std::vector<bool> iopt(6, false);
  std::vector<int> lused(6, 0);
  int n_new = 0;
  moiety(iopt, lused, 1, n_new);
  check("T4 moiety count=5", n_new == 5);
  check("T4 all marked", iopt[1] && iopt[2] && iopt[3] && iopt[4] && iopt[5]);
  // H must be last in the lused list
  bool h_last = true;
  for (int i = 1; i <= n_new; ++i)
    if (nat[lused[i]] == 1 && i != n_new) h_last = false;
  check("T4 H last", h_last);
}

static void t_distance() {
  int z[2] = {6, 6};
  base_setup(2, z, nullptr);
  coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
  coord[1][2] = 3.0; coord[2][2] = 4.0; coord[3][2] = 0.0;
  check("T5 distance 3-4-5", std::fabs(distance(1, 2) - 5.0) < 1e-12);
  // torsion sanity via dihed path: all-in-plane -> 0 or pi
  double tor = torsion(1, 1, 1, 2);  // degenerate; just must not crash
  (void)tor;
}

static void t_methanol() {
  // C(1) bonded to O(2), H(3), H(4), H(5); O(2) bonded to H(6).
  int z[6] = {6, 8, 1, 1, 1, 1};
  base_setup(6, z, nullptr);
  bond(1, 2);
  bond(1, 3);
  bond(1, 4);
  bond(1, 5);
  bond(2, 6);
  std::vector<int> start_res(4, -200);
  int ires = 0, nfrag = 1;
  ligand(ires, start_res, nfrag);
  // C has 1 C? no: n_C=1 (only C itself), n_O=1 -> "MOH" (methanol)
  check("T6 MOH residue", txtatm[1].substr(14, 3) == "MOH");
  check("T6 C element", txtatm[1].substr(11, 3) == " C ");
  check("T6 O element", txtatm[2].substr(11, 3) == " O ");
  check("T6 ires advanced", ires >= 1);
}

int main() {
  t_water();
  t_phosphate();
  t_nheavy();
  t_moiety();
  t_distance();
  t_methanol();
  if (failures == 0) {
    std::cout << "ALL PASS\n";
    return 0;
  }
  std::cerr << failures << " FAILURES\n";
  return 1;
}
