// modchg.cpp — C++ translation of MOPAC 2016 "modchg.F90".
//
//   MODCHG prints the charge due to (a) backbone residue atoms, and
//   (b) all side-chain atoms in each residue.
//   build_res_start_etc works out which residue each atom belongs to, and the
//   location of the first atom in each residue.
#include "modchg.h"
#include "find_salt_bridges.h"
#include "set_up_dentate.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

namespace {

// Fortran fixed-length character semantics: txtatm entries are *26 characters.
// seg26 returns the 1-based inclusive slice [a..b] of a 26-char string.
std::string seg26(const std::string& s, int a, int b) {
  std::string t = s;
  if ((int)t.size() < 26) t.resize(26, ' ');
  if (a < 1) a = 1;
  if (b > 26) b = 26;
  if (b < a) return std::string();
  return t.substr(a - 1, b - a + 1);
}

}  // namespace

void build_res_start_etc() {
  using namespace molkst_C;
  using namespace MOZYME_C;
  using namespace common_arrays_C;

  int i, j;
  std::string residue;
  std::vector<std::string> ren_name(numat + 1, std::string(9, ' '));

  nres = 0;
  if ((int)at_res.size() < numat + id + 1) at_res.resize(numat + id + 1);
  if ((int)res_start.size() < numat + id + 2) res_start.resize(numat + id + 2, 0);
  for (i = 1; i <= numat; ++i) {
    residue = seg26(txtatm[i], 18, 26);
    for (j = 1; j <= nres; ++j)
      if (ren_name[j] == residue) break;
    if (j > nres) {
      nres = nres + 1;
      ren_name[j] = residue;
      res_start[j] = i;
    }
    at_res[i] = j;
  }
}

void modchg() {
  using namespace molkst_C;
  using namespace MOZYME_C;
  using namespace common_arrays_C;

  static bool first = true;
  char nam_het[600][4], het_num[600][5];
  double res_charge = 0.0;
  int i, j, n_cat, n_ani, n_het = 0;
  int cations[50], anions[50];

  std::vector<double> work(numat + id + 1, 0.0);
  build_res_start_etc();
  if (q.empty()) {
    std::fprintf(stdout, "%s\n",
                 "Density matrix is not available.  Run a 1SCF.");
    return;
  }
  std::vector<bool> l_used(numat + id + 1, false);

  // WORK collects charge contributions from atoms.  The negative addresses
  // refer to the residue backbone, the positive addresses refer to the
  // side-chain.  (Backbone atoms are -NH-CH-CO-)
  for (i = 1; i <= numat; ++i) {
    j = at_res[i];
    work[j] = work[j] + q[i];
    l_used[i] = (j > 0);
  }
  std::fprintf(stdout, "\n%9s\n\n", " NET CHARGE ON RESIDUES");
  std::fprintf(stdout, "%s\n", "      Residue         Charge  Anion or");
  std::fprintf(stdout, "%s\n", "                              Cation?");
  n_cat = 0;
  n_ani = 0;
  for (i = 1; i <= nres; ++i) {
    if (res_start[i] < 1) continue;
    if (work[i] > 0.5) {
      std::fprintf(stdout, "      %s%+13.3f %s\n",
                   seg26(txtatm[res_start[i]], 18, 26).c_str(), work[i],
                   " CATION");
      n_cat = (n_cat < 50) ? n_cat + 1 : 50;
      cations[n_cat - 1] = res_start[i];
    } else if (work[i] < -0.5) {
      std::fprintf(stdout, "      %s%13.3f %s\n",
                   seg26(txtatm[res_start[i]], 18, 26).c_str(), work[i],
                   " ANION");
      n_ani = (n_ani < 50) ? n_ani + 1 : 50;
      anions[n_ani - 1] = res_start[i];
    } else {
      if (work[i] < 0. && work[i] > -0.0005) work[i] = 0.;
      if (seg26(txtatm[res_start[i]], 18, 20) != "HOH" ||
          std::fabs(work[i]) > 0.1)
        std::fprintf(stdout, "      %s%+13.3f\n",
                     seg26(txtatm[res_start[i]], 18, 26).c_str(), work[i]);
    }
    res_charge = res_charge + work[i];
  }
  work.assign(numat + id + 1, 0.0);
  for (i = 1; i <= numat; ++i) {
    if (!l_used[i]) {
      for (j = 1; j <= n_het; ++j) {
        if (seg26(txtatm[i], 18, 20) == nam_het[j - 1] &&
            seg26(txtatm[i], 23, 26) == het_num[j - 1])
          break;
      }
      if (j > n_het) {
        n_het = n_het + 1;
        if (n_het > 600) break;  // Fortran: exit the i-loop
        std::string s1 = seg26(txtatm[i], 18, 20);
        std::string s2 = seg26(txtatm[i], 23, 26);
        std::memcpy(nam_het[n_het - 1], s1.c_str(), 3);
        nam_het[n_het - 1][3] = '\0';
        std::memcpy(het_num[n_het - 1], s2.c_str(), 4);
        het_num[n_het - 1][4] = '\0';
      }
      work[j] = work[j] + q[i];
    }
  }
  if (n_het > 0) {
    for (i = 1; i <= n_het; ++i) {
      if (std::fabs(work[i]) < 0.1) continue;
      if (first) {
        first = false;
        std::fprintf(stdout, "\n%6s\n%2s\n\n",
                     "NET CHARGES ON HETERO-GROUPS  (|Charges| less than 0.1 "
                     "not printed)",
                     "");
        std::fprintf(stdout, "%s\n", "    HETERO-GROUP      Charge ");
        std::fprintf(stdout, "%s\n", " ");
      }
      std::fprintf(stdout, "      %s%14.3f\n",
                   seg26(txtatm[res_start[i]], 18, 26).c_str(), work[i]);
    }
  }
  std::fprintf(stdout, "\n");

  //  Print out salt bridges
  if (n_cat > 0 && n_ani > 0) {
    set_up_dentate();
    find_salt_bridges(cations, anions, n_cat, n_ani);
  }
  std::fprintf(stdout, "\n\n");
}
