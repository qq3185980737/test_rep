// find_salt_bridges.cpp — C++ translation of geochk.F90
// "subroutine find_salt_bridges" (lines 3106-3490).  Identifies all
// ionizable sites that can be used for making salt bridges and rewrites the
// SITE keyword accordingly.  1-based Fortran indexing is kept.
#include "find_salt_bridges.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "geochk.h"
#include "common_arrays_C.h"
#include "distance.h"
#include "molkst_C.h"
#include "mopend.h"
#include "reada.h"

namespace mc = common_arrays_C;
namespace mk = molkst_C;
namespace ch = chanel_C;

static std::size_t len_trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? 0 : e + 1;
}

static std::string fsub(const std::string& s, int a, int b) {
  if (a < 1) a = 1;
  if (b > (int)s.size()) b = (int)s.size();
  if (a > b) return "";
  return s.substr(a - 1, b - a + 1);
}

static void setsub(std::string& s, int a, int b, const std::string& v) {
  if (a < 1) a = 1;
  if (b > (int)s.size()) s.resize(b, ' ');
  int n = b - a + 1;
  std::string t = v;
  if ((int)t.size() > n) t.resize(n);
  while ((int)t.size() < n) t += ' ';
  s.replace(a - 1, n, t);
}

static int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// txtatm(i)(22:22) in 'A'..'Z'?
static bool chain_letter_ok(int i) {
  char c = fsub(mc::txtatm[i], 22, 22)[0];
  return c >= 'A' && c <= 'Z';
}

void find_salt_bridges(int in_cat[], int in_ani[], int n_cat, int n_ani) {
  int i, j, k, l, m, n, n_cations, n_anions, n_C, n_H, n_O, C;
  int pairs[3][201] = {{0}};          // pairs(2,200), 1-based
  int n_pairs, n_salt;
  int salt_bridges[3][201] = {{0}};   // salt_bridges(2,200), 1-based
  int i_length, ii, jj, nh_cat, nh_ani;
  double r, cutoff;
  std::vector<double> Rab(201, 0.0);
  double R_min;
  std::vector<double> R_sorted(201, 0.0);
  std::string bits(10, ' ');

  std::vector<int> cations(mk::numat + 1, 0), anions(mk::numat + 1, 0);
  if (n_cat == 0) {
    i = index1(mk::keywrd, "SALT=");
    if (i != 0) {
      cutoff = reada(mk::keywrd, i + 4);
    } else {
      cutoff = 4.0;
    }
    nh_cat = 0;
    nh_ani = 0;
  } else {
    cutoff = 8.0;
    nh_cat = 1;
    nh_ani = -1;
  }

  // Locate all Arg guanidine "C" atoms and all -C(R)H-NH2 "N" atoms.
  n_cations = 0;
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] == 6) {
      if (mc::nbonds[i] == 3) {
        if (mc::nat[mc::ibonds[1][i]] == 7 && mc::nat[mc::ibonds[2][i]] == 7 &&
            mc::nat[mc::ibonds[3][i]] == 7) {
          // Make sure it has the correct charge.
          m = 0;
          for (j = 1; j <= 3; ++j) {
            n = mc::ibonds[j][i];
            for (k = 1; k <= mc::nbonds[n]; ++k) {
              if (mc::nat[mc::ibonds[k][n]] == 1) m = m + 1;
            }
          }
          if (m == 4 + nh_cat) {
            if (n_cat > 0) {
              for (n = 1; n <= n_cat; ++n) {
                if (fsub(mc::txtatm[in_cat[n]], 18, 26) == fsub(mc::txtatm[i], 18, 26)) break;
              }
            } else {
              n = 0;
            }
            if (n <= n_cat && fsub(mc::txtatm[i], 18, 20) != "UNK" && chain_letter_ok(i)) {
              n_cations = n_cations + 1;
              cations[n_cations] = i;
            }
          }
        }
      }
    }
  }
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] == 7) {
      if (mc::nbonds[i] == 3 + nh_cat) {
        // Test for -NH2, e.g. N-terminus and lysine.
        n_C = 0;
        n_H = 0;
        for (j = 1; j <= 3 + nh_cat; ++j) {
          if (mc::nat[mc::ibonds[j][i]] == 6 &&
              fsub(mc::txtatm[mc::ibonds[j][i]], 18, 20) != "UNK" &&
              chain_letter_ok(mc::ibonds[j][i])) {
            n_C = n_C + 1;
            C = mc::ibonds[j][i];
          }
          if (mc::nat[mc::ibonds[j][i]] == 1) n_H = n_H + 1;
        }
        // Make sure it has the correct charge.
        if (n_C == 1 && n_H == 2 + nh_cat) {
          for (j = 1; j <= mc::nbonds[C]; ++j) {
            if (mc::nat[mc::ibonds[j][C]] == 8) break;
            if (mc::nat[mc::ibonds[j][C]] == 7 && mc::ibonds[j][C] != i) break;
          }
          if (j > mc::nbonds[C]) {
            if (n_cat > 0) {
              for (n = 1; n <= n_cat; ++n) {
                if (fsub(mc::txtatm[in_cat[n]], 18, 26) == fsub(mc::txtatm[i], 18, 26)) break;
              }
            } else {
              n = 0;
            }
            if (n <= n_cat && fsub(mc::txtatm[C], 18, 20) != "UNK" && chain_letter_ok(C)) {
              n_cations = n_cations + 1;
              cations[n_cations] = C;
            }
          }
        }
      } else if (mc::nbonds[i] == 2 + nh_cat && fsub(mc::txtatm[i], 18, 20) == "HIS") {
        // Check for Histidine.
        for (j = 1; j <= mc::nbonds[i]; ++j) {
          k = mc::ibonds[j][i];
          if (mc::nbonds[k] < 3) continue;
          for (l = 1; l <= mc::nbonds[k]; ++l) {
            if (mc::nat[mc::ibonds[l][k]] == 7 && mc::ibonds[l][k] != i) break;
          }
          if (l <= mc::nbonds[k]) break;
        }
        if (j > mc::nbonds[i]) continue;
        if (n_cat > 0) {
          for (n = 1; n <= n_cat; ++n) {
            if (fsub(mc::txtatm[in_cat[n]], 18, 26) == fsub(mc::txtatm[i], 18, 26)) break;
          }
        } else {
          n = 0;
        }
        // Avoid a "double count".
        for (k = 1; k <= n_cations; ++k) {
          if (fsub(mc::txtatm[cations[k]], 18, 26) == fsub(mc::txtatm[mc::ibonds[j][i]], 18, 26)) break;
        }
        if (n <= n_cat && k > n_cations) {
          n_cations = n_cations + 1;
          cations[n_cations] = mc::ibonds[j][i];
        }
      }
    }
  }

  // Locate all -COOH "C" atoms.
  n_anions = 0;
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] == 6) {
      if (mc::nbonds[i] == 3) {
        n_C = 0;
        n_O = 0;
        for (j = 1; j <= 3; ++j) {
          if (mc::nat[mc::ibonds[j][i]] == 6) n_C = n_C + 1;
          if (mc::nat[mc::ibonds[j][i]] == 8) n_O = n_O + 1;
        }
        if (n_C == 1 && n_O == 2) {
          m = 0;
          ii = 0;
          for (j = 1; j <= 3; ++j) {
            if (mc::nat[mc::ibonds[j][i]] == 8) {
              l = mc::ibonds[j][i];
              ii = ii + mc::nbonds[l];
              for (k = 1; k <= mc::nbonds[l]; ++k) {
                if (mc::nat[mc::ibonds[k][l]] == 1) m = m + 1;
              }
            }
          }
          if (m == 1 + nh_ani && ii < 4) {
            if (n_ani > 0) {
              for (n = 1; n <= n_ani; ++n) {
                if (fsub(mc::txtatm[in_ani[n]], 18, 26) == fsub(mc::txtatm[i], 18, 26)) break;
              }
            } else {
              n = 0;
            }
            if (n <= n_ani && fsub(mc::txtatm[i], 18, 20) != "UNK" && chain_letter_ok(i)) {
              n_anions = n_anions + 1;
              anions[n_anions] = i;
            }
          }
        }
      }
    }
  }

  // Find interatomic distances and weight them by type.
  n_pairs = 0;
  for (i = 1; i <= n_cations; ++i) {
    for (j = 1; j <= n_anions; ++j) {
      r = distance(cations[i], anions[j]);
      if (r < cutoff + 2.0) {
        m = cations[i];
        n = anions[j];
        R_min = 100.0;
        for (k = 1; k <= mc::nbonds[m]; ++k) {
          if (mc::nat[mc::ibonds[k][m]] != 7) continue;
          for (l = 1; l <= mc::nbonds[n]; ++l) {
            if (mc::nat[mc::ibonds[l][n]] != 8) continue;
            r = distance(mc::ibonds[k][m], mc::ibonds[l][n]);
            if (r < R_min) {
              R_min = r;
              ii = mc::ibonds[k][m];
              jj = mc::ibonds[l][n];
            }
          }
        }
        if (R_min < cutoff) {
          n_pairs = n_pairs + 1;
          pairs[1][n_pairs] = ii;
          pairs[2][n_pairs] = jj;
          Rab[n_pairs] = R_min;
        }
      }
    }
  }

  // Sort interatomic distances.
  n_salt = 0;
  while (true) {
    if (n_pairs == 0) break;
    r = 1.e10;
    for (i = 1; i <= n_pairs; ++i) {
      if (r > Rab[i]) {
        r = Rab[i];
        j = i;
      }
    }
    n_salt = n_salt + 1;
    salt_bridges[1][n_salt] = pairs[1][j];
    salt_bridges[2][n_salt] = pairs[2][j];
    R_sorted[n_salt] = r;
    n_pairs = n_pairs - 1;
    for (k = j; k <= n_pairs; ++k) {
      Rab[k] = Rab[k + 1];
      pairs[1][k] = pairs[1][k + 1];
      pairs[2][k] = pairs[2][k + 1];
    }
  }
  if (n_cat == 0) {
    // Eliminate duplicates.
    for (i = 1; i <= n_salt; ++i) {
      if (salt_bridges[1][i] == 0) continue;
      for (j = n_salt; j >= i + 1; --j) {
        if (salt_bridges[1][j] == 0) continue;
        m = salt_bridges[1][i];
        n = salt_bridges[2][i];
        ii = salt_bridges[1][j];
        jj = salt_bridges[2][j];
        if (salt_bridges[1][j] == salt_bridges[1][i] ||
            salt_bridges[2][j] == salt_bridges[2][i]) {
          salt_bridges[1][j] = 0;
          salt_bridges[2][j] = 0;
        } else if (fsub(mc::txtatm[m], 18, 26) == fsub(mc::txtatm[jj], 18, 26)) {
          salt_bridges[1][j] = 0;
          salt_bridges[2][j] = 0;
        } else {
          for (ii = 1; ii <= 2; ++ii) {
            m = salt_bridges[ii][i];
            n = salt_bridges[ii][j];
            for (k = 1; k <= mc::nbonds[m]; ++k) {
              if (mc::nat[mc::ibonds[k][m]] != 6) continue;
              for (l = 1; l <= mc::nbonds[n]; ++l) {
                if (mc::nat[mc::ibonds[l][n]] != 6) continue;
                if (mc::ibonds[k][m] == mc::ibonds[l][n]) {
                  salt_bridges[1][j] = 0;
                  salt_bridges[2][j] = 0;
                }
              }
            }
            if (salt_bridges[1][j] == 0) break;
          }
        }
      }
    }
  }

  // Build new SITE keyword.
  mk::line = " ";
  if (n_salt > 0) {
    std::fprintf(stdout, "\n\n%19s\n", "Salt Bridges Found");
    std::fprintf(stdout, "\n%4s%10s%28s%16s\n", " No.", "Cationic site", "Anionic site", "Dist. (Angstroms)");
    if (n_cat == 0) update_txtatm(true, true);
  } else {
    if (n_cat == 0) {
      mopend("SALT option used in SITE command, but no Salt Bridges found");
      mk::moperr = false;
      return;
    } else {
      std::fprintf(stdout, "\n\n%19s\n", "No Salt Bridges Found");
    }
  }
  j = 0;
  for (i = 1; i <= n_salt; ++i) {
    if (salt_bridges[1][i] == 0) continue;
    k = salt_bridges[1][i];
    for (m = 26; m >= 23; --m) {
      if (fsub(mc::txtatm[k], m, m) == " ") break;
    }
    m = (23 > m + 1) ? 23 : (m + 1);
    l = salt_bridges[2][i];
    for (n = 26; n >= 23; --n) {
      if (fsub(mc::txtatm[l], n, n) == " ") break;
    }
    n = (23 > n + 1) ? 23 : (n + 1);
    i_length = (int)len_trim(mk::line);
    if (i_length > 900) break;
    std::string cat_txt = std::string(1, fsub(mc::txtatm[k], 22, 22)[0]) +
                          fsub(mc::txtatm[k], m, 26) + "(+)";
    std::string ani_txt = std::string(1, fsub(mc::txtatm[l], 22, 22)[0]) +
                          fsub(mc::txtatm[l], n, 26) + "(-)";
    setsub(mk::line, i_length + 1, i_length + 1 + (int)cat_txt.size() + 1 + (int)ani_txt.size() - 1,
           "," + cat_txt + "," + ani_txt);
    j = j + 1;
    std::fprintf(stdout, "%7d%2s%-9s%4s%-9s%10.2f\n", j,
                 ("(" + mc::txtatm[k] + ")").c_str(), cat_txt.c_str(),
                 ("(" + mc::txtatm[l] + ")").c_str(), ani_txt.c_str(), R_sorted[i]);
    if (fsub(mc::txtatm[k], 22, 22)[0] < 'A' || fsub(mc::txtatm[k], 22, 22)[0] > 'Z') {
      mk::line = "Chain letter for the cation is not in the range 'A' to 'Z'";
      mopend(mk::line.substr(0, len_trim(mk::line)));
      std::fprintf(stdout, "%10s\n", ("Chain letter = \"" + std::string(1, fsub(mc::txtatm[k], 22, 22)[0]) + "\"").c_str());
      return;
    }
    if (fsub(mc::txtatm[l], 22, 22)[0] < 'A' || fsub(mc::txtatm[l], 22, 22)[0] > 'Z') {
      mk::line = "Chain letter for the anion is not in the range 'A' to 'Z'";
      mopend(mk::line.substr(0, len_trim(mk::line)));
      std::fprintf(stdout, "%10s\n", ("Chain letter = \"" + std::string(1, fsub(mc::txtatm[l], 22, 22)[0]) + "\"").c_str());
      return;
    }
  }
  if (n_cat > 0) return;
  mk::line = mk::line.substr(1);

  // Delete the word "SALT".
  k = index1(mk::keywrd, "(SALT") + index1(mk::keywrd, ",SALT") + 1;
  i = index1(mk::keywrd, "(SALT=") + index1(mk::keywrd, ",SALT=");
  if (i > 0) {
    for (j = i; j <= (int)len_trim(mk::keywrd); ++j) {
      if (fsub(mk::keywrd, j, j) == "," || fsub(mk::keywrd, j, j) == ")") break;
    }
    mk::keywrd.replace(i, mk::keywrd.size() - i, mk::keywrd.substr(j - 1));
    if (fsub(mk::keywrd, i, i) == ",") {
      std::string t = mk::keywrd.substr(i);
      while (!t.empty() && t[0] == ' ') t = t.substr(1);
      mk::keywrd.replace(i - 1, mk::keywrd.size() - (i - 1), t);
    }
  } else {
    i = index1(mk::keywrd, "(SALT") + index1(mk::keywrd, ",SALT") + 1;
    mk::keywrd.replace(i - 1, mk::keywrd.size() - (i - 1), mk::keywrd.substr(i + 3));
    if (fsub(mk::keywrd, i, i) == ",") {
      std::string t = mk::keywrd.substr(i);
      while (!t.empty() && t[0] == ' ') t = t.substr(1);
      mk::keywrd.replace(i - 1, mk::keywrd.size() - (i - 1), t);
    }
  }

  // Check for duplicate sites in keywrd.
  m = index1(mk::keywrd, " SITE");
  m = index1(mk::keywrd.substr(m - 1), "(") + m;
  j = index1(mk::keywrd.substr(m - 1), ") ") + m;
  if (fsub(mk::keywrd, m, m) == ",") m = m + 1;
  while (true) {
    k = index1(mk::keywrd.substr(m - 1, j - m), ",");
    if (k == 0) k = index1(mk::keywrd.substr(m - 1, j - m + 1), ") ");
    if (k == 0) break;
    k = k + m - 2;
    bits = mk::keywrd.substr(m - 1, k - m + 1);
    if (k - m < 3) break;
    l = index1(mk::line, bits.substr(0, len_trim(bits)));
    if (l > 0) {
      // Remove the comma.
      if (fsub(mk::line, l - 1, l - 1) == ",") {
        mk::line = mk::line.substr(0, l - 2) + mk::line.substr(l + k - m);
      } else {
        mk::line = mk::line.substr(0, l - 1) + mk::line.substr(l + k - m - 1);
      }
    }
    m = k + 2;
    if (fsub(mk::keywrd, m, m) == ",") m = m + 1;
  }
  // Insert the new keyword.
  i = index1(mk::keywrd, " SITE");
  if (fsub(mk::keywrd, i + 7, i + 7) != ")") {
    mk::keywrd = mk::keywrd.substr(0, i + 7 - 1) +
                 mk::line.substr(0, len_trim(mk::line)) + "," +
                 mk::keywrd.substr(i + 7);
  } else {
    mk::keywrd = mk::keywrd.substr(0, i + 7 - 1) +
                 mk::line.substr(0, len_trim(mk::line)) +
                 mk::keywrd.substr(i + 6);
  }
}
