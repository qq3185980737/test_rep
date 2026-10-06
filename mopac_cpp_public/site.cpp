// site.cpp — C++ translation of geochk.F90 "subroutine site" (lines
// 2274-3105).  Ionizes / de-ionizes residues and adds or removes hydrogen
// atoms as requested by keyword SITE.  1-based Fortran indexing is kept.
#include "geochk.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "add_hydrogen_atoms.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "distance.h"
#include "elemts_C.h"
#include "mopend.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "reada.h"
#include "set_up_dentate.h"

namespace mc = common_arrays_C;
namespace mk = molkst_C;
namespace pc = parameters_C;
namespace ch = chanel_C;
namespace ec = elemts_C;

// Fortran len_trim.
static std::size_t len_trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? 0 : e + 1;
}

// Fortran substring s(a:b), 1-based inclusive.
static std::string fsub(const std::string& s, int a, int b) {
  if (a < 1) a = 1;
  if (b > (int)s.size()) b = (int)s.size();
  if (a > b) return "";
  return s.substr(a - 1, b - a + 1);
}

// Fortran write s(a:b) = v (1-based inclusive), padding short v with blanks.
static void setsub(std::string& s, int a, int b, const std::string& v) {
  if (a < 1) a = 1;
  if (b > (int)s.size()) s.resize(b, ' ');
  int n = b - a + 1;
  std::string t = v;
  if ((int)t.size() > n) t.resize(n);
  while ((int)t.size() < n) t += ' ';
  s.replace(a - 1, n, t);
}

static int nint(double x) { return static_cast<int>(std::floor(x + 0.5)); }

// Fortran "c == ' '" for a fixed-length character variable: all blanks.
static bool is_blank(const std::string& s) {
  for (char c : s) if (c != ' ') return false;
  return true;
}

// 1-based Fortran index().
static int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(i)(:6), natoms, txtatm(i)(15:)
// prefix6: chars 1..6 of the source label; suffix12: chars 15..26.
static std::string fmt_txtatm_H(const std::string& prefix6, int natoms,
                                const std::string& suffix12) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%6s%5d 3H%-12s", prefix6.c_str(), natoms,
                suffix12.c_str());
  return std::string(buf);
}

void site(const std::vector<bool>& neutral, const std::vector<char>& chain,
          const std::vector<int>& res, std::vector<std::vector<char>>& charge,
          int nres, int max_sites, std::string& allkey) {
  int i_atom, i_charge;
  int i, j, k, l, m, n, ii, o[3], jres;
  std::vector<int> metals(2, 0);
  int nmetals = 0, nadd, ndel;
  bool l_res, bug = false;
  std::string res_txt(14, ' '), txt(26, ' '), num(1, ' ');
  std::vector<std::string> changes(mk::numat + 1, std::string(30, ' '));

  i = index1(mk::keywrd, " SITE=");
  j = index1(mk::keywrd.substr(i - 1), ") ") + i;
  k = index1(allkey, "SALT");
  if (k != 0) {
    k = index1(allkey, " SITE");
    l = index1(allkey.substr(k - 1), ") ") + k - 1;
    allkey = " \"" + allkey.substr(k, l - k) + "\" => " +
             mk::keywrd.substr(i - 1, j - i + 1);
  }

  // Check for duplicated residues.
  for (i = 1; i <= nres; ++i) {
    for (j = i + 1; j <= nres; ++j) {
      if (chain[i] == chain[j] && res[i] == res[j]) {
        num = std::string(1, (char)('1' + (int)std::floor(std::log10(res[j] * 1.01))));
        char b[32];
        std::snprintf(b, sizeof b, "%c%d", chain[i], res[i]);
        setsub(txt, 1, 26, b);
        mk::line = "Residue " + txt.substr(1) + " in chain " +
                   std::string(1, chain[i]) + " occurs more than once.";
        mopend(mk::line.substr(0, len_trim(mk::line)));
        k = index1(allkey, "=>");
        if (k > 0) {
          m = 0;
          ii = 1;
          for (l = 1; l <= 2; ++l) {
            n = index1(allkey.substr(ii - 1, k - ii), txt.substr(1, len_trim(txt)));
            if (n > 0) {
              m = m + 1;
              ii = ii + n;
            }
          }
          if (m == 1) {
            std::fprintf(stdout,
                         "%10s\n%10s\n%10s\n%2s%s\n%10s\n%2s%s\n%10s\n",
                         "This error was caused by option SALT within keyword SITE creating"
                         " a salt-bridge that involved residue " + txt.substr(1, len_trim(txt)) + ".",
                         "Original SITE keyword:",
                         "", allkey.substr(2, k - 3 - 2 + 1).c_str(),
                         "Expanded SITE keyword:",
                         "", allkey.substr(k + 3).c_str(),
                         "To correct this error, use two \"SITE\" keywords");
          } else {
            std::fprintf(stdout, "%10s\n", "Delete the faulty residue.");
          }
          return;
        }
      }
    }
  }

  changes.assign(mk::numat + 1, std::string(30, ' '));
  nadd = 0;
  ndel = 0;
  for (i = 1; i <= mk::natoms; ++i)
    for (int kk = 1; kk <= 3; ++kk) mc::geo[kk][i] = mc::coord[kk-1][i];
  i_charge = 100;
  for (i_atom = 1; i_atom <= mk::numat; ++i_atom) {
    l_res = false;
    if (nres > 0) {
      for (jres = 1; jres <= nres; ++jres) {
        if (fsub(mc::txtatm[i_atom], 22, 22)[0] != chain[jres]) continue;
        i = nint(reada(mc::txtatm[i_atom], 23));
        if (i == res[jres] && charge[jres][1] != '*') break;
      }
      l_res = (jres <= nres);
      if (l_res) {
        if (charge[jres][1] == '+') i_charge = 1;
        else if (charge[jres][1] == '-') i_charge = -1;
        else i_charge = 0;
      } else {
        i_charge = 100;  // "not l_res"
      }
    }

    // Deliberately add or remove hydrogen atoms as required by keyword SITE.
    if (neutral[5] || neutral[6] || l_res) {  // Must be -Arg(+)- or -Arg-
      // Check for Arg or Arg(+).
      if (mc::nat[i_atom] == 6 && mc::nbonds[i_atom] == 3) {
        if (mc::nat[mc::ibonds[1][i_atom]] == 7 &&
            mc::nat[mc::ibonds[2][i_atom]] == 7 &&
            mc::nat[mc::ibonds[3][i_atom]] == 7) {  // C bonded to three N
          i = 0;
          for (j = 1; j <= 3; ++j) {
            k = mc::ibonds[j][i_atom];
            for (l = 1; l <= mc::nbonds[k]; ++l) {
              m = mc::ibonds[l][k];
              if (mc::nat[m] == 1) i = i + 1;
            }
          }
          if (i_charge != 100 && (i == 4 || i == 5)) {
            charge[jres][1] = charge[jres][2];
            charge[jres][2] = charge[jres][3];
          }
          if (i != 5 && (neutral[6] || i_charge == 0) ||
              i != 4 && (neutral[5] || i_charge == 1)) continue;
          if (neutral[6] || i_charge == 0) {
            // Found the system -NH-C(NH2)2, now delete a hydrogen atom.
            for (j = 1; j <= 3; ++j) {
              if (mc::nbonds[mc::ibonds[j][i_atom]] == 3) {
                l = mc::ibonds[j][i_atom];
                k = 0;
                if (mc::nat[mc::ibonds[1][l]] == 1) k = k + 1;
                if (mc::nat[mc::ibonds[2][l]] == 1) k = k + 1;
                if (mc::nat[mc::ibonds[3][l]] == 1) k = k + 1;
                if (k == 2) {
                  if (mc::nat[mc::ibonds[1][l]] == 1) k = mc::ibonds[1][l];
                  else k = mc::ibonds[2][l];
                  mc::nat[k] = 99;
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
                  break;
                }
              }
            }
          } else if (neutral[5] || i_charge == 1) {
            // Found the system -NH-C(NH2)(NH), now add a hydrogen atom.
            for (j = 1; j <= 3; ++j) {
              if (mc::nbonds[mc::ibonds[j][i_atom]] == 2) {  // N with two bonds
                l = mc::ibonds[j][i_atom];
                add_sp2_H(mc::ibonds[1][l], l, mc::ibonds[2][l]);
                mk::numat = mk::numat + 1;
                mc::nbonds[i_atom] = mc::nbonds[i_atom] + 1;
                mc::ibonds[mc::nbonds[i_atom]][i_atom] = mk::numat;
                nadd = nadd + 1;
                setsub(changes[nadd], 1, 15, fsub(mc::txtatm[l], 12, 26));
                mc::txtatm[mk::natoms] =
                    fmt_txtatm_H(fsub(mc::txtatm[l], 1, 6), mk::natoms,
                                 fsub(mc::txtatm[l], 15, 26));
                break;
              }
            }
          }
        }
      }
    }
    if (neutral[1] || neutral[2] || l_res) {
      if (mc::nat[i_atom] == 6 &&
          (mc::nbonds[i_atom] == 3 || (mc::nbonds[i_atom] > 2 && l_res))) {
        // Check for -COOH or -COO(-) or -CH2OH or -CH2O(-).
        k = 0;
        l = 0;
        for (i = 1; i <= 3; ++i) {
          j = mc::ibonds[i][i_atom];
          if (mc::nat[j] == 8) {
            if (mc::nbonds[j] == 1) {
              k = k + 1;
            } else if (mc::nbonds[j] == 2 &&
                       (mc::nat[mc::ibonds[1][j]] == 1 ||
                        mc::nat[mc::ibonds[2][j]] == 1)) {
              k = k + 1;
            }
            l = j;
          }
        }
        // k is the number of oxygen atoms that are ionizable or are ionized.
        if (k > 0) {
          if (k == 2 || k == 1 && (mc::nbonds[i_atom] == 4 ||
              (l_res && mc::nbonds[i_atom] == 3 && distance(i_atom, l) > 1.3))) {
            if (i_charge != 100 && i_charge < 1) {
              charge[jres][1] = charge[jres][2];
              charge[jres][2] = charge[jres][3];
            }
            if (neutral[2] || i_charge == -1) {  // Make it -COO(-) or -CH2O(-)
              for (i = 1; i <= mc::nbonds[i_atom]; ++i) {
                j = mc::ibonds[i][i_atom];
                o[1] = mc::nat[mc::ibonds[1][j]];
                if (mc::nbonds[j] == 2) o[2] = mc::nat[mc::ibonds[2][j]];
                if (mc::nat[j] == 8 && mc::nbonds[j] == 2 &&
                    (o[1] == 1 || o[2] == 1)) {
                  if (mc::nat[mc::ibonds[1][j]] == 1) mc::nat[mc::ibonds[1][j]] = 99;
                  if (mc::nat[mc::ibonds[2][j]] == 1) mc::nat[mc::ibonds[2][j]] = 99;
                  mc::ibonds[mc::nbonds[j]][j] = 0;
                  mc::nbonds[j] = mc::nbonds[j] - 1;
                  break;
                }
              }
            } else if (neutral[1] || i_charge == 0) {
              k = 0;
              l = 0;
              for (i = 1; i <= mc::nbonds[i_atom]; ++i) {
                j = mc::ibonds[i][i_atom];
                if (mc::nat[j] == 8 && mc::nbonds[j] == 2) {
                  if (mc::nat[mc::ibonds[1][j]] == 1 || mc::nat[mc::ibonds[2][j]] == 1)
                    k = k + 1;
                }
                if (mc::nat[j] == 8 && mc::nbonds[j] == 1) {
                  l = l + 1;
                  o[l] = j;
                }
              }
              if (k == 0 && (l == 2 && mc::nbonds[i_atom] == 3 ||
                             l == 1 && mc::nbonds[i_atom] == 4)) {
                // -COO(-) or -CH2O(-) group identified.  Now add a hydrogen.
                add_sp_H(o[1], i_atom, o[2]);
                mc::nbonds[i_atom] = mc::nbonds[i_atom] + 1;
                mc::ibonds[mc::nbonds[i_atom]][i_atom] = mk::numat;
                nadd = nadd + 1;
                setsub(changes[nadd], 1, 15, fsub(mc::txtatm[j], 12, 26));
                mc::txtatm[mk::natoms] =
                    fmt_txtatm_H(fsub(mc::txtatm[j], 1, 6), mk::natoms,
                                 fsub(mc::txtatm[j], 15, 26));
              }
            }
          }
        }
      }
    }
    if (neutral[3] || neutral[4] || l_res) {  // Must be -NH3(+) or -NH2
      if (mc::nat[i_atom] == 7) {
        // Exclude -NH-C-(NH2) type.
        l = 0;
        for (j = 1; j <= mc::nbonds[i_atom]; ++j) {
          if (mc::nat[mc::ibonds[j][i_atom]] == 6) {
            i = mc::ibonds[j][i_atom];
            l = 0;
            for (k = 1; k <= mc::nbonds[i]; ++k) {
              if (mc::nat[mc::ibonds[k][i]] == 7) l = l + 1;
            }
            if (l == 3) break;
          }
        }
        if (l == 3) continue;
        i = 0;
        l = 0;
        for (j = 1; j <= mc::nbonds[i_atom]; ++j) {
          if (mc::nat[mc::ibonds[j][i_atom]] == 1) i = i + 1;
          for (k = 1; k <= mc::nbonds[mc::ibonds[j][i_atom]]; ++k) {
            m = mc::ibonds[k][mc::ibonds[j][i_atom]];
            if (mc::nat[m] == 7 && m != i_atom) l = 1;
          }
        }
        if (l == 1) continue;
        if ((i == 2 && mc::nbonds[i_atom] == 3 || i == 3 && mc::nbonds[i_atom] == 4) &&
            i_charge != 100) {
          charge[jres][1] = charge[jres][2];
          charge[jres][2] = charge[jres][3];
        }
        if (i == 2 && (neutral[3] || i_charge == 1) && mc::nbonds[i_atom] == 3) {
          // Selectively exclude -CO-NH2 if all -NH2 are flagged for ionization.
          if (neutral[3]) {
            for (i = 1; i <= mc::nbonds[i_atom]; ++i) {
              j = mc::ibonds[i][i_atom];
              if (mc::nat[j] == 6) break;
            }
            if (i <= mc::nbonds[i_atom]) {
              for (i = 1; i <= mc::nbonds[j]; ++i) {
                if (mc::nat[mc::ibonds[i][j]] == 8) break;
              }
              if (i <= mc::nbonds[j]) continue;
            }
          }
          // Nitrogen bonded to two hydrogen atoms, add a hydrogen atom.
          add_sp3_H(mc::ibonds[1][i_atom], i_atom, mc::ibonds[2][i_atom],
                    mc::ibonds[3][i_atom]);
          mc::nbonds[i_atom] = mc::nbonds[i_atom] + 1;
          mc::ibonds[mc::nbonds[i_atom]][i_atom] = mk::numat;
          mk::numat = mk::numat + 1;
          nadd = nadd + 1;
          setsub(changes[nadd], 1, 15, fsub(mc::txtatm[i_atom], 12, 26));
          mc::txtatm[mk::natoms] =
              fmt_txtatm_H(fsub(mc::txtatm[i_atom], 1, 6), mk::natoms,
                           fsub(mc::txtatm[i_atom], 15, 26));
        } else if (i == 3 && (neutral[4] || i_charge == 0) && mc::nbonds[i_atom] == 4) {
          // Nitrogen bonded to three hydrogen atoms, delete a hydrogen atom.
          for (i = 1; i <= 4; ++i) {
            if (mc::nat[mc::ibonds[i][i_atom]] == 1) {
              j = mc::ibonds[i][i_atom];
              mc::nat[j] = 99;
              mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
              mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
              break;
            }
          }
        }
      }
    }
    if (neutral[7] || neutral[8] || l_res) {  // Must be -His(+)- or -His-
      // Check for His or His(+).
      if (mc::nat[i_atom] == 6 && mc::nbonds[i_atom] == 3) {
        j = 0;
        k = 0;
        for (i = 1; i <= 3; ++i) {
          if (mc::nat[mc::ibonds[i][i_atom]] == 7) j = j + 1;
          if (mc::nat[mc::ibonds[i][i_atom]] == 1) k = k + 1;
        }
        if (j == 2 && k == 1) {  // C bonded to two N and one H, therefore His
          i = 0;
          for (j = 1; j <= 3; ++j) {
            k = mc::ibonds[j][i_atom];
            for (l = 1; l <= mc::nbonds[k]; ++l) {
              m = mc::ibonds[l][k];
              if (mc::nat[m] == 1) i = i + 1;
            }
          }
          if (i != 2 && (neutral[8] || i_charge == 0) ||
              i != 1 && (neutral[7] || i_charge == 1)) continue;
          if (i_charge != 100) {
            charge[jres][1] = charge[jres][2];
            charge[jres][2] = charge[jres][3];
          }
          if (neutral[8] || i_charge == 0) {
            // Found a His(+), now delete a hydrogen atom.
            bool outer_done = false;
            for (j = 1; j <= 3 && !outer_done; ++j) {
              k = mc::ibonds[j][i_atom];
              for (l = 1; l <= mc::nbonds[k]; ++l) {
                m = mc::ibonds[l][k];
                if (mc::nat[m] == 1) {
                  mc::nat[m] = 99;
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
                  outer_done = true;
                  break;
                }
              }
            }
          } else if (neutral[7] || i_charge == 1) {
            // Found the system -C-N-C-, now add a hydrogen atom.
            for (j = 1; j <= 3; ++j) {
              if (mc::nbonds[mc::ibonds[j][i_atom]] == 2) {  // N with two bonds
                l = mc::ibonds[j][i_atom];
                add_sp2_H(mc::ibonds[1][l], l, mc::ibonds[2][l]);
                mk::numat = mk::numat + 1;
                mc::nbonds[l] = mc::nbonds[l] + 1;
                mc::ibonds[mc::nbonds[l]][l] = mk::numat;
                nadd = nadd + 1;
                setsub(changes[nadd], 1, 15, fsub(mc::txtatm[l], 12, 26));
                mc::txtatm[mk::natoms] =
                    fmt_txtatm_H(fsub(mc::txtatm[l], 1, 6), mk::natoms,
                                 fsub(mc::txtatm[l], 15, 26));
                break;
              }
            }
          }
        }
      }
    }
    if (neutral[9]) {  // Check for sulfate.  If S-O-H exists, delete the "H".
      if (mc::nat[i_atom] == 8 && mc::nbonds[i_atom] == 2) {
        for (j = 1; j <= 2; ++j) {
          if (mc::nat[mc::ibonds[j][i_atom]] == 16) break;  // S attached to the O
        }
        if (j < 3) {
          for (k = 1; k <= 2; ++k) {
            if (mc::nat[mc::ibonds[k][i_atom]] == 1) break;  // H attached to the O
          }
          if (k < 3) {
            l = mc::ibonds[j][i_atom];
            m = 0;
            for (n = 1; n <= mc::nbonds[l]; ++n) {
              if (mc::nat[mc::ibonds[n][l]] == 8) m = m + 1;
            }
            if (m == 4) {
              mc::nat[mc::ibonds[k][i_atom]] = 99;
              mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
              mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
            }
          }
        }
      }
    }
    if (neutral[10]) {
      // Check for phosphate.  Aim for [PO4](-), e.g., CH3-PO4-H.
      if (mc::nat[i_atom] == 8 && mc::nbonds[i_atom] == 2) {
        for (j = 1; j <= 2; ++j) {
          if (mc::nat[mc::ibonds[j][i_atom]] == 15) break;  // P attached to the O
        }
        if (j < 3) {
          if (mc::nbonds[mc::ibonds[j][i_atom]] == 4) {
            // Oxygen atoms on phosphorus should be attached to only two ligands.
            m = mc::ibonds[j][i_atom];
            l = 0;
            for (k = 1; k <= 4; ++k) {
              if (mc::nat[mc::ibonds[k][m]] == 8 && mc::nbonds[mc::ibonds[k][m]] == 2) {
                if (mc::nat[mc::ibonds[1][mc::ibonds[k][m]]] != 99 &&
                    mc::nat[mc::ibonds[2][mc::ibonds[k][m]]] != 99) l = l + 1;
              }
            }
            if (l > 2) {
              for (k = 1; k <= 2; ++k) {
                if (mc::nat[mc::ibonds[k][i_atom]] == 1) break;  // H attached to the O
              }
              if (k < 3) {
                l = mc::ibonds[j][i_atom];
                m = 0;
                for (n = 1; n <= mc::nbonds[l]; ++n) {
                  if (mc::nat[mc::ibonds[n][l]] == 8) m = m + 1;
                }
                if (m == 4) {
                  mc::nat[mc::ibonds[k][i_atom]] = 99;
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
                }
              }
            }
          }
        }
      }
    }
    if (l_res) {
      if (charge[jres][1] != '*') {
        // Check for sulfate.
        if (mc::nat[i_atom] == 16 && mc::nbonds[i_atom] == 4) {
          k = 2;
          for (i = 1; i <= 4; ++i) {
            j = mc::ibonds[i][i_atom];
            if (mc::nat[j] != 8) break;
            if (mc::nbonds[j] == 1) k = k - 1;
          }
          if (i == 5) {
            // Found a sulfate.
            i_charge = (i_charge < 0) ? i_charge : 0;
            if (charge[jres][2] == '-') i_charge = -2;
            // Work out the change in the number of ionized atoms.
            i_charge = i_charge - k;
            if (i_charge < 0) {
              // Now delete hydrogens, as necessary.
              for (i = 1; i <= -i_charge; ++i) {
                for (l = 1; l <= 4; ++l) {
                  j = mc::ibonds[l][i_atom];
                  if (mc::nbonds[j] == 2) {
                    if (mc::nat[mc::ibonds[2][j]] == 1) break;
                  }
                }
                if (l < 5) {
                  mc::nat[mc::ibonds[2][j]] = 99;
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
                }
              }
            } else if (i_charge > 0) {
              // Now add hydrogens, as necessary.
              l = 1;
              for (i = 1; i <= 4; ++i) {
                j = mc::ibonds[i][i_atom];
                if (mc::nbonds[j] == 1) {
                  for (k = 1; k <= 4; ++k) {
                    if (k != i) break;
                  }
                  add_sp_H(j, i_atom, mc::ibonds[k][i_atom]);
                  mk::numat = mk::numat + 1;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] + 1;
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = mk::numat;
                  nadd = nadd + 1;
                  setsub(changes[nadd], 1, 15, fsub(mc::txtatm[j], 12, 26));
                  l = l + 1;
                  if (l > i_charge) break;
                }
              }
            }
            charge[jres][1] = '*';
          }
        }
        // Check for phosphate.
        if (mc::nat[i_atom] == 15 && mc::nbonds[i_atom] == 4) {
          k = 1;
          for (i = 1; i <= 4; ++i) {
            j = mc::ibonds[i][i_atom];
            if (mc::nat[j] != 8) break;
            if (mc::nbonds[j] == 1) k = k - 1;
          }
          if (i == 5) {
            // Found a phosphate.
            i_charge = (i_charge < 0) ? i_charge : 0;
            if (charge[jres][2] == '-') i_charge = -2;
            i_charge = i_charge - k;
            if (i_charge < 0) {
              // Now delete hydrogens, as necessary.
              for (i = 1; i <= -i_charge; ++i) {
                for (l = 1; l <= 4; ++l) {
                  j = mc::ibonds[l][i_atom];
                  if (mc::nbonds[j] == 2) {
                    if (mc::nat[mc::ibonds[2][j]] == 1) break;
                  }
                }
                if (l < 5) {
                  mc::nat[mc::ibonds[2][j]] = 99;
                  charge[jres][1] = charge[jres][2];
                  charge[jres][2] = charge[jres][3];
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = 0;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] - 1;
                }
              }
            } else if (i_charge > 0) {
              // Now add hydrogens, as necessary.
              l = 1;
              for (i = 1; i <= 4; ++i) {
                j = mc::ibonds[i][i_atom];
                if (mc::nbonds[j] == 1) {
                  for (k = 1; k <= 4; ++k) {
                    if (k != i) break;
                  }
                  add_sp_H(j, i_atom, mc::ibonds[k][i_atom]);
                  mk::numat = mk::numat + 1;
                  nadd = nadd + 1;
                  mc::nbonds[i_atom] = mc::nbonds[i_atom] + 1;
                  mc::ibonds[mc::nbonds[i_atom]][i_atom] = mk::numat;
                  setsub(changes[nadd], 1, 15, fsub(mc::txtatm[j], 12, 26));
                  mc::txtatm[mk::natoms] =
                      fmt_txtatm_H(fsub(mc::txtatm[j], 1, 6), mk::natoms,
                                   fsub(mc::txtatm[j], 15, 26));
                  charge[jres][1] = charge[jres][2];
                  charge[jres][2] = charge[jres][3];
                  l = l + 1;
                  if (l > i_charge) break;
                }
              }
            } else {
              charge[jres][1] = charge[jres][2];
              charge[jres][2] = charge[jres][3];
            }
          }
        }
      }
    }
  }

  // Now check for explicit atoms.
  i = index1(mk::keywrd, " SITE=");
  j = index1(mk::keywrd.substr(i - 1), ") ") + i;
  while (true) {
    k = index1(mk::keywrd.substr(i - 1), ") ") + i;
    j = index1(mk::keywrd.substr(i - 1, k - i), "\"") + i;
    if (j == i) break;
    l = index1(mk::keywrd.substr(j - 1), "\"") + j;
    mk::line = mk::keywrd.substr(j - 1, l - 1 - j);
    if (fsub(mk::line, 1, 1) == "[") {
      // Atom defined using Jmol format.
      n = index1(mk::keywrd.substr(j - 1, l - 1 - j), ".") + j;
      k = index1(mk::keywrd.substr(j - 1, l - 1 - j), ":") + j;
      m = index1(mk::keywrd.substr(j - 1, l - 1 - j), "]") + j;
      mk::line = mk::keywrd.substr(n - 1, l - 1 - n + 1) +
                 mk::keywrd.substr(j, m - 2 - j + 1) +
                 mk::keywrd.substr(k - 1, 1) +
                 mk::keywrd.substr(m - 1, k - 2 - m + 1);
    }
    mk::keywrd.replace(j - 1, mk::keywrd.size() - (j - 1),
                       mk::line.substr(0, len_trim(mk::line)) +
                       mk::keywrd.substr(l - 2));
    i = j + (int)len_trim(mk::line) + 2;
  }
  i = index1(mk::keywrd, " SITE=");
  j = index1(mk::keywrd.substr(i - 1), ") ") + i;
  mk::line = mk::keywrd.substr(i + 6 - 1, j - (i + 6));
  while (true) {
    i = index1(mk::line, "\"");
    if (i == 0) break;
    for (j = i + 1; j <= (int)len_trim(mk::line); ++j) {
      if (fsub(mk::line, j, j) == "\"") break;
    }
    std::string res_txt = mk::line.substr(i, j - 1 - i);
    if (fsub(mk::line, j + 2, j + 2) == "+") i_charge = 1;
    if (fsub(mk::line, j + 2, j + 2) == "0") i_charge = 0;
    if (fsub(mk::line, j + 2, j + 2) == "-") i_charge = -1;
    mk::line = mk::line.substr(j);
    m = 0;
    for (k = 1; k <= (int)len_trim(res_txt); ++k) {
      if (fsub(res_txt, k, k) != " ") {
        m = m + 1;
        setsub(res_txt, m, m, fsub(res_txt, k, k));
      }
    }
    res_txt.resize(m + 1, ' ');
    for (i = 1; i <= mk::numat; ++i) {
      std::string txt = mc::txtatm[i];
      n = 0;
      for (k = 13; k <= 26; ++k) {
        if (fsub(txt, k, k) != " ") {
          n = n + 1;
          setsub(txt, n, n, fsub(txt, k, k));
        }
      }
      txt.resize(n + 1, ' ');
      n = (n > m) ? n : m;
      for (k = 1; k <= n; ++k) {
        if (fsub(res_txt, k, k) != fsub(txt, k, k) &&
            fsub(res_txt, k, k) != "*") break;
      }
      if (k > n) break;
    }
    j = mc::nat[i];
    switch (j) {
      case 7:  // Select for Nitrogen in different coordination numbers.
        switch (mc::nbonds[i]) {
          case 2:
            add_sp2_H(mc::ibonds[1][i], i, mc::ibonds[2][i]);
            mk::numat = mk::numat + 1;
            mc::nbonds[i] = mc::nbonds[i] + 1;
            mc::ibonds[mc::nbonds[i]][i] = mk::numat;
            nadd = nadd + 1;
            setsub(changes[nadd], 1, 15, fsub(mc::txtatm[i], 12, 26));
            mc::txtatm[mk::natoms] =
                fmt_txtatm_H(fsub(mc::txtatm[i], 1, 6), mk::natoms,
                             fsub(mc::txtatm[i], 15, 26));
            break;
          case 3:
            if (i_charge == 1) {
              // Three-coordinate nitrogen, make into a quaternary nitrogen.
              add_a_sp3_hydrogen_atom_ext(i, mc::ibonds[1][i], mc::ibonds[2][i],
                                          mc::ibonds[3][i], 1.1, metals, nmetals);
              mk::natoms = mk::natoms + 1;
              mc::nbonds[i] = mc::nbonds[i] + 1;
              mc::ibonds[mc::nbonds[i]][i] = mk::numat;
              for (int kk = 1; kk <= 3; ++kk) mc::geo[kk][mk::numat] = mc::coord[kk-1][mk::numat];
              nadd = nadd + 1;
              setsub(changes[nadd], 1, 15, fsub(mc::txtatm[i], 12, 26));
              mc::txtatm[mk::numat] =
                  fmt_txtatm_H(fsub(mc::txtatm[i], 1, 6), mk::natoms,
                               fsub(mc::txtatm[i], 15, 26));
            } else if (i_charge == 0) {
              // Three-coordinate nitrogen, neutralize (assume pyridinium).
              for (j = 1; j <= mc::nbonds[i]; ++j) {
                if (mc::nat[mc::ibonds[j][i]] == 1) break;
              }
              if (j > mc::nbonds[i]) {
                std::fprintf(stdout, "\n %s\n",
                             " A nitrogen atom defined by keyword SITE is not bonded to any hydrogen atoms.");
                std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                             "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                             0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                             "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                bug = true;
              }
              mc::nat[mc::ibonds[j][i]] = 99;
              mc::ibonds[mc::nbonds[i]][i] = 0;
              mc::nbonds[i] = mc::nbonds[i] - 1;
            } else {
              std::fprintf(stdout, "\n %s\n",
                           " A nitrogen atom defined by keyword SITE cannot have a negative charge.");
              std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                           "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                           0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                           "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
              bug = true;
            }
            break;
          case 4:
            if (i_charge == 0) {
              // Quaternary nitrogen, neutralize.
              for (j = 1; j <= 4; ++j) {
                if (mc::nat[mc::ibonds[j][i]] == 1) break;
              }
              if (j > 4) {
                std::fprintf(stdout, "\n %s\n",
                             " A nitrogen atom defined by keyword SITE is not bonded to any hydrogen atoms.");
                std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                             "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                             0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                             "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                bug = true;
              }
              mc::nat[mc::ibonds[j][i]] = 99;
              mc::ibonds[mc::nbonds[i]][i] = 0;
              mc::nbonds[i] = mc::nbonds[i] - 1;
            } else if (i_charge == 1) {
              std::fprintf(stdout, "\n %s\n",
                           " A nitrogen atom defined by keyword SITE is already bonded to four atoms.");
              std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                           "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                           0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                           "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
            }
            break;
        }
        break;
      case 8:  // Select for Oxygen in different coordination numbers.
        // label 99: retry after deleting an H on a 3-coordinate oxygen.
        {
        bool goto99 = false;
        do {
          goto99 = false;
          switch (mc::nbonds[i]) {
            case 1:
              if (i_charge == 0 || i_charge == 1) {
                j = mc::ibonds[1][i];
                for (k = 1; k <= mc::nbonds[j]; ++k) {  // search for -COO structure
                  if (mc::ibonds[k][j] != i && mc::nat[mc::ibonds[k][j]] == 8) break;
                }
                if (k > mc::nbonds[j]) {
                  for (k = 1; k <= mc::nbonds[j]; ++k) {  // search for anything other than atom i
                    if (mc::ibonds[k][j] != i) break;
                  }
                }
                k = mc::ibonds[k][j];
                add_sp_H(i, j, k);
                mk::numat = mk::numat + 1;
                nadd = nadd + 1;
                setsub(changes[nadd], 1, 15, fsub(mc::txtatm[i], 12, 26));
                if (i_charge == 0) {
                  mc::nbonds[i] = 2;
                  mc::ibonds[2][i] = mk::natoms;
                }
              } else if (i_charge == -1) {
                std::fprintf(stdout, "\n %s\n",
                             " An oxygen atom defined by keyword SITE is already only bonded to one atom.");
                std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                             "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                             0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                             "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                bug = true;
              }
              break;
            case 2:
              if (i_charge == -1) {
                for (j = 1; j <= 2; ++j) {
                  if (mc::nat[mc::ibonds[j][i]] == 1) break;
                }
                if (j > 2) {
                  std::fprintf(stdout, "\n %s\n",
                               " An oxygen atom defined by keyword SITE is not bonded to any hydrogen atoms.");
                  std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                               "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                               0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                               "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                  bug = true;
                } else {
                  mc::nat[mc::ibonds[j][i]] = 99;
                  mc::ibonds[mc::nbonds[i]][i] = 0;
                  mc::nbonds[i] = mc::nbonds[i] - 1;
                }
              } else if (i_charge == 0) {
                std::fprintf(stdout, "\n %s\n",
                             " An oxygen atom defined by keyword SITE is already bonded to two atoms.");
                std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                             "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                             0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                             "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                bug = true;
              } else {
                add_sp2_H(mc::ibonds[1][i], i, mc::ibonds[2][i]);
                mk::numat = mk::numat + 1;
                mc::nbonds[i] = mc::nbonds[i] + 1;
                mc::ibonds[mc::nbonds[i]][i] = mk::numat;
              }
              break;
            case 3:
              if (i_charge == -1 || i_charge == 0) {
                for (j = 1; j <= mc::nbonds[i]; ++j) {
                  if (mc::nat[mc::ibonds[j][i]] == 1) break;
                }
                if (j > mc::nbonds[i]) {
                  std::fprintf(stdout, "\n %s\n",
                               " An oxygen atom defined by keyword SITE is not bonded to any hydrogen atoms.");
                  std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                               "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                               0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                               "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                  bug = true;
                } else {
                  mc::nat[mc::ibonds[j][i]] = 99;
                  mc::nbonds[i] = 2;
                  if (j == 2) mc::ibonds[2][i] = mc::ibonds[3][i];
                  goto99 = true;  // goto 99: re-test with nbonds=2
                }
              } else {
                std::fprintf(stdout, "\n %s\n",
                             " An oxygen atom defined by keyword SITE is already bonded to three atoms.");
                std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                             "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                             0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                             "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
                bug = true;
              }
              break;
          }
        } while (goto99);
        }
        break;
      default:
        std::fprintf(stdout, "\n %s\n",
                     " An atom defined by keyword SITE is not a nitrogen or oxygen.");
        std::fprintf(stdout, " %s%4d%8.3f%8.3f%8.3f %s%s%s\n",
                     "  Faulty atom = \"" + mc::txtatm[i] + "\"",
                     0, mc::coord[0][i], mc::coord[1][i], mc::coord[2][i],
                     "  1.00  0.00      PROT", ec::elemnt[mc::nat[i]].c_str(), "\"");
        bug = true;
        break;
    }
  }

  // Check that all sites have been recognized, whether they are modified or not.
  k = 0;
  for (j = 1; j <= nres; ++j) {
    if (charge[j][1] != '*') {
      if (k == 0) {
        k = index1(mk::keywrd, " SITE=");
        l = index1(mk::keywrd.substr(k - 1), ") ") + k;
        if (!is_blank(allkey))
          std::fprintf(stdout, "\n %s\n",
                       ("Faulty SITE Keyword: '" + allkey.substr(0, len_trim(allkey)) + "'").c_str());
      }
      k = index1(mk::keywrd, " SITE=");
      if (j == 1) {
        for (k = k; k <= (int)len_trim(mk::keywrd); ++k) {
          if (fsub(mk::keywrd, k, k) == "(") break;
        }
        k = k + 1;
      } else {
        for (i = 1; i <= j - 1; ++i) {
          k = index1(mk::keywrd.substr(k, l - k - 1), ",") + k + 1;
        }
      }
      i = index1(mk::keywrd.substr(k, l - k - 1), ",");
      if (i == 0) i = index1(mk::keywrd.substr(k, l - k), ") ");
      i = i + k - 1;
      char b[128];
      std::snprintf(b, sizeof b, " Site:%3d, '%s', was not properly recognized.", j,
                    mk::keywrd.substr(k - 1, i - k + 1).c_str());
      mopend(std::string(b));
      std::fprintf(stdout, "\n%10s\n",
                   "(Check that the residue is not duplicated and can be ionized or de-ionized)");
      if (mk::keywrd.find(" 0SCF") != std::string::npos) mk::moperr = false;
    }
  }
  if (bug) {
    if (!is_blank(allkey)) {
      mopend("FAULTY SITE KEYWORD: '" + allkey.substr(0, len_trim(allkey)) + "'");
    } else {
      mopend("FAULTY SITE KEYWORD");
    }
    std::fprintf(stdout, "\n%10s\n%10s\n",
                 "(Check that the atom is N or O and can be ionized or de-ionized)", "");
    if (mk::keywrd.find(" 0SCF") != std::string::npos) mk::moperr = false;
  }
  j = 0;
  for (i = 1; i <= mk::natoms - mk::id; ++i) {
    if (mc::nat[i] != 99) {
      j = j + 1;
      for (int kk = 1; kk <= 3; ++kk) mc::geo[kk][j] = mc::geo[kk][i];
      for (int kk = 1; kk <= 3; ++kk) mc::coord[kk-1][j] = mc::geo[kk][i];
      mc::nat[j] = mc::nat[i];
      mc::atmass[j] = pc::ams[mc::nat[i]];
      mc::txtatm[j] = mc::txtatm[i];
    } else {
      ndel = ndel + 1;
      setsub(changes[ndel], 16, 30, fsub(mc::txtatm[i], 12, 26));
    }
  }
  mk::numat = j;
  mk::natoms = mk::numat;
  mc::labels.resize(mk::numat + 1);
  for (i = 1; i <= mk::numat; ++i) mc::labels[i] = mc::nat[i];
  set_up_dentate();
  update_txtatm(true, false);
  j = (nadd > ndel) ? nadd : ndel;
  if (j > 0) {
    std::fprintf(stdout, "\n%4s%s\n", "", "Changes in Ionization caused by keyword SITES");
    std::fprintf(stdout, "%4s%s\n", "", ("Keyword:" + allkey.substr(0, len_trim(allkey))).c_str());
    std::fprintf(stdout, "\n%17s\n%11s%s\n", "Hydrogen Atoms",
                 "", "Added                Deleted");
    for (i = 1; i <= j; ++i) {
      std::fprintf(stdout, "%5s%-15s%s\n", "",
                   fsub(changes[i], 1, 15).c_str(),
                   fsub(changes[i], 16, 30).c_str());
    }
  }
}
