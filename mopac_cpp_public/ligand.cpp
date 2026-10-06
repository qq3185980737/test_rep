// ligand.cpp — C++ translation of "ligand.F90" (43053 B).
// Contains ligand(), moiety(), nheavy(), identify_hexose(), inc_res().
#include "ligand.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "dihed.h"
#include "distance.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "MOZYME_C.h"
#include "reada.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace MOZYME_C;

namespace {

// 1-based Fortran Index.
int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// 1-based inclusive slice.
std::string sl(const std::string& s, int f, int l) {
  if (f < 1) f = 1;
  if (l > (int)s.size()) l = (int)s.size();
  if (f > l) return "";
  return s.substr(f - 1, l - f + 1);
}

int nint(double x) { return static_cast<int>(std::floor(x + 0.5)); }

std::size_t len_trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? 0 : e + 1;
}

std::string trim(const std::string& s) {
  std::size_t b = s.find_first_not_of(' ');
  if (b == std::string::npos) return "";
  std::size_t e = s.find_last_not_of(' ');
  return s.substr(b, e - b + 1);
}

// PDB-style record: 6-char tag, 5-digit serial, element (3), residue (5),
// residue number (6) -> exactly 26 chars.
std::string hetatm26(int i, const std::string& el, const std::string& res,
                     int ires) {
  char b[64];
  std::snprintf(b, sizeof b, "%-6s%5d %-3s%-5s%6d", "HETATM", i, el.c_str(),
                res.c_str(), ires);
  std::string s = b;
  if ((int)s.size() < 26) s += std::string(26 - s.size(), ' ');
  return s;
}

// a6,i5,a3,a6,i6 -> no separating blank.
std::string atom26(const std::string& tag, int i, const std::string& el,
                   const std::string& res, int resn) {
  char b[64];
  std::snprintf(b, sizeof b, "%-6s%5d%-3s%-6s%6d", tag.c_str(), i, el.c_str(),
                res.c_str(), resn);
  std::string s = b;
  if ((int)s.size() < 26) s += std::string(26 - s.size(), ' ');
  return s;
}

// a6,i5,a3,a6 (no residue number).
std::string atom26nr(const std::string& tag, int i, const std::string& el,
                     const std::string& res) {
  char b[64];
  std::snprintf(b, sizeof b, "%-6s%5d%-3s%-6s", tag.c_str(), i, el.c_str(),
                res.c_str());
  std::string s = b;
  if ((int)s.size() < 26) s += std::string(26 - s.size(), ' ');
  return s;
}

// Element label for PDB: 2-char symbol, blank-padded to 3 with leading space
// if the symbol is one letter (F90: el = cap_elemnt; if el(2:2)==" " el =
// " "//el(1:1)).
std::string pdb_el(int z) {
  std::string el = cap_elemnt[z];
  if (el.size() >= 2 && el[1] == ' ') el = " " + el.substr(0, 1);
  return el;
}

bool het_elem(int z) {  // select case (6:9, 15:17, 34:35, 53)
  return (z >= 6 && z <= 9) || (z >= 15 && z <= 17) || (z >= 34 && z <= 35) ||
         z == 53;
}

}  // namespace

// ---------------------------------------------------------------------------
void inc_res(int& ires, const std::vector<int>& start_res, int& nfrag) {
  if (start_res[nfrag] != -200) {
    ires = start_res[nfrag];
    nfrag = nfrag + 1;
  }
  ires = ires + 1;
}

// ---------------------------------------------------------------------------
int nheavy(int icc) {
  int n = 0;
  for (int i = 1; i <= nbonds[icc]; ++i)
    if (nat[ibonds[i][icc]] > 1) n = n + 1;
  return n;
}

// ---------------------------------------------------------------------------
void ligand(int& ires, const std::vector<int>& start_res, int& nfrag) {
  const int maxlive = 2000;
  int i, j, k, l, m, n, ii, jj, kk, mm, nn, res, change_no, max_comments;
  std::vector<int> live(maxlive + 1, 0), n_ele(107, 0), inres;
  int ninres, i_panic, nlive, n_C, n_O, n_N, n_S, n_P;
  std::vector<bool> l_used(natoms + 1, false);
  bool attached, l_chain, first;
  std::string het(3, ' '), het_group(80, ' '), num(1, ' '), num1(1, ' '),
      el, het2(3, ' ');

  j = 0;
  for (i = 1; i <= ncomments; ++i) {
    if (index1(all_comments[i], "REMARK   3") == 0) {
      j = j + 1;
      all_comments[j] = all_comments[i];
    }
  }
  ncomments = j;
  max_comments = (int)all_comments.size();
  first = true;
  change_no = 0;

  for (i = 1; i <= natoms - id; ++i) {
    if (sl(txtatm[i], 26, 26) != " ") continue;
    //
    //  Residue number is missing, therefore not a residue.
    //
    //  Test for specific molecules
    //
    if (labels[i] == 15) {
      if (nbonds[i] == 4) {
        bool all_oxygen = true;
        for (j = 1; j <= 4; ++j)
          if (labels[ibonds[j][i]] != 8) { all_oxygen = false; break; }
        if (all_oxygen) {
          // Phosphate
          inc_res(ires, start_res, nfrag);
          txtatm[i] = hetatm26(i, " P ", "PO4", ires);
          for (j = 1; j <= 4; ++j)
            txtatm[ibonds[j][i]] = hetatm26(ibonds[j][i], " O ", "PO4", ires);
        }
      }
    } else if (labels[i] == 16) {
      if (nbonds[i] == 4) {
        bool all_oxygen = true;
        for (j = 1; j <= 4; ++j)
          if (labels[ibonds[j][i]] != 8) { all_oxygen = false; break; }
        if (all_oxygen) {
          // Sulfate
          inc_res(ires, start_res, nfrag);
          txtatm[i] = hetatm26(i, " S ", "SO4", ires);
          for (j = 1; j <= 4; ++j)
            txtatm[ibonds[j][i]] = hetatm26(ibonds[j][i], " O ", "SO4", ires);
        }
      }
    } else if (labels[i] == 8 && nheavy(i) == 0) {
      // Water
      inc_res(ires, start_res, nfrag);
      txtatm[i] = hetatm26(i, " O ", "HOH", ires);
      for (j = 1; j <= nbonds[i]; ++j)
        txtatm[ibonds[j][i]] = hetatm26(ibonds[j][i], " H ", "HOH", ires);
    } else if (labels[i] == 8) {
      // Ethylene glycol and glycerol
      // H-O-CH2-CH2-O-H and HOCH2-HCOH-H2COH
      if (nheavy(i) == 1) {
        for (ii = 1; ii <= nbonds[i]; ++ii) {
          j = ibonds[ii][i];
          if (nat[j] != 6) continue;
          if (nheavy(j) != 2) continue;  // Oxygen joined to carbon(1)
          for (jj = 1; jj <= nbonds[j]; ++jj) {
            k = ibonds[jj][j];
            if (nat[k] != 6) continue;
            if (nheavy(k) == 1) {
              // Found methanol
            } else if (nheavy(k) == 2) {
              // Carbon(1) joined to oxygen and carbon(2)
              for (kk = 1; kk <= nbonds[k]; ++kk) {
                l = ibonds[kk][k];
                if (nat[l] == 8 && l != j) l = l;
                if (nat[l] != 8) continue;
                if (nheavy(l) == 1) {
                  // Found O-C-C-O
                  inc_res(ires, start_res, nfrag);
                  txtatm[i] = hetatm26(i, " O ", "EDO", ires);
                  txtatm[j] = hetatm26(i, " C ", "EDO", ires);
                  txtatm[k] = hetatm26(i, " C ", "EDO", ires);
                  txtatm[l] = hetatm26(i, " O ", "EDO", ires);
                }
              }
            } else if (nheavy(k) == 3) {
              for (kk = 1; kk <= nbonds[k]; ++kk) {
                l = ibonds[kk][k];
                if (nat[l] != 8) continue;
                if (nheavy(l) != 1) continue;
                for (mm = 1; mm <= nbonds[k]; ++mm) {
                  m = ibonds[mm][k];
                  if (nat[m] != 6) continue;
                  if (m == j) continue;
                  if (nheavy(m) != 2) continue;  // Carbon(3)
                  for (nn = 1; nn <= nbonds[m]; ++nn) {
                    n = ibonds[nn][m];
                    if (nat[n] != 8) continue;
                    // Found O-C-CO-C-O
                    inc_res(ires, start_res, nfrag);
                    txtatm[i] = hetatm26(i, " O ", "GOL", ires);
                    txtatm[j] = hetatm26(i, " C ", "GOL", ires);
                    txtatm[k] = hetatm26(i, " C ", "GOL", ires);
                    txtatm[l] = hetatm26(i, " O ", "GOL", ires);
                    txtatm[m] = hetatm26(i, " C ", "GOL", ires);
                    txtatm[n] = hetatm26(i, " O ", "GOL", ires);
                  }
                }
              }
            }
          }
        }
      }
    }
    // label 1000: common tail for the "specific molecule" tests
    n_ele.assign(107, 0);
    ninres = 0;
    i_panic = 0;
    l_used.assign(natoms + 1, false);
    if (sl(txtatm[i], 26, 26) == " " && nat[i] != 1) {
      ninres = 1;
      inres.assign(natoms + 1, 0);
      inres[1] = i;
      l_used[i] = true;
      nlive = nbonds[i];
      for (int q = 1; q <= nlive; ++q) live[q] = ibonds[q][i];
      if (!het_elem(nat[i])) nlive = 0;
      inc_res(ires, start_res, nfrag);
      while (true) {  // outer_loop
        l = live[1];
        if (nlive == 0) break;
        if (l == 0) break;
        if (l_used[l] || sl(txtatm[l], 11, 14) != "    ") {
          //  Atom L is already labeled
          live[1] = live[nlive];
          nlive = nlive - 1;
        } else {
          //  Assign atom L to the hetero group.
          l_used[l] = true;
          for (ii = 1; ii <= nbonds[l]; ++ii) {
            jj = ibonds[ii][l];
            if (sl(txtatm[jj], 11, 14) != "    ") {
              if (het_elem(nat[i])) i_panic = jj;
            }
          }
          if (het_elem(nat[i])) {
            ninres = ninres + 1;
            inres[ninres] = l;
          }
          if (nbonds[l] != 0) {
            //  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
            for (ii = 2; ii <= nbonds[l]; ++ii) {
              j = ibonds[ii][l];
              bool dup = false;
              for (jj = 1; jj <= nlive; ++jj)
                if (live[jj] == j) { dup = true; break; }
              if (dup) continue;
              if (!l_used[j] && sl(txtatm[j], 11, 14) == "    ") {
                nlive = nlive + 1;
                if (nlive > maxlive) return;
                live[nlive] = j;
              }
            }
            live[1] = ibonds[1][l];
          } else {
            if (nlive == 0) break;
            live[1] = live[nlive];
            nlive = nlive - 1;
          }
        }
      }
      if (i_panic == 0) {
        n_ele.assign(107, 0);
        for (j = 1; j <= ninres; ++j) {
          k = inres[j];
          n_ele[nat[k]] = n_ele[nat[k]] + 1;
        }
        n_C = n_ele[6];
        n_N = n_ele[7];
        n_O = n_ele[8];
        n_S = n_ele[16];
        n_P = n_ele[15];
        attached = false;
        res = -10000;
        n = 0;
        for (k = 1; k <= ninres; ++k) {  // loop_k
          m = inres[k];
          for (l = 1; l <= nbonds[m]; ++l) {
            n = ibonds[l][m];
            if (sl(txtatm[n], 26, 26) != " ") {
              attached = true;
              res = nint(reada(txtatm[n], 23));
              goto loop_k_done;
            }
          }
        }
      loop_k_done:
        if (n > 0 && res != -10000) {
          if (sl(txtatm[n], 13, 16) == " N  " && nat[m] == 6) {
            //  Inserting a fragment at the start of a chain - so move all
            //  residue numbers in the rest of the chain up by 1
            k = res;
            while (true) {
              j = nint(reada(txtatm[n], 23));
              if (j == k) {
                char b[8];
                std::snprintf(b, sizeof b, "%4d", k + 1);
                txtatm[n].replace(22, 4, b);
              } else if (j == k + 1) {
                k = k + 1;
                char b[8];
                std::snprintf(b, sizeof b, "%4d", k + 1);
                txtatm[n].replace(22, 4, b);
              } else {
                break;
              }
              n = n + 1;
              if (n > natoms) break;
            }
            k = res;
            het = sl(allres[k], 1, 3);
            allres[k] = "UNK";
            while (true) {
              het2 = sl(allres[k + 1], 1, 3);
              if (het2 == "   ") break;
              allres[k + 1] = het;
              het = het2;
              k = k + 1;
            }
          }
        }
        //  Formula of hetero group is C(n_C) N(n_N) O(n_O)
        //  and atoms are in inres(:ninres)
        //
        //  Now assign "obvious" hetero groups
        het = "HET";
        line = "Not identified;";
        for (k = 3; k <= 84; ++k) {
          if (n_ele[k] > 0) {
            l = n_ele[k];
            line += cap_elemnt[k];
            if (l != 1) line += std::to_string(l);
          }
        }
        het_group = trim(line);
        //
        //  Make sure all oxygen atoms are counted
        kk = ninres;
        for (ii = 1; ii <= kk; ++ii) {
          k = inres[ii];
          if (nat[k] != 6) continue;
          jj = ibonds[nbonds[k] + 1][k];
          if (jj > 0) {
            if (nat[jj] == 8) {
              if (distance(k, jj) < 1.9) {
                for (k = 1; k <= ninres; ++k)
                  if (inres[k] == jj) break;
                if (k > ninres) {
                  n_O = n_O + 1;
                  ninres = ninres + 1;
                  inres[ninres] = jj;
                }
              }
            }
          }
        }
        switch (n_C) {
          case 0: {
            j = 0;
            l = 0;
            for (k = 1; k <= ninres; ++k)
              if (nat[inres[k]] != 1) { j = j + 1; l = inres[k]; }
            if (ninres == 1) {
              j = 0;
              switch (nat[i]) {
                case 3: case 4: case 11: case 12:
                case 19: case 20: case 21: case 22: case 23: case 24:
                case 25: case 26: case 27: case 28: case 29: case 30:
                case 37: case 38: case 39: case 40: case 41: case 42:
                case 43: case 44: case 45: case 46: case 47: case 48:
                case 55: case 56: case 57: case 58: case 59: case 60:
                case 61: case 62: case 63: case 64: case 65: case 66:
                case 67: case 68: case 69: case 70: case 71: case 72:
                case 73: case 74: case 75: case 76: case 77: case 78:
                case 79: case 80:
                  j = 1;
                  break;
                default:
                  break;
              }
              if (j == 1 || nbonds[i] == 0) {
                het = " " + cap_elemnt[nat[inres[1]]];
                if (het.size() >= 3 && het[2] >= 'a' && het[2] <= 'z')
                  het[2] = het[2] + ('A' - 'a');
                het_group = "Isolated element";
              } else if (nbonds[i] > 0) {
                het = " " + cap_elemnt[nat[inres[1]]];
                if (het.size() >= 3 && het[2] >= 'a' && het[2] <= 'z')
                  het[2] = het[2] + ('A' - 'a');
                het_group = "Covalently bound element";
              }
              k = ibonds[1][i];
              if (nbonds[i] == 1) {
                if (trim(txtatm[k]) != "") {
                  el = pdb_el(nat[i]);
                  if (nat[i] == 1) {
                    txtatm[i] = txtatm[k].substr(0, 12) + " H" +
                                txtatm[k].substr(14);
                  } else {
                    //  A single element is bonded to something.
                    if (trim(txtatm1[i]) == "") {
                      j = nint(reada(txtatm[k], 21)) + 1;
                      char b[16];
                      std::snprintf(b, sizeof b, "%8d", j);
                      line = b;
                    } else {
                      line = txtatm1[i].substr(14);
                    }
                    if (cap_elemnt[nat[i]].size() >= 2 &&
                        cap_elemnt[nat[i]][1] == ' ')
                      txtatm[i] = txtatm[k].substr(0, 12) + " " +
                                  cap_elemnt[nat[i]].substr(0, 1) + trim(line);
                    else
                      txtatm[i] = txtatm[k].substr(0, 12) + cap_elemnt[nat[i]] +
                                  txtatm[k].substr(14);
                  }
                  if ((int)txtatm[i].size() < 26)
                    txtatm[i] += std::string(26 - txtatm[i].size(), ' ');
                  continue;  // cycle i_loop
                }
              }
            }
            if (ninres == 2 && n_O == 1 && n_ele[1] == 1) {
              het = " OH";
              het_group = "Hydroxide ion";
            }
            if (ninres == 3 && n_O == 1 && n_ele[1] == 2) {
              het = "H2O";
              het_group = "Complexed water";
            }
            if (nat[l] == 7) {
              if (nheavy(l) > 1) {
                for (k = 1; k <= nbonds[l]; ++k) {
                  j = ibonds[k][l];
                  if (nat[j] != 6) continue;
                  for (m = 1; m <= nbonds[j]; ++m)
                    if (nat[ibonds[m][j]] == 8) break;
                  if (m <= nbonds[j]) {
                    txtatm[l] = txtatm[j].substr(0, 6) +
                                sl(txtatm[l], 7, 11) + "  N" +
                                txtatm[j].substr(14);
                    continue;  // cycle i_loop
                  }
                }
              } else {
                het = "NH3";
                het_group = "Ammonia";
              }
            }
            if (n_N == 0 && n_O == 2 && n_S == 0 && n_P == 0) {
              het = "O-O";
              het_group = "Hydrogen peroxide?";
            }
            if (n_N == 1 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "NO3";
              het_group = "Nitrate";
            }
            if (n_P == 1 && n_O == 4 && n_S == 0) {
              het = "PO4";
              het_group = "Phosphate";
            }
            break;
          }
          case 1:
            if (n_N == 0 && n_O == 1 && n_S == 0 && n_P == 0) {
              het = "MOH";
              het_group = "Methanol";
            }
            if (n_N == 0 && n_O == 2 && n_S == 0 && n_P == 0) {
              het = "FMT";
              het_group = "Formic acid";
            }
            if (n_N == 0 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "CO3";
              het_group = "Carbonate ion";
            }
            if (n_N == 1 && n_O == 0 && n_S == 0 && n_P == 0) {
              het = "CYN";
              het_group = "Cyanide ion";
            }
            if (n_N == 0 && n_O == 0 && n_ele[1] == 3 && n_S == 0 &&
                n_P == 0) {
              het = "CH3";
              het_group = "Methyl group";
            }
            break;
          case 2:
            if (n_N == 0 && n_O == 1 && n_S == 0 && n_P == 0) {
              het = "EOH";
              het_group = "Ethanol";
            } else if (n_N == 0 && n_O == 2 && n_S == 0 && n_P == 0) {
              for (ii = 1; ii <= ninres; ++ii)
                if (nat[inres[ii]] == 6) break;
              if (nheavy(inres[ii]) == 2) {
                het = "EDO";
                het_group = "1,2-Ethanediol ";
              } else {
                het = "ACY";
                het_group = "Acetic acid ";
              }
            }
            break;
          case 3:
            if (n_N == 0 && n_O == 2 && n_S == 0 && n_P == 0) {
              het = "PGO";
              het_group = "S-1,2-Propanediol";
            }
            if (n_N == 0 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "GOL";
              het_group = "Glycerol";
            }
            if (n_N == 0 && n_O == 4 && n_S == 0 && n_P == 0) {
              het = "MLI";
              het_group = "Malonate ion (-)";
            }
            if (n_N == 2 && n_O == 0 && n_S == 0 && n_P == 0) {
              het = "IMD";
              het_group = "Imidazole";
            }
            break;
          case 4:
            if (n_N == 0 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "PEG";
              het_group = "Di(hydroxyethyl) ether";
            }
            if (n_N == 0 && n_O == 6 && n_S == 0 && n_P == 0) {
              het = "TLA";
              het_group = "L(+)-Tartaric acid";
            }
            if (n_N == 1 && n_O == 0 && n_S == 0 && n_P == 0) {
              het = "NTB";
              het_group = "N-Tertiary butyl";
            }
            if (n_N == 1 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "TRS";
              het_group = "2-Amino-2-hydroxymethyl-propane-1,3-diol";
            }
            break;
          case 5:
            if (n_N == 1 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "PCA";
              het_group = "Pyroglutamic acid";
            }
            if (n_N == 1 && n_O == 7 && n_S == 0 && n_P == 0) {
              het = "10E";
              het_group = "4-Amino-3-methylbut-2-en-1-yl diphosphate";
            }
            break;
          case 6:
            if (n_N == 0 && n_O == 1 && n_P == 1) {
              het = "HXP";
              het_group = "Hexyl bonded to phosphorus";
            }
            if (n_N == 0 && n_O == 2 && n_S == 0 && n_P == 0) {
              het = "RCO";
              het_group = "Resorcinol";
            }
            if (n_N == 0 && n_O == 9 && n_P == 0) {
              het = "SGA";
              het_group = "O3-Sulfonylgalactose";
            }
            if (n_N == 0 && n_O == 5 && n_S == 0 && n_P == 0) {
              het = "FUC";
              het_group = "Alpha-L-fucose";
            }
            if (n_N == 0 && n_O == 6 && n_S == 0 && n_P == 0) {
              het = "HEX";
              het_group = "A hexose, e.g. glucose or mannose";
              identify_hexose(ninres, inres, het, het_group);
            }
            if (n_N == 0 && n_O == 2 && n_S == 0 && n_P == 0) {
              het = "MPD";
              het_group = "(4S)-2-Methyl-2,4-pentanediol";
            }
            if (n_N == 0 && n_O == 4 && n_S == 0 && n_P == 0) {
              het = "PGE";
              het_group = "Triethylene glycol";
            }
            if (n_N == 1 && n_O == 4 && n_P == 0) {
              het = "MES";
              het_group = "2-(N-Morpholino)-ethanesulfonic acid";
            }
            break;
          case 7:
            if (n_N == 1 && n_O == 1 && n_S == 0 && n_P == 0) {
              het = "DBZ";
              het_group = "Benzoylamino (H6C7NO)";
            }
            if (n_N == 1 && n_O == 3 && n_S == 0 && n_P == 0) {
              het = "TAM";
              het_group = "Tris(hydroxyethyl) aminomethane";
            }
            break;
          case 8:
            if (n_N == 0 && n_O == 0 && n_S == 0 && n_P == 0) {
              het = "OCT";
              het_group = "N-Octane";
            }
            if (n_N == 0 && n_O == 5 && n_S == 0 && n_P == 0) {
              het = "PG4";
              het_group = "Tetraethylene glycol";
            }
            if (n_N == 0 && n_O == 1 && n_S == 0 && n_P == 0) {
              het = "0VT";
              het_group = "6-Methyl-5-hepten-2-one";
            }
            if (n_N == 1 && (n_O > 3 && n_O < 7) && n_S == 0 && n_P == 0) {
              het = "NAG";
              het_group = "N-Acetyl-D-glucosamine";
            }
            if (n_N == 1 && n_O == 0 && n_S == 0 && n_P == 0) {
              het = "TBA";
              het_group = "Tetrabutylammonium ion";
            }
            if (n_N == 1 && n_O == 2 && n_P == 0) {
              het = "AES";
              het_group = "4-(2-Aminoethyl)benzenesulfonyl fluoride";
            }
            if (n_N == 2 && n_O == 4 && n_P == 0) {
              het = "EPE";
              het_group = "Hepes (C8 H18 N2 O4 S)";
            }
            break;
          case 9:
            if (n_N == 1 && n_O == 2 && n_ele[17] == 1 && n_S == 0 &&
                n_P == 0) {
              het = "GM5";
              het_group = "4-Chlorocinnamylhydroxamate";
            }
            if (n_N == 3 && n_O == 13 && n_S == 0 && n_P == 3) {
              het = "DCP";
              het_group = "2'-deoxycytidine-5'-triphosphate";
            }
            break;
          case 10:
            if (n_N == 0 && n_O == 6 && n_S == 0 && n_P == 0) {
              het = "TSA";
              het_group = "Endo-oxabicyclic transition state analogue";
            }
            if (n_N == 1 && n_O == 1 && n_S == 0 && n_P == 0) {
              het = "T55";
              het_group = "8-Methylnonanoic acid";
            }
            if (n_N == 1 && n_O == 3 && n_S == 0) {
              het = "VPF";
              het_group = "C22 H32 N3 O8 P (incomplete)";
            }
            if (n_N == 5 && n_O == 10 && n_S == 0 && n_P == 2) {
              het = "ADP";
              het_group = "Adenosine-5'-diphosphate";
            }
            if (n_N == 6 && n_O == 12 && n_S == 0) {
              het = "ANP";
              het_group = "Phosphoaminophosphonic acid-adenylate ester";
            }
            break;
          case 13:
            if (n_N == 2 && n_O == 0 && n_S == 0 && n_P == 0) {
              het = "THA";
              het_group = "Tacrine";
            }
            break;
          case 14:
            if (n_N == 0 && n_O == 4 && n_S == 0 && n_P == 0) {
              het = "HKA";
              het_group = "3-Methoxy-4-phenoxybenzoic acid";
            }
            break;
          case 15:
            if (n_N == 1 && n_O == 5 && n_S == 0) {
              het = "1CT";
              het_group = "Oxo-trifluoro methyl phenyl ethoxy phenyl phosphonic acid";
            }
            if (n_N == 2 && n_O == 12 && n_S == 0 && n_P == 2) {
              het = "5GW";
              het_group = "5-Phenyluridine 5'-(trihydrogen diphosphate)";
            }
            break;
          case 16:
            if (n_N == 2 && n_O == 11 && n_S == 0 && n_P == 0) {
              het = "NA2";
              het_group = "N-Acetyl-D-glucosamine, dimer";
            }
            if (n_N == 0 && n_O == 1 && n_S == 0 && n_P == 0) {
              het = "BOM";
              het_group = "Hexadecyl-10,12-dien-1-ol";
            }
            break;
          case 17:
            if (n_N == 4 && n_O == 9 && n_S == 0 && n_P == 1) {
              het = "FMN";
              het_group = "Riboflavin monophosphate ";
            }
            if (n_N == 2 && n_O == 2) {
              het = "641";
              het_group = "Oxopyrrolidine-3-carboxamide";
            }
            break;
          case 18:
            if (n_N == 2 && n_O == 0) {
              het = "HUX";
              het_group = "(-)-Huprine X (C18 H19 N2 Cl)";
            }
            break;
          case 19:
            if (n_N == 0 && n_O == 2) {
              het = "ASD";
              het_group = "4-Androstene-3-17-dione";
            }
            if (n_N == 4 && n_O == 2) {
              het = "PNT";
              het_group = "1,5-Bis(4-amidinophenoxy)pentane";
            }
            break;
          case 20:
            if (n_N == 0 && n_O == 1) {
              het = "ARC";
              het_group = "3,7,11,15-Tetramethyl-hexadecan-1-ol";
            }
            if (n_N == 0 && n_O == 10) {
              het = "BHE";
              het_group = "Galactopyranoside";
            }
            if (n_N == 4 && n_O == 2) {
              het = "DID";
              het_group = "4,4'[1,6-Hexanediylbis(oxy)]bisbenzenecarboximidamide";
            }
            if (n_N == 7 && n_O == 9) {
              het = "BT5";
              het_group = "Biotinyl-5-amp";
            }
            break;
          case 21:
            if (n_N == 2 && n_O == 17 && n_P == 2) {
              het = "2GW";
              het_group = "5-Phenyl-uridine-5'-alpha-D-galactosyl-diphosphate";
            }
            if (n_N == 7 && n_O == 14) {
              het = "NAD";
              het_group = "Nicotinamide-adenine-dinucleotide";
            }
            if (n_N == 7 && n_O == 17 && n_P > 0) {
              het = "NAP";
              het_group = "Nicotinamide-adenine-dinucleotide phosphate";
            }
            if (n_N == 10 && n_O == 1) {
              het = "32G";
              het_group = "C21 H30 N10 O S";
            }
            break;
          case 22:
            if (n_N == 2 && n_O == 2) {
              het = "CZM";
              het_group = "3,3'-Me2-salophen";
            }
            if (n_N == 3 && n_O == 8) {
              het = "VPF";
              het_group = "C22 H32 N3 O8 P";
            }
            break;
          case 23:
            if (n_N == 0 && n_O == 11) {
              het = "CM5";
              het_group = "5-Cyclohexyy-1-pentyl-beta-d-maltoside";
            }
            if (n_N == 3 && n_O == 3) {
              het = "4A2";
              het_group = "C23 H17 F4 N3 O3";
            }
            break;
          case 26:
            if (n_N == 0 && n_O == 8) {
              het = "0DV";
              het_group = "Fusicoccin H";
            }
            break;
          case 34:
            if (n_N == 4 && n_O == 4) {
              het = "HEM";
              het_group = "Heme ring";
            }
            break;
          default:
            if (n_C >= 35) {
              het = "BIG";
              het_group = "Large organic compound " + trim(het_group);
            }
            break;
        }
        if (n_ele[15] > 3) {
          het = "NUC";
          het_group = "Nucleic acid (DNA or RNA type)";
        }
        l = index1(keywrd, " XENO");
        if (l != 0) {
          int idx_xeno = l;
          l_chain = (index1(keywrd, " CHAINS=(") == 0);
          j = index1(keywrd.substr(idx_xeno - 1), ") ");
          if (j != 0) {
            line = " " + keywrd.substr(idx_xeno + 4, j - 5);
            while (true) {
              for (l = 1; l <= 10; ++l) {
                line = trim(line.substr(1));
                std::string c1 = line.substr(0, 1);
                if (c1 == "(" || c1 == "," || c1 == ";" || c1 == ")") break;
              }
              if (trim(line) == "") break;
              num = line.substr(1, 1);
              k = 0;
              if (l_chain) {
                if (num >= "A" && num <= "Z") {
                  if (sl(txtatm1[i], 22, 22) >= "A" &&
                      sl(txtatm1[i], 22, 22) <= "Z") {
                    if (num != sl(txtatm1[i], 22, 22)) k = -1000;
                  }
                }
              }
              k = nint(reada(line, 1)) + k;
              if (k == ires) {
                for (l = 1; l <= 10; ++l) {
                  if (line.substr(0, 1) == "=") break;
                  line = trim(line.substr(1));
                }
                if (len_trim(het_group) > 50)
                  het_group = "Defined using keyword XENO" + het_group.substr(61);
                else
                  het_group = "Defined using keyword XENO" + het_group.substr(14);
                if (first) {
                  first = false;
                  std::fprintf(stdout,
                               "\n      Ligand names that have been changed\n");
                  std::fprintf(stdout,
                               "      Residue No.  Calc'd name   XENO name\n");
                }
                change_no = change_no + 1;
                std::fprintf(stdout, "%3d%10d  %1s%7s%9s\n", change_no, ires,
                             num.c_str(), het.c_str(),
                             line.substr(1, 3).c_str());
                het = line.substr(1, 3);
                break;
              }
            }
          }
        }
        if (nat[i] != 1 || nbonds[i] != 1) {
          if (ncomments < max_comments && index1(keywrd, " RESID") != 0) {
            ncomments = ncomments + 1;
            std::string cm = "*REMARK   3   " + het + " = " +
                             het_group.substr(0, 47) + "res: " +
                             std::to_string(ires);
            all_comments[ncomments] = cm;
            if (chanel_C::log) {
              // ilog unit not mapped to a C++ stream yet (log==false in tests).
              if (het == "HET")
                std::fprintf(stdout, " (Includes atom number: %5d)\n",
                             inres[1]);
            }
          }
          for (j = 1; j <= ninres; ++j) {
            k = inres[j];
            el = pdb_el(nat[k]);
            if (attached)
              txtatm[k] = atom26("ATOM  ", k, el, het, res);
            else
              txtatm[k] = atom26("HETATM", k, el, het, ires);
          }
        }
      } else {
        for (j = 1; j <= ninres; ++j) {
          k = inres[j];
          el = pdb_el(nat[k]);
          txtatm[k] = atom26("HETATM", k, el, "UNK", ires);
        }
        inc_res(ires, start_res, nfrag);
      }
    }
  }
  //
  //  Label any remaining atoms
  //
  for (i = 1; i <= natoms - id; ++i) {
    if (sl(txtatm[i], 26, 26) == " ") {
      el = pdb_el(nat[i]);
      txtatm[i] = atom26nr("HETATM", i, el, "UNK");
    }
  }
}

// ---------------------------------------------------------------------------
void identify_hexose(int ninres, const std::vector<int>& inres,
                     std::string& nam, std::string& name) {
  static const char* hexose_aldose[8] = {
      "Allose", "Altrose", "Glucose",  "Mannose",
      "Gulose", "Idose",   "Galactose", "Talose"};
  static const char* hexose_ketose[4] = {"Psicose", "Fructose", "Sorbose",
                                         "Tagatose"};
  int i, j, k, l, m, C1 = 0, Cn, Cm, Cp, On, Hn;
  std::vector<int> backbone(7, 0), chiral(7, 0);
  bool aldose;
  double torsion = 0.0;

  // First, check that all carbon atoms have four ligands
  for (i = 1; i <= ninres; ++i) {
    if (nat[inres[i]] == 6) {
      if (nbonds[inres[i]] != 4) {
        if (nbonds[inres[i]] < 3) return;
        j = inres[i];
        k = ibonds[4][j];
        if (k == 0) return;
        if (nat[k] != 8) return;
      }
    }
  }
  aldose = true;
  //
  //  Locate C1
  //
  for (i = 1; i <= ninres; ++i) {
    if (nat[inres[i]] == 6) {
      C1 = inres[i];
      k = 0;
      for (j = 1; j <= 4; ++j)
        if (nat[ibonds[j][C1]] == 8) k = k + 1;
      if (k == 2) break;
    }
  }
  for (i = 1; i <= nbonds[C1]; ++i) {
    if (nat[ibonds[i][C1]] == 6) {
      j = ibonds[i][C1];
      if (nheavy(j) == 2) {
        l = 0;
        for (k = 1; k <= nbonds[j]; ++k)
          if (nat[ibonds[k][j]] == 8) l = l + 1;
        if (l == 1) {
          C1 = j;
          aldose = false;
          break;
        }
      }
    }
  }
  backbone[1] = C1;
  Cn = C1;
  Cm = Cn;
  //
  //  Now locate C2 - C6
  //
  for (i = 2; i <= 6; ++i) {
    for (j = 1; j <= 4; ++j) {
      l = ibonds[j][Cn];
      if (l == 0) return;
      if (nat[l] == 6 && l != Cm) break;
    }
    backbone[i] = l;
    Cm = backbone[i - 1];
    Cn = l;
  }
  chiral.assign(7, 0);
  //
  //  Work out chirality of C1
  //
  On = 0;
  for (i = 1; i <= 4; ++i) {
    j = ibonds[i][C1];
    if (j == 0) return;
    if (nat[j] == 8 && nheavy(j) == 2) {
      //  Make sure that the atom the oxygen is attached to is not in the ring
      for (l = 1; l <= nbonds[j]; ++l) {
        m = ibonds[l][j];
        if (m == 0) return;
        if (nat[m] > 1 && m != C1) break;
      }
      for (l = 1; l <= ninres; ++l)
        if (inres[l] == m) break;
      if (l > ninres) On = j;
    }
    if (nat[j] == 8 && nheavy(j) == 1) On = j;
    if (nat[j] == 1) Hn = j;
  }
  if (On == 0) return;  //  Failed to find oxygen
  Cp = backbone[2];
  dihed(coord, Hn, Cp, C1, On, torsion);
  if (torsion > pi) torsion = torsion - 2 * pi;
  if (torsion < 0) chiral[1] = 1;
  //
  //  Work out chirality of C2 - C5
  //
  for (k = 2; k <= 5; ++k) {
    Cn = backbone[k];
    Cp = backbone[k + 1];
    l = 0;
    for (i = 1; i <= 4; ++i)
      if (nat[ibonds[i][Cn]] == 8) l = l + 1;
    for (i = 1; i <= 4; ++i) {
      j = ibonds[i][Cn];
      if (l == 1) {
        if (nat[j] == 8) On = j;
        if (nat[j] == 1) Hn = j;
      } else {
        if (nat[j] == 8 && nheavy(j) == 2) Hn = j;
        if (nat[j] == 8 && nheavy(j) == 1) On = j;
      }
    }
    dihed(coord, Hn, Cn, Cp, On, torsion);
    if (torsion > pi) torsion = torsion - 2 * pi;
    if (torsion > 0) chiral[k] = 1;
  }
  //
  //  Determine the hexose
  //
  if (aldose) {
    i = chiral[2] + 2 * chiral[3] + 4 * chiral[4] + 1;
    name = hexose_aldose[i - 1];
  } else {
    i = chiral[3] + 2 * chiral[4] + 1;
    name = hexose_ketose[i - 1];
  }
  if (chiral[5] == 0)
    name = "D-" + name;
  else
    name = "L-" + name;
  if (aldose) {
    if (chiral[1] == 0)
      name = "alpha-" + name;
    else
      name = "beta-" + name;
  } else {
    if (chiral[2] == 0)
      name = "alpha-" + name;
    else
      name = "beta-" + name;
  }
  std::string key = name;
  nam = "HEX";
  if (key == "alpha-D-Allose") nam = "ALO";
  else if (key == "alpha-D-Altrose") nam = "ALT";
  else if (key == "alpha-D-Glucose") nam = "GLC";
  else if (key == "alpha-D-Mannose") nam = "MAN";
  else if (key == "alpha-D-Gulose") nam = "GUL";
  else if (key == "alpha-D-Idose") nam = "IDO";
  else if (key == "alpha-D-Galactose") nam = "GAL";
  else if (key == "alpha-D-Talose") nam = "TAL";
  else if (key == "beta-D-Allose") nam = "BAL";
  else if (key == "beta-D-Altrose") nam = "BAT";
  else if (key == "beta-D-Glucose") nam = "BGC";
  else if (key == "beta-D-Mannose") nam = "BMA";
  else if (key == "beta-D-Gulose") nam = "BGU";
  else if (key == "beta-D-Idose") nam = "BID";
  else if (key == "beta-D-Galactose") nam = "BGA";
  else if (key == "beta-D-Talose") nam = "BTA";
  else if (key == "alpha-D-Psicose") nam = "ADP";
  else if (key == "alpha-D-Fructose") nam = "ADF";
  else if (key == "alpha-D-Sorbose") nam = "ADS";
  else if (key == "alpha-D-Tagatose") nam = "ADT";
}

// ---------------------------------------------------------------------------
void moiety(std::vector<bool>& iopt, std::vector<int>& lused, int istart,
            int& n_new) {
  const int natomr = 800;
  int i2, i3, iatom, j, k, l, ninbit, nlive;
  std::vector<int> live(201, 0);
  std::vector<int> inres(natomr + 1, 0);

  iatom = istart;
  //
  //  The first atom identified in the moiety is atom IATOM.
  //
  iopt[iatom] = true;
  //
  //  NOW TO WORK OUT THE ATOMS IN THE MOIETY
  //
  nlive = nbonds[iatom];
  if (nlive == 0) {
    ninbit = 1;
    inres[1] = iatom;
  } else {
    ninbit = 0;
    for (i2 = 1; i2 <= nlive; ++i2) live[i2] = ibonds[i2][iatom];
    while (true) {
      l = live[1];
      if (iopt[l] || l == iatom) {
        if (nlive < 1) break;
        live[1] = live[nlive];
        nlive = nlive - 1;
      } else {
        iopt[l] = true;
        ninbit = ninbit + 1;
        if (ninbit > natomr) {
          std::fprintf(stdout, " There are more than%4d atoms in moiety \n",
                       natomr);
          std::fprintf(stdout, " Atoms in moiety\n");
          for (l = 1; l <= natomr; ++l)
            std::fprintf(stdout, " %s%5d", elemnt[nat[inres[l]]].c_str(),
                         inres[l]);
          mopend("Too many atoms in moiety");
          return;
        }
        inres[ninbit] = l;
        if (nbonds[l] != 0) {
          //  THERE IS AT LEAST ONE ATOM ATTACHED TO THE 'LIVE' ATOM
          for (i2 = 2; i2 <= nbonds[l]; ++i2) {
            j = ibonds[i2][l];
            bool dup = false;
            for (i3 = 1; i3 <= nlive; ++i3)
              if (live[i3] == j) { dup = true; break; }
            if (dup) continue;
            nlive = nlive + 1;
            live[nlive] = j;
          }
          live[1] = ibonds[1][l];
        } else {
          if (nlive == 0) break;
          live[1] = live[nlive];
          nlive = nlive - 1;
        }
      }
    }
  }
  //
  //  Check that atoms are not counted twice
  //
  for (j = 1; j <= ninbit; ++j)
    for (k = j + 1; k <= ninbit; ++k)
      if (inres[j] == inres[k]) inres[j] = 0;
  //
  //  Put hydrogen atoms at the end of the list
  //
  k = ninbit;
  for (j = 1; j <= ninbit; ++j) {
    l = inres[j];
    if (l != 0) {
      if (nat[l] == 1) {
        k = k + 1;
        inres[k] = l;
        inres[j] = 0;
      }
    }
  }
  ninbit = k;
  if (iatom != 0) {
    n_new = n_new + 1;
    lused[n_new] = iatom;
    iopt[iatom] = true;
  }
  for (i2 = 1; i2 <= ninbit; ++i2) {
    j = inres[i2];
    if (j != 0 && j != iatom) {
      n_new = n_new + 1;
      lused[n_new] = j;
    }
  }
}
