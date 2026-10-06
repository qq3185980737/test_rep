// geochk_txtatm.cpp — C++ translation of geochk.F90 subprograms:
//   update_txtatm   (lines 1821-1983),
//   rectify_sequence (lines 1984-2067),
//   compare_sequence (lines 1725-1820),
//   write_sequence   (lines 2068-2273).
// 1-based Fortran indexing is kept for all module arrays.
#include "geochk.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "lewis.h"
#include "mopend.h"
#include "molkst_C.h"
#include "reada.h"
#include "set_up_dentate.h"

namespace mc = common_arrays_C;
namespace mk = molkst_C;
namespace mz = MOZYME_C;

static int iw_out = 6;  // channel_C::iw (stdout in tests)

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

// Fortran "c == ' '" for a fixed-length character variable: all blanks.
static bool is_blank(const std::string& s) {
  for (char c : s) if (c != ' ') return false;
  return true;
}

// Fortran write s(a:b) = v  (1-based inclusive), padding short v with blanks.
static void setsub(std::string& s, int a, int b, const std::string& v) {
  if (a < 1) a = 1;
  if (b > (int)s.size()) { s.resize(b, ' '); }
  int n = b - a + 1;
  std::string t = v;
  if ((int)t.size() > n) t.resize(n);
  while ((int)t.size() < n) t += ' ';
  s.replace(a - 1, n, t);
}

static int nint(double x) { return static_cast<int>(std::floor(x + 0.5)); }

// ===========================================================================
//  subroutine update_txtatm(output, sort)
// ===========================================================================
void update_txtatm(bool output, bool sort) {
  if (mk::maxtxt != 26) return;

  bool pdb, update_chain, L_all;
  int i, j, k, l, H_Z;

  if (mk::keywrd.find(" RESID") != std::string::npos ||
      mk::keywrd.find(" ADD-H") != std::string::npos ||
      mk::keywrd.find(" SITE=") != std::string::npos ||
      mk::keywrd.find(" RESEQ") != std::string::npos || output) {
    H_Z = 1;
  } else {
    H_Z = 0;
  }
  pdb = (mk::keywrd.find(" PDBOUT") != std::string::npos);
  L_all = (output && mk::keywrd.find(" RESID") == std::string::npos);

  if (L_all && !sort && mk::numat == mk::numat_old) {
    // Do nothing!
    for (i = 1; i <= mk::numat; ++i) {
      if (is_blank(mc::txtatm1[i])) mc::txtatm1[i] = mc::txtatm[i];
      else                          mc::txtatm[i] = mc::txtatm1[i];
    }
    goto label_99;
  }

  // Add text to TXTATM to label hydrogen atoms and to add chain letter.
  set_up_dentate();
  check_CVS(false);
  check_h(i);
  update_chain = (mk::keywrd.find(" CHAINS=(") == std::string::npos);
  for (i = 1; i <= mk::numat; ++i) {
    if (sort) {
      if (mc::nat[i] > H_Z) {
        if (L_all || update_chain) {
          for (j = 1; j <= mk::numat_old; ++j) {
            if (std::abs(mc::coord[0][i] - mc::coorda[0][j]) > 0.01) continue;
            if (std::abs(mc::coord[1][i] - mc::coorda[1][j]) > 0.01) continue;
            if (std::abs(mc::coord[2][i] - mc::coorda[2][j]) > 0.01) continue;
            if (L_all) {
              if (!is_blank(mc::txtatm1[j])) mc::txtatm[i] = mc::txtatm1[j];
            } else {
              if (!is_blank(mc::txtatm1[j])) setsub(mc::txtatm[i], 22, 22,
                                                    fsub(mc::txtatm1[j], 22, 22));
            }
            break;
          }
        }
      } else {
        if (mc::nbonds[i] > 0) mc::txtatm[i] = mc::txtatm[mc::ibonds[1][i]];
      }
    } else {
      if (L_all) {
        if (mc::nat[i] > 1) {
          if (!is_blank(mc::txtatm1[i])) mc::txtatm[i] = mc::txtatm1[i];
        } else if (mc::nbonds[i] > 0) {
          mc::txtatm[i] = mc::txtatm[mc::ibonds[1][i]];
        }
      } else {
        if (update_chain && !is_blank(mc::txtatm1[i]))
          setsub(mc::txtatm[i], 22, 22, fsub(mc::txtatm1[i], 22, 22));
      }
    }
  }

label_99:
  // Number hydrogen atoms, if more than one on a heavy atom.
  std::vector<int> n_H(mk::numat + 1, 0), nn_H(mk::numat + 1, 0);
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] > H_Z) {
      j = 0;
      for (k = 1; k <= mc::nbonds[i]; ++k) {
        if (mc::nat[mc::ibonds[k][i]] == 1) j = j + 1;
      }
      n_H[i] = j;
    }
  }
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] <= H_Z && mc::nbonds[i] > 0) {
      k = mc::ibonds[1][i];
      // txtatm(i) = txtatm(k)(:12)//" H"//txtatm(k)(15:)
      mc::txtatm[i] = fsub(mc::txtatm[k], 1, 12) + " H" + fsub(mc::txtatm[k], 15, 26);
      // If a hydrogen atom is attached to an unlabeled carbon atom, make the
      // hydrogen atom a terminal hydrogen.
      if (mc::nat[k] == 6 && fsub(mc::txtatm[i], 15, 15) == " ")
        setsub(mc::txtatm[i], 15, 15, "T");
      if (n_H[k] == 1) {
        setsub(mc::txtatm[i], 13, 14, " H");
      } else {
        nn_H[k] = nn_H[k] + 1;
        setsub(mc::txtatm[i], 13, 14, std::string(1, (char)(nn_H[k] + '0')) + "H");
      }
    }
    if (mk::numat == mk::numat_old) {
      if (is_blank(mc::txtatm1[i])) mc::txtatm1[i] = mc::txtatm[i];
    }
  }

  // Check for duplicate hydrogen atom labels.
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] == 1) {
      l = 1;
      for (j = i + 1; j <= mk::numat; ++j) {
        if (mc::nat[j] == 1) {
          if (fsub(mc::txtatm[i], 13, 26) == fsub(mc::txtatm[j], 13, 26)) {
            for (k = 13; k <= 16; ++k) {
              if (fsub(mc::txtatm[j], k, k) == " ") {
                l = l + 1;
                setsub(mc::txtatm[j], k, k, std::string(1, (char)(l + '0')));
                if (mk::numat == mk::numat_old) {
                  if (is_blank(mc::txtatm1[j])) mc::txtatm1[j] = mc::txtatm[j];
                }
                break;
              }
            }
          }
        }
        if (l == 9) break;
      }
      if (l > 1) {
        for (k = 13; k <= 16; ++k) {
          if (fsub(mc::txtatm[i], k, k) == " ") {
            setsub(mc::txtatm[i], k, k, "1");
            if (mk::numat == mk::numat_old) {
              if (is_blank(mc::txtatm1[i])) mc::txtatm1[i] = mc::txtatm[i];
            }
            break;
          }
        }
      }
    }
  }

  // Add atom numbering, using PDB format.
  j = 1;
  for (i = 1; i <= mk::numat; ++i) {
    // write(txtatm(i),'(a6,i5,a15)')txtatm(i)(:6), i+j-1, txtatm(i)(12:)
    char buf[32];
    std::snprintf(buf, sizeof buf, "%6s%5d%15s",
                  fsub(mc::txtatm[i], 1, 6).c_str(), i + j - 1,
                  fsub(mc::txtatm[i], 12, 26).c_str());
    mc::txtatm[i] = buf;
    if (pdb && i == mc::breaks[j]) j = j + 1;
  }
}

// ===========================================================================
//  subroutine rectify_sequence
// ===========================================================================
void rectify_sequence() {
  if (mk::maxtxt != 26) return;
  int i, i_lower, res2, res1, res3, res, ii, jj, num_min, num_max;
  char chain1, chain2, chain;
  std::vector<char> used(mk::numat + 1, ' ');

  mc::txtatm1.resize(mk::numat + 1);
  for (i = 1; i <= mk::numat; ++i) mc::txtatm1[i] = mc::txtatm[i];
  for (i = 1; i <= mk::numat; ++i) used[i] = fsub(mc::txtatm1[i], 22, 22)[0];
  used[1] = 'a';
  mc::labels.resize(mk::numat + 1);
  for (i = 1; i <= mk::numat; ++i) mc::labels[i] = mc::nat[i];
  for (i = 1; i <= mk::numat; ++i)
    for (int k = 1; k <= 3; ++k) mc::geo[k][i] = mc::coord[k-1][i];

  res1 = nint(reada(mc::txtatm[1], 23));
  chain1 = 'a';
  ii = 1;
  for (i_lower = 2; i_lower <= mk::numat; ++i_lower) {
    chain2 = used[i_lower];
    if (chain2 == 'a') continue;  // atom already used
    res2 = nint(reada(mc::txtatm1[i_lower], 23));

    // Look for a fault at the end of the first bit of chain.
    for (res3 = res1; res3 <= res1 + 1; ++res3) {
      num_min = (i_lower - 50 > 1) ? (i_lower - 50) : 1;
      num_max = (i_lower + 50 < mk::numat) ? (i_lower + 50) : mk::numat;
      while (true) {
        for (jj = num_min; jj <= num_max; ++jj) {
          res = nint(reada(mc::txtatm1[jj], 23));
          chain = used[jj];
          if (res == res3 && chain == chain1 && used[jj] != 'a') break;
        }
        if (jj > num_max) break;
        ii = ii + 1;
        for (int k = 1; k <= 3; ++k) mc::coord[k-1][ii] = mc::geo[k][jj];
        mc::txtatm[ii] = mc::txtatm1[jj];
        mc::nat[ii] = mc::labels[jj];
        used[jj] = 'a';
      }
    }
    // Look for a fault at the start of the next bit of chain.
    for (res3 = res2 - 1; res3 <= res2; ++res3) {
      num_min = (i_lower - 50 > 1) ? (i_lower - 50) : 1;
      num_max = (i_lower + 50 < mk::numat) ? (i_lower + 50) : mk::numat;
      while (true) {
        for (jj = num_min; jj <= num_max; ++jj) {
          res = nint(reada(mc::txtatm1[jj], 23));
          chain = used[jj];
          if (res == res3 && chain == chain1 && used[jj] != 'a') break;
        }
        if (jj > num_max) break;
        ii = ii + 1;
        for (int k = 1; k <= 3; ++k) mc::coord[k-1][ii] = mc::geo[k][jj];
        mc::txtatm[ii] = mc::txtatm1[jj];
        mc::nat[ii] = mc::labels[jj];
        used[jj] = 'a';
      }
    }
    chain1 = chain2;
    res1 = res2;
  }
}

// ===========================================================================
//  subroutine compare_sequence(n_new)
// ===========================================================================
void compare_sequence(int n_new) {
  int i_atom = 0, i, j, i_delta = 0, new_res, old_res, previous = -200;
  int mbreaks, loop;
  std::string old, new_;
  bool first = true;
  std::string num;

  mbreaks = 1;
  i_atom = 0;
  for (loop = 1; loop <= mk::numat; ++loop) {
    // Find first non-hydrogen atom.
    while (true) {
      i_atom = i_atom + 1;
      if (i_atom > mk::numat) break;
      if (i_atom == mc::breaks[mbreaks]) mbreaks = mbreaks + 1;
      if (mc::nat[i_atom] != 1) break;
    }
    if (i_atom > mk::numat) break;
    new_res = nint(reada(mc::txtatm[i_atom], 23));
    new_ = fsub(mc::txtatm[i_atom], 18, 20);
    if (mk::maxtxt == 14) {
      old_res = nint(reada(mc::txtatm1[i_atom], 12));
      old = fsub(mc::txtatm1[i_atom], 8, 10);
    } else {
      old_res = nint(reada(mc::txtatm1[i_atom], 23));
      old = fsub(mc::txtatm1[i_atom], 18, 20);
    }
    if (new_res - i_delta != old_res) i_delta = new_res - old_res;
    if (new_res != previous && n_new == 0) {
      if (old != new_) {
        if (first) {
          std::fprintf(stdout, "\n%16s\n%7s\n", "Residue names that have changed",
                       "Original residue name   Calculated residue name");
          first = false;
          mk::line = "XENO=(";
        }
        for (i = 1; i <= 20; ++i) {
          if (old == mz::tyres[i]) break;
        }
        num = (i == 21) ? " " : mz::tyr[i];
        std::fprintf(stdout, "%14d%s  %-3s      %11d%s  %-3s\n", old_res,
                     (" " + std::string(1, mc::chains[mbreaks - 1])).c_str(), old.c_str(),
                     new_res, (" " + std::string(1, mc::chains[mbreaks - 1])).c_str(),
                     new_.c_str());
        if (i < 21) {
          j = static_cast<int>(len_trim(mk::line)) + 1;
          if (j > 200) break;
          char cbuf[32];
          if (old_res > -1) {
            int w = 1 + (int)std::floor(std::log10(old_res + 0.05));
            std::snprintf(cbuf, sizeof cbuf, "%c%*d=%c,", mc::chains[mbreaks - 1], w,
                          old_res, mz::tyr[i][0]);
          } else {
            int w = 1 + (int)std::floor(std::log10(-old_res + 0.05)) + 1;
            std::snprintf(cbuf, sizeof cbuf, "%c%*d=%c,", mc::chains[mbreaks - 1], w,
                          old_res, mz::tyr[i][0]);
          }
          mk::line.append(cbuf);
        }
      }
      previous = new_res;
    }
  }
  if (first && n_new == 0) {
    std::fprintf(stdout, "\n%s\n",
                 "        Calculated and original residue sequences agree perfectly");
  } else {
    j = static_cast<int>(len_trim(mk::line));
    if (j > 6 && n_new == 0) {
      setsub(mk::line, j, j, ")");
      mk::line = "(Use the XENO keyword to re-define unrecognized residues.)";
      std::fprintf(stdout, "\n%2s%s\n", "", mk::line.c_str());
    }
  }
}

// ===========================================================================
//  subroutine write_sequence
// ===========================================================================
void write_sequence() {
  static bool prt = true;
  static int icalcn = -50;
  int i, nfrag, jj, ii, ires, kl, ku, irold, l, j, k, charge, ifrag;
  std::string chain, ch;

  if (!prt) return;
  if (icalcn == mk::numcal) {
    return;
  } else {
    icalcn = mk::numcal;
    prt = false;
  }

  // First pass: write out residue names (4-char: 3 letters + charge).
  ii = 1;
  ifrag = 0;
  for (nfrag = 1; nfrag <= 100; ++nfrag) {
    if (ii > mk::numat) break;
    if (is_blank(mc::txtatm[ii])) break;
    chain = std::string(1, fsub(mc::txtatm[ii], 22, 22)[0]);
    ires = nint(reada(mc::txtatm[ii], 23));
    charge = mz::ions[ii];
    irold = ires;
    j = 0;
    for (ii = ii + 1; ii <= mk::numat; ++ii) {
      if (mc::nat[ii] != 1 && fsub(mc::txtatm[ii], 22, 22)[0] != chain[0]) break;
      jj = nint(reada(mc::txtatm[ii], 23));
      if (mc::nat[ii - 1] != 1) j = ii - 1;
      charge = charge + mz::ions[ii];
      if (ires == jj) continue;                 // residue is same
      charge = charge - mz::ions[ii];
      if (fsub(mc::txtatm[ii], 14, 14) == "H") continue;  // ignore hydrogen
      if (ires + 1 != jj) break;                // residue is not contiguous
      // Take residue name from the previous atom.
      mz::allres[ires] = fsub(mc::txtatm[ii - 1], 18, 20);
      if (charge == 1) setsub(mz::allres[ires], 4, 4, "+");
      else if (charge == -1) setsub(mz::allres[ires], 4, 4, "-");
      charge = 0;
      ires = ires + 1;
      if (ires > mz::maxres) {
        char b[64];
        std::snprintf(b, sizeof b, " Maximum residue number allowed: %d", mz::maxres);
        mopend(std::string(b));
        return;
      }
    }
    if (j != 0) {
      if (fsub(mc::txtatm[j], 18, 20) != "   ") mz::allres[ires] = fsub(mc::txtatm[j], 18, 20);
    }
    if (charge == 1) setsub(mz::allres[ires], 4, 4, "+");
    else if (charge == -1) setsub(mz::allres[ires], 4, 4, "-");

    for (j = irold; j <= ires; ++j) {
      for (k = 1; k <= 20; ++k) {
        if (fsub(mz::allres[j], 1, 3) == mz::tyres[k]) break;
      }
      if (k < 21) break;
    }
    if (j > ires || (j == ires && j == 0)) continue;
    ifrag = ifrag + 1;
    if (ifrag == 1) {
      std::fprintf(stdout, "\n%16s\n", ("RESIDUE SEQUENCE IN PROTEIN Chain: " + chain).c_str());
    } else {
      std::fprintf(stdout, "\n%16s%2d%s\n", "RESIDUE SEQUENCE IN PROTEIN FRAGMENT:", ifrag,
                   (" Chain: " + chain).c_str());
    }
    if (irold < 0) {
      jj = 100;
      i = (irold + jj) % 10;
      if (i == 0) i = 10;
      kl = irold;
      ku = (ires < kl - i + 10) ? ires : (kl - i + 10);
      mk::line = " ";
      l = 6 * (i - 1) + 1;
      j = ((kl + jj) / 10) * 10 - jj;
      if (i == 10) j = j - 10;
      std::fprintf(stdout, "%8s", "");
      for (k = 1; k <= 10; ++k) std::fprintf(stdout, "%6d", k - 10);
      std::fprintf(stdout, "\n%5d%2s", j, fsub(mk::line, 1, l).c_str());
      for (k = kl; k <= ku; ++k) std::fprintf(stdout, "%-4s  ", mz::allres[k].c_str());
      std::fprintf(stdout, "\n");
      while (true) {
        kl = ku + 1;
        if (kl > (ires < 0 ? ires : 0)) break;
        j = j + 10;
        ku = (ires < ku + 10) ? ires : (ku + 10);
        std::fprintf(stdout, "%5d%3s", j, "");
        for (k = kl; k <= ku; ++k) std::fprintf(stdout, "%-4s  ", mz::allres[k].c_str());
        std::fprintf(stdout, "\n");
      }
      std::fprintf(stdout, "\n");
    }
    kl = (irold > 1) ? irold : 1;
    i = kl % 10;
    if (i == 0) i = 10;
    ku = (ires < kl - i + 10) ? ires : (kl - i + 10);
    mk::line = " ";
    l = 6 * (i - 1) + 1;
    j = (kl / 10) * 10;
    if (i == 10) j = j - 10;
    std::fprintf(stdout, "%8s", "");
    for (k = 1; k <= 10; ++k) std::fprintf(stdout, "%6d", k);
    std::fprintf(stdout, "\n%5d%2s", j, fsub(mk::line, 1, l).c_str());
    for (k = kl; k <= ku; ++k) std::fprintf(stdout, "%-4s  ", mz::allres[k].c_str());
    std::fprintf(stdout, "\n");
    while (true) {
      kl = ku + 1;
      if (kl > ires) break;
      j = j + 10;
      ku = (ires < ku + 10) ? ires : (ku + 10);
      std::fprintf(stdout, "%5d%3s", j, "");
      for (k = kl; k <= ku; ++k) std::fprintf(stdout, "%-4s  ", mz::allres[k].c_str());
      std::fprintf(stdout, "\n");
    }
  }

  // Second pass: one-letter residue codes.
  ii = 1;
  ifrag = 0;
  for (nfrag = 1; nfrag <= 100; ++nfrag) {
    j = 0;
    if (ii > mk::numat) break;
    if (is_blank(mc::txtatm[ii])) break;
    chain = std::string(1, fsub(mc::txtatm[ii], 22, 22)[0]);
    ires = nint(reada(mc::txtatm[ii], 23));
    irold = ires;
    for (ii = ii + 1; ii <= mk::numat; ++ii) {
      if (mc::nat[ii] != 1 && fsub(mc::txtatm[ii], 22, 22)[0] != chain[0]) break;
      jj = nint(reada(mc::txtatm[ii], 23));
      if (mc::nat[ii - 1] != 1) j = ii - 1;
      if (ires == jj) continue;
      if (fsub(mc::txtatm[ii], 14, 14) == "H") continue;
      if (ires + 1 != jj) break;
      mz::allres[ires] = fsub(mc::txtatm[ii - 1], 18, 20);
      ires = ires + 1;
    }
    if (j != 0) {
      if (fsub(mc::txtatm[j], 18, 20) != "   ") mz::allres[ires] = fsub(mc::txtatm[j], 18, 20);
    }
    for (j = irold; j <= ires; ++j) {
      for (k = 1; k <= 20; ++k) {
        if (fsub(mz::allres[j], 1, 3) == mz::tyres[k]) break;
      }
      if (k < 21) break;
    }
    if (j > ires || (j == ires && j == 0)) continue;
    for (i = irold; i <= ires; ++i) {
      for (k = 1; k <= 20; ++k) {
        if (mz::tyres[k] == fsub(mz::allres[i], 1, 3)) break;
      }
      if (k > 0 && k < 21) mz::allr[i] = mz::tyr[k];
      else if (!is_blank(mz::allres[i])) mz::allr[i] = "X";
      else mz::allr[i] = " ";
    }
    ifrag = ifrag + 1;
    if (ifrag == 1) {
      std::fprintf(stdout, "\n%16s\n", ("RESIDUE SEQUENCE IN PROTEIN Chain: " + chain).c_str());
    } else {
      std::fprintf(stdout, "\n%16s%2d%s\n", "RESIDUE SEQUENCE IN PROTEIN FRAGMENT:", ifrag,
                   (" Chain: " + chain).c_str());
    }
    jj = 100;
    i = (irold + jj) % 10;
    if (irold < 0) i = i - 10;
    if (i == 0) i = 10;
    kl = irold;
    if (i < 0) ku = (ires < kl - i + 40) ? ires : (kl - i + 40);
    else       ku = (ires < kl - i + 50) ? ires : (kl - i + 50);
    mk::line = " ";
    j = ((kl + jj) / 10) * 10 - jj;
    if (i == 10) j = j - 10;
    if (i == 1) {
      ch = "10";
    } else {
      ch = " ";
      if (i < 0) {
        ch += (char)(9 + i + '0');
        if (ch[1] == '0') { ch[0] = '1'; i = 12; }
        else i = 2 - i;
      } else {
        ch += (char)(11 - i + '0');
      }
    }
    if (i != 0) {
      std::fprintf(stdout, "%5d%s", j, fsub(mk::line, 1, i + 1).c_str());
      for (k = kl; k <= ku; ++k) {
        if ((k - kl) % 10 == 0) std::fprintf(stdout, " ");
        std::fprintf(stdout, "%c", mz::allr[k][0]);
      }
      std::fprintf(stdout, "\n");
    }
    while (true) {
      kl = ku + 1;
      if (kl > ires) break;
      j = j + 50;
      ku = (ires < kl + 49) ? ires : (kl + 49);
      std::fprintf(stdout, "%5d%2s", j, "");
      for (k = kl; k <= ku; ++k) {
        if ((k - kl) % 10 == 0) std::fprintf(stdout, " ");
        std::fprintf(stdout, "%c", mz::allr[k][0]);
      }
      std::fprintf(stdout, "\n");
    }
  }
}
