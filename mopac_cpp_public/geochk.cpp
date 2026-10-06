// geochk.cpp — C++ translation of geochk.F90 "subroutine geochk" (lines
// 1-1489).  Main driver: Lewis-structure check, ionization identification,
// system-charge calculation and optional residue re-sequencing.
// 1-based Fortran indexing is kept for all module arrays.
#include "geochk.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "add_hydrogen_atoms.h"
#include "chanel_C.h"
#include "chkion.h"
#include "chklew.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "extvdw_for_MOZYME.h"
#include "find_salt_bridges.h"
#include "findn1.h"
#include "geout.h"
#include "lewis.h"
#include "ligand.h"
#include "lyse.h"
#include "mod_atomradii.h"
#include "molkst_C.h"
#include "mopend.h"
#include "names.h"
#include "newflg.h"
#include "parameters_C.h"
#include "pdbout.h"
#include "reada.h"
#include "reseq.h"
#include "timer.h"
#include "upcase.h"
#include "web_message.h"
#include "xyzint.h"

namespace mc = common_arrays_C;
namespace mk = molkst_C;
namespace mz = MOZYME_C;
namespace ch = chanel_C;
namespace pc = parameters_C;
namespace ec = elemts_C;
namespace mr = mod_atomradii;

// ---- file-local Fortran semantics helpers ---------------------------------
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

static bool is_blank(const std::string& s) {
  for (char c : s) if (c != ' ') return false;
  return true;
}

static int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

static int nint(double x) { return static_cast<int>(std::floor(x + 0.5)); }

// ion_names(-6..6) from geochk.F90 data statement.
static const char* ion_names[13] = {
    "Less than -5", "Penta-anion ", "Tetra-anion ", "Tri-anion   ",
    "Di-anion    ", "Anion       ", "(not used)  ", "Cation      ",
    "Di-cation   ", "Tri-cation  ", "Tetra-cation", "Penta-cation",
    "More than +5"};

void geochk() {
  int i, j, k, l, ii, jj, m, n, io, kk, ibad, ichrge, irefq, ires, nfrag;
  int uni_res, mres, maxtxt_store, nn1, n_new, j2, mbreaks, alloc_stat, jbad = 0;
  int n1 = 0, new_ = 0;
  int max_sites = 400;
  std::string padding(40, ' '), txtatm_1(26, ' '), txtatm_2(26, ' '), tmp(130, ' ');
  std::string num(1, ' ');
  std::vector<std::vector<std::string>> Lewis_formatted;
  std::vector<char> atom_charge;
  std::vector<double> radius;
  static std::vector<int> numbon(4, 0);   // 1-based [1..3]
  static int num_ions[13] = {0};          // -6..6 -> [0..12]
  std::vector<int> near_ions(101, 0);
  std::vector<double> r_ions(101, 0.0);
  std::vector<std::string> new_name(max_sites + 1, "   ");
  std::vector<std::string> old_name(max_sites + 1, "   ");
  std::vector<char> new_chain(max_sites + 1, ' ');
  std::vector<int> new_res(max_sites + 1, 0);
  std::vector<std::vector<char>> charge(max_sites + 1, std::vector<char>(4, '*'));
  std::vector<int> iopt(mk::maxatoms + 1, 0);   // geochk.F90: iopt(maxatoms) (local here)
  bool debug = false, let = false, lres = false, lreseq = false, times = false;
  bool opend = false, charges = false, l_protein = false, done = false;
  bool lsite = false, ter = false, residues = false, lbreaks = false, first = false;
  std::vector<bool> neutral(101, false);
  std::vector<bool> l_names(max_sites + 1, false);
  std::vector<int> mb;
  std::vector<int> nnbonds;
  std::vector<std::vector<int>> iibonds;
  double sum = 0.0;
  std::vector<double> r_ions_save;

  // ---- allocations --------------------------------------------------------
  mz::ions.assign(mk::maxatoms + 1, 0);
  mz::iz.assign(mk::maxatoms + 1, 0);
  mz::ib.assign(mk::maxatoms + 1, 0);
  mb.assign(mk::maxatoms + 1, 0);
  mz::at_res.assign(mk::maxatoms + 1, 0);
  atom_charge.assign(mk::maxatoms + 1, ' ');
  radius.assign(mk::maxatoms + 1, 0.0);
  bool* ioptl = new bool[mk::maxatoms + 1];
  for (i = 0; i <= mk::maxatoms; ++i) ioptl[i] = false;

  l_protein = false;
  for (i = 1; i <= mk::maxatoms; ++i) mz::ib[i] = 0;
  for (i = 1; i <= mk::maxatoms; ++i) mz::at_res[i] = 0;
  mres = 0;
  for (i = 1; i <= mk::maxatoms; ++i) mz::ions[i] = 0;
  j = 0;
  lbreaks = (mc::breaks.size() > 1 && mc::breaks[1] != -300);
  if (!lbreaks) mbreaks = 0;
  for (i = 1; i <= mk::ncomments; ++i) {
    if (mc::all_comments[i].find("REMARK   2") == std::string::npos) {
      j = j + 1;
      mc::all_comments[j] = mc::all_comments[i];
    }
  }
  mk::ncomments = j;

  done = false;
  n_new = 0;
  mk::line = mk::keywrd;

  // ---- XENO keyword parsing ----------------------------------------------
  i = index1(mk::line, " XENO");
  if (i != 0) {
    while (true) {
      if (fsub(mk::line, i, i) == "(") break;
      mk::line = mk::line.substr(1);
    }
    j = index1(mk::line.substr(i - 1), ") ");
    if (j == 0) {
      mopend("Closing parenthesis for XENO not found");
      delete[] ioptl;
      return;
    }
    mk::line = "," + mk::line.substr(i, j - 1);
    while (true) {
      if (is_blank(mk::line)) break;
      for (j = 1; j <= 100; ++j) {
        if (fsub(mk::line, 1, 1) == "," || fsub(mk::line, 1, 1) == ";" ||
            fsub(mk::line, 1, 1) == ")") break;
        mk::line = mk::line.substr(1);
      }
      if (fsub(mk::line, 1, 1) == ")") break;
      mk::line = mk::line.substr(1);
      k = index1(mk::line, "=");
      if (k == 0) {
        mopend("Equals sign (\"=\") expected in XENO but not found");
        delete[] ioptl;
        return;
      }
      l = k + 2;
      j = k + 4;
      n_new = n_new + 1;
      new_chain[n_new] = fsub(mk::line, 1, 1)[0];
      if (new_chain[n_new] >= '0' && new_chain[n_new] <= '9') new_chain[n_new] = 'A';
      new_res[n_new] = nint(reada(mk::line, 1));
      if (fsub(mk::line, l, l) == "," || fsub(mk::line, l, l) == ";" ||
          fsub(mk::line, l, l) == ")") {
        // A one-letter residue name has been detected.
        new_name[n_new] = fsub(mk::line, k + 1, k + 1);
      } else if (fsub(mk::line, j, j) == "," || fsub(mk::line, j, j) == ";" ||
                 fsub(mk::line, j, j) == ")") {
        // A three-letter residue name has been detected.
        new_name[n_new] = fsub(mk::line, k + 1, k + 3);
      } else {
        break;
      }
    }
  }
  while (true) {
    i = index1(mk::keywrd, " xeno");
    if (i == 0) break;
    mk::keywrd.replace(i - 1, 7, " XENO=(");
  }
  for (i = 1; i <= n_new; ++i) {
    for (j = i + 1; j <= n_new; ++j) {
      if (new_chain[i] == new_chain[j] && new_res[i] == new_res[j] &&
          new_name[i] != new_name[j]) {
        char b[128];
        std::snprintf(b, sizeof b, " XENO residue %c%4d occurs more than once.  Remove extra definition.",
                      new_chain[i], new_res[i]);
        mopend(std::string(b));
        delete[] ioptl;
        return;
      }
    }
    if (fsub(new_name[i], 2, 3) == "  ") {
      for (j = 1; j <= 20; ++j) {
        if (mz::tyr[j] == fsub(new_name[i], 1, 1)) {
          new_name[i] = mz::tyres[j];
          break;
        }
      }
    }
  }

  // ---- Add or remove hydrogen atoms, as necessary -------------------------
  i = index1(mk::keywrd, " SITE=(IONIZE)");
  if (i > 0) {
    mk::line = "SITE=(COO,NH3,ARG(+),SO4,PO4)";
    mk::keywrd = mk::keywrd.substr(0, i) + mk::line +
                 mk::keywrd.substr(i + 14);
  }
  i = index1(mk::keywrd, " SITE=");
  lsite = false;

  i = index1(mk::keywrd, " SITE=");
  if (i != 0 && mk::keywrd.find(" ADD-H") == std::string::npos) {
    lewis(true);
    while (true) {
      i = index1(mk::keywrd, " SITE=");
      if (i == 0) break;
      if (fsub(mk::keywrd, i + 7, i + 7) == ")") break;
      if (i != 0 && mk::keywrd.find(" ADD-H") == std::string::npos) {
        if (mk::moperr) { delete[] ioptl; return; }
        j = index1(mk::keywrd.substr(i - 1), ") ") + i;
        mk::allkey = mk::keywrd.substr(i - 1, j - i + 1);
        for (k = 1; k <= 10; ++k) neutral[k] = false;
        if (index1(mk::keywrd.substr(i - 1, j - i), "\"") == 0) {
          if (mk::keywrd.substr(i - 1, j - i).find("COOH") != std::string::npos) neutral[1] = true;
          else if (mk::keywrd.substr(i - 1, j - i).find("COO") != std::string::npos) neutral[2] = true;
          if (mk::keywrd.substr(i - 1, j - i).find("NH3") != std::string::npos) neutral[3] = true;
          else if (mk::keywrd.substr(i - 1, j - i).find("NH2") != std::string::npos) neutral[4] = true;
          if (mk::keywrd.substr(i - 1, j - i).find("ARG(+)") != std::string::npos) neutral[5] = true;
          else if (mk::keywrd.substr(i - 1, j - i).find("ARG") != std::string::npos) neutral[6] = true;
          if (mk::keywrd.substr(i - 1, j - i).find("HIS(+)") != std::string::npos) neutral[7] = true;
          else if (mk::keywrd.substr(i - 1, j - i).find("HIS") != std::string::npos) neutral[8] = true;
          if (mk::keywrd.substr(i - 1, j - i).find("SO4") != std::string::npos) neutral[9] = true;
          if (mk::keywrd.substr(i - 1, j - i).find("PO4") != std::string::npos) neutral[10] = true;
        }
        for (k = 1; k <= 10; ++k) {
          if (neutral[k]) break;
        }
        if (k > 10) {
          if (mk::keywrd.substr(i - 1, j - i).find("SALT") != std::string::npos) {
            find_salt_bridges(numbon.data(), numbon.data(), 0, 0);
            if (mk::moperr) { delete[] ioptl; return; }
            if (mk::keywrd.substr(i - 1, j - i).find("SALT") != std::string::npos) {
              i = index1(mk::keywrd.substr(i - 1, j - i), "SALT") + i - 2;
              for (j = i + 1; j <= (int)len_trim(mk::keywrd); ++j) {
                if (fsub(mk::keywrd, j, j + 1) == ") " || fsub(mk::keywrd, j, j) == ",") break;
              }
              mk::keywrd = mk::keywrd.substr(0, i) + mk::keywrd.substr(j - 1);
              continue;
            }
          }
          // Now check for specific residues.
          i = index1(mk::keywrd, " SITE=");
          j = index1(mk::keywrd.substr(i - 1), ") ") + i;
          for (i = i; i <= (int)len_trim(mk::keywrd); ++i) {
            if (fsub(mk::keywrd, i, i) == "(") break;
          }
          j = j - 1;
          m = 0;
          for (int r = 1; r <= max_sites; ++r)
            for (int c = 1; c <= 3; ++c) charge[r][c] = '*';
          while (true) {
            i = i + 1;
            while (true) {
              if (fsub(mk::keywrd, i, i) != "\"") break;
              for (i = i + 1; i <= (int)len_trim(mk::keywrd); ++i) {
                if (fsub(mk::keywrd, i, i) == ")") break;
              }
              i = i + 2;
            }
            if (i >= j) break;
            m = m + 1;
            new_chain[m] = fsub(mk::keywrd, i, i)[0];
            num = std::string(1, (char)('1' + (int)std::floor(std::log10(m * 1.0))));
            if (new_chain[m] < 'A' || new_chain[m] > 'Z') {
              std::fprintf(stdout, "\n%10s\n%10s%s\n",
                           " There is a fault in the SITE keyword",
                           " Chain letter for residue", num.c_str());
            }
            i = i + 1;
            k = (int)((unsigned char)fsub(mk::keywrd, i, i)[0] - (unsigned char)'0');
            if ((k < 0 || k > 9) && fsub(mk::keywrd, i, i) != "-") {
              std::fprintf(stdout, "\n%10s\n%10s%s\n",
                           " There is a fault in the SITE keyword",
                           " The entry for site", num.c_str());
            }
            if (fsub(mk::keywrd, i, i) == "-") {
              ii = -1;
              k = 0;
            } else {
              ii = 1;
            }
            while (true) {
              i = i + 1;
              if (fsub(mk::keywrd, i, i) < "0" || fsub(mk::keywrd, i, i) > "9") break;
              k = k * 10 + (int)((unsigned char)fsub(mk::keywrd, i, i)[0] - (unsigned char)'0');
            }
            new_res[m] = k * ii;
            i = i + 1;
            charge[m][1] = fsub(mk::keywrd, i, i)[0];
            if (charge[m][1] != '+' && charge[m][1] != '-' && charge[m][1] != '0') {
              std::fprintf(stdout, "\n%10s\n%10s%4d%s\n",
                           " There is a fault in  the SITE keyword",
                           " Charge on site", m, " is not '+', '0', or '-'");
            }
            i = i + 1;
            charge[m][2] = fsub(mk::keywrd, i, i)[0];
            if (charge[m][2] == '+' || charge[m][2] == '-' || charge[m][2] == '0') {
              i = i + 1;
            } else {
              charge[m][2] = '*';
            }
            i = i + 1;
          }
        } else {
          m = 0;
        }
        lsite = true;
        update_txtatm(true, true);
        i = mk::numat;
        site(neutral, new_chain, new_res, charge, m, max_sites, mk::allkey);
        for (ii = i; ii <= mk::numat; ++ii) mc::l_atom[ii] = true;
        update_txtatm(true, true);
        if (mk::moperr) { delete[] ioptl; return; }
      }
      i = index1(mk::keywrd, " SITE=");
      if (i == 0) break;
      for (j = i + 1; j <= (int)len_trim(mk::keywrd); ++j) {
        if (fsub(mk::keywrd, j, j + 1) == ") ") break;
      }
      mk::keywrd = mk::keywrd.substr(0, i) + mk::keywrd.substr(j);
    }
  }

  // ---- Store charge, if present ------------------------------------------
  for (i = 1; i <= mk::natoms; ++i) {
    atom_charge[i] = fsub(mc::txtatm[i], 2, 2)[0];
    // Prevent atom number being mis-read as a charge.
    if (fsub(mc::txtatm[i], 2, 2) != "+" && fsub(mc::txtatm[i], 2, 2) != "-" &&
        fsub(mc::txtatm[i], 2, 2) != "0")
      atom_charge[i] = ' ';
  }

  // ---- Assign logicals using keywrd --------------------------------------
  lres = (index1(mk::keywrd, " RESI") + index1(mk::keywrd, " NEWGEO") +
          index1(mk::keywrd, " RESEQ") != 0);
  if (!lres) lres = (index1(mk::keywrd, " PDBOUT") != 0 && mk::maxtxt != 26);
  if (!lres) lres = (index1(mk::keywrd, " ADD-H") != 0 && index1(mk::keywrd, " NORESEQ") == 0);
  lreseq = (index1(mk::keywrd, " NORESEQ") == 0 && index1(mk::keywrd, " RESEQ") != 0);
  if (lreseq && mk::maxtxt != 26 && index1(mk::keywrd, "RESID") == 0) {
    mk::line = "RESEQ only works when the atom labels are in PDB format";
    mopend(mk::line.substr(0, len_trim(mk::line)));
    std::fprintf(stdout, "%s\n", "(Before using RESEQ, run a job using keyword RESIDUES to add PDB atom labels.)");
    delete[] ioptl;
    return;
  }
  let = (index1(mk::keywrd, " 0SCF") + index1(mk::keywrd, " LET") +
         index1(mk::keywrd, " RESEQ") + index1(mk::keywrd, " GEO-OK") != 0);
  times = (index1(mk::keywrd, " TIMES") != 0);
  if (times) timer(" START OF GEOCHK");
  debug = (index1(mk::keywrd, " GEOCHK") != 0);
  if (index1(mk::keywrd, " CHARGE=") != 0) {
    irefq = nint(reada(mk::keywrd, index1(mk::keywrd, " CHARGE=")));
  } else {
    irefq = 0;
  }
  extvdw_for_MOZYME(radius, mr::atom_radius_covalent);
  if (mk::moperr) { delete[] ioptl; return; }

  int large;
  if (index1(mk::keywrd, " LARGE") != 0) {
    large = 1000000;
  } else {
    large = 20;
  }

  // ---- Work out what atoms are bonded to each other ------------------------
  lewis(true);
  if (mk::moperr) { delete[] ioptl; return; }
  nnbonds.assign(mk::numat + 1, 0);
  iibonds.assign(16, std::vector<int>(mk::numat + 1, 0));
  for (i = 1; i <= mk::numat; ++i) nnbonds[i] = mc::nbonds[i];
  for (i = 1; i <= mk::numat; ++i)
    for (j = 1; j <= 15; ++j) iibonds[j][i] = mc::ibonds[j][i];
  if (mk::moperr) {
    if (index1(mk::keywrd, " GEO-OK") == 0) {
      std::fprintf(stdout, " GEOMETRY CONTAINS FAULTS. TO CONTINUE CALCULATION SPECIFY \"GEO-OK\"\n");
      goto label_1100;
    } else {
      mk::moperr = false;
    }
  }

  // ---- Zero out ions ------------------------------------------------------
  for (i = 1; i <= mk::numat; ++i) ioptl[i] = false;
  findn1(n1, ioptl, io);
  if (lreseq) {
    // ---- Resequence the atoms ---------------------------------------------
    // First, delete all bonds between ATOMs and HETATMs.
    for (i = 1; i <= mk::numat; ++i) {
      if (fsub(mc::txtatm[i], 1, 4) != "ATOM") continue;
      for (j = 1; j <= mc::nbonds[i]; ++j) {
        k = mc::ibonds[j][i];
        if (fsub(mc::txtatm[k], 1, 4) != "ATOM") mc::ibonds[j][i] = 0;
      }
      l = 0;
      for (j = 1; j <= mc::nbonds[i]; ++j) {
        k = mc::ibonds[j][i];
        if (k > 0) {
          l = l + 1;
          mc::ibonds[l][i] = k;
        }
      }
      mc::nbonds[i] = l;
    }
    if (index1(mk::keywrd, "RESID") != 0)
      for (i = 1; i <= mk::numat; ++i) mc::txtatm1[i] = " ";
    new_ = 0;
    for (i = 1; i <= mk::maxatoms; ++i) mz::iz[i] = -1000;
    while (true) {
      if (n1 != 0) reseq(ioptl, mz::iz.data(), n1, new_, io);
      if (mk::moperr) goto label_1100;
      findn1(n1, ioptl, io);
      if (n1 == 0) break;
    }
    if (new_ != mk::numat) {
      for (i = 1; i <= mk::numat; ++i) {
        if (mc::nat[i] != 1 && !ioptl[i]) {
          // Identify all non-protein molecules in the system.
          std::vector<bool> vb(mk::maxatoms + 1, false);
          for (ii = 1; ii <= mk::numat; ++ii) vb[ii] = ioptl[ii];
          moiety(vb, mz::iz, i, new_);
          for (ii = 1; ii <= mk::numat; ++ii) ioptl[ii] = vb[ii];
        }
      }
      for (i = 1; i <= mk::numat; ++i) {
        if (mc::nat[i] == 1 && !ioptl[i]) {
          // Identify all hydrogens attached to residue-like species.
          new_ = new_ + 1;
          mz::iz[new_] = i;
        }
      }
      if (new_ != mk::numat) {
        std::fprintf(stdout, " THERE IS A FAULT IN RESEQ\n");
        std::fprintf(stdout, "  Number of atoms found in data-set:  %5d\n", mk::numat);
        std::fprintf(stdout, "  Number of atoms after re-sequencing:%5d\n", new_);
        if (new_ < mk::numat) {
          std::fprintf(stdout, "Atoms missing (Use original numbering system)\n");
          for (i = 1; i <= mk::numat; ++i) ioptl[i] = true;
          for (i = 1; i <= new_; ++i) ioptl[mz::iz[i]] = false;
          for (i = 1; i <= mk::numat; ++i) {
            if (ioptl[i]) std::fprintf(stdout, "%5d\n", i);
          }
        }
        mopend("THERE IS A FAULT IN RESEQ");
        goto label_1100;
      }
      // Unconditionally, convert geometry into Cartesian coordinates.
      for (i = 1; i <= mk::numat; ++i)
        for (k = 1; k <= 3; ++k) mc::geo[k][i] = mc::coord[k-1][i];
      for (i = 1; i <= mk::numat; ++i) {
        if (mc::na[i] != 0) break;
      }
      if (i <= mk::numat) {
        mopend("SOME COORDINATES WERE IN INTERNAL. THESE HAVE BEEN CHANGED TO CARTESIAN");
        mk::moperr = false;
      }
      for (i = 1; i <= mk::maxatoms; ++i) mc::na[i] = 0;
    }
    l = 1;
    for (i = 1; i <= mk::numat; ++i) {
      mc::nfirst[i] = l;
      j = mz::iz[i];
      mb[j] = i;
      for (k = 1; k <= 3; ++k) mc::geo[k][i] = mc::coord[k-1][j];
      mc::labels[i] = mc::nat[j];
      mc::nlast[i] = mc::nfirst[i] + pc::natorb[mc::labels[i]] - 1;
      l = mc::nlast[i] + 1;
    }
    for (i = 1; i <= mk::numat; ++i) mc::nat[i] = mc::labels[i];
    for (i = 1; i <= mk::numat; ++i)
      for (k = 1; k <= 3; ++k) mc::coord[k-1][i] = mc::geo[k][i];
    done = true;
    mk::natoms = mk::numat;

    // Rearrange atoms to suit the new numbering system.
    for (i = 1; i <= mk::numat; ++i) {
      j = mz::iz[i];
      mz::ib[i] = (mc::nbonds[j] < 4) ? mc::nbonds[j] : 4;
      l = (mc::nbonds[i] < 4) ? mc::nbonds[i] : 4;
      for (k = 1; k <= l; ++k) {
        mc::ibonds[k + 4][i] = mb[mc::ibonds[k][i]];
      }
    }
    for (i = 1; i <= mk::numat; ++i) {
      mc::nbonds[i] = mz::ib[i];
      for (k = 1; k <= mc::nbonds[i]; ++k) {
        mc::ibonds[k][i] = mc::ibonds[k + 4][mz::iz[i]];
      }
    }
  }
  mz::noccupied = 0;
  if (index1(mk::keywrd, " RESEQ") + index1(mk::keywrd, " SITE=") +
          index1(mk::keywrd, " ADD-H") + index1(mk::keywrd, " 0SCF") == 0) {
    // ---- Examine the geometry: identify the Lewis elements ----------------
    for (i = 1; i <= 3; ++i) numbon[i] = 0;
    chklew(mb, numbon, l, large, debug);
    if (mk::moperr) { delete[] ioptl; goto label_1100; }
    l = 0;
    for (i = 1; i <= mk::numat; ++i) l = l + std::abs(mz::iz[i]);
    if (l != 0) {
      // There are ions.  Identify them.
      chkion(mb, numbon[2], atom_charge);
      if (lreseq) mk::moperr = false;
    }
    for (i = 1; i <= mk::numat; ++i) mz::ions[i] = nint(pc::tore[mc::nat[i]]);
    for (i = 1; i <= mz::Lewis_tot; ++i) {
      if (mz::Lewis_elem[1][i] > 0) {
        mz::noccupied = mz::noccupied + 1;
        j = mz::Lewis_elem[1][i];
        if (mz::Lewis_elem[2][i] > 0) {
          k = mz::Lewis_elem[2][i];
          mz::ions[k] = mz::ions[k] - 1;  // one electron from a bond
          mz::ions[j] = mz::ions[j] - 1;  // one electron from a bond
        } else {
          mz::ions[j] = mz::ions[j] - 2;  // two electrons from a lone pair
        }
      }
    }
    mz::nvirtual = 0;
    for (i = 1; i <= mz::Lewis_tot; ++i) {
      if (mz::Lewis_elem[2][i] > 0) mz::nvirtual = mz::nvirtual + 1;
    }
  } else {
    mk::nvar = 0;
    for (i = 1; i <= mk::numat; ++i) {
      for (j = 1; j <= mk::numat_old; ++j) {
        if (std::abs(mc::coord[0][i] - mc::coorda[0][j]) > 0.1) continue;
        if (std::abs(mc::coord[1][i] - mc::coorda[1][j]) > 0.1) continue;
        if (std::abs(mc::coord[2][i] - mc::coorda[2][j]) > 0.1) continue;
        break;
      }
      if (j <= mk::numat_old) {
        for (k = 1; k <= 3; ++k) {
          if (mc::lopt[k][j] == 1) {
            mk::nvar = mk::nvar + 1;
            mc::loc[1][mk::nvar] = i;
            mc::loc[2][mk::nvar] = k;
          }
        }
      } else {
        for (k = 1; k <= 3; ++k) {
          mk::nvar = mk::nvar + 1;
          mc::loc[1][mk::nvar] = i;
          mc::loc[2][mk::nvar] = k;
        }
      }
    }
  }

  // ---- RESIDUE handling (lres) -------------------------------------------
  if (lres) {
    for (i = 1; i <= mk::maxatoms; ++i) mc::txtatm[i] = " ";
    for (i = 1; i <= mk::maxatoms; ++i) mz::angles[i].assign(4, 0.0);
    for (i = 1; i <= mk::maxatoms; ++i) mz::allres[i] = " ";
    if (done) {
      for (i = 1; i <= mk::numat; ++i) ioptl[i] = false;
      findn1(n1, ioptl, io);
    }
    l_protein = (n1 != 0);
    for (i = 1; i <= mk::numat; ++i) mz::ib[i] = -100000;
    nfrag = 0;
    ires = 0;
    uni_res = 0;
    mz::odd_h = true;
    // Break all intra-chain bonds, so that the residues can easily be identified.
    lyse();
    for (i = 1; i <= mk::maxatoms; ++i) mz::allr[i] = " ";
    while (true) {
      nfrag = nfrag + 1;
      if (mz::start_res[nfrag] != -200) ires = mz::start_res[nfrag];
      if (nfrag > 99) {
        mopend("STRUCTURE UNRECOGNIZABLE");
        opend = false;  // inquire(unit=iarc, opened=opend): simplified
        goto label_1100;
      }
      names(ioptl, mz::ib.data(), n1, ires, nfrag, io, uni_res, mres);
      if (mk::moperr) { delete[] ioptl; return; }
      nn1 = n1;
      findn1(n1, ioptl, io);
      if (!l_protein) nfrag = 0;
      if (n1 == 0) break;
      if (n1 == nn1) ioptl[n1] = true;
      if (!lbreaks) {
        mbreaks = mbreaks + 1;
        mc::breaks[mbreaks] = n1;
        // Find the last atom that has been defined.
        for (i = 1; i <= mk::numat; ++i) {
          if (!ioptl[i]) break;
        }
        if (i - 1 > 0)
          for (k = 1; k <= 3; ++k) mc::break_coords[k - 1][mbreaks] = mc::coord[k - 1][i - 1];
      }
    }
    // Re-evaluate all residues.
    j = 1;
    mz::allres[j] = fsub(mc::txtatm[1], 18, 20);
    for (i = 2; i <= mk::natoms; ++i) {
      if (is_blank(mc::txtatm[i])) break;
      if (mc::nat[i] != 1 && fsub(mc::txtatm[i], 23, 26) != fsub(mc::txtatm[i - 1], 23, 26)) {
        j = j + 1;
        mz::allres[j] = fsub(mc::txtatm[i], 18, 20);
      }
    }
    ires = j;
    for (i = 1; i <= mk::natoms; ++i) iopt[i] = mz::ib[i];
    // Label the atoms in any non-protein molecules in the system.
    nfrag = nfrag + 1;
    if (mz::start_res[nfrag] == -200) {
      if (!l_protein) ires = 0;
    } else {
      ires = mz::start_res[nfrag];
    }
    ligand(ires, mz::start_res, nfrag);
    // If ligands are present, set a break at the end of the protein.
    if (!lbreaks) {
      mbreaks = mbreaks + 1;
      mc::breaks[mbreaks] = n1;
      for (i = 1; i <= mk::numat; ++i) {
        if (!ioptl[i]) break;
      }
      if (i - 1 > 0)
        for (k = 1; k <= 3; ++k) mc::break_coords[k - 1][mbreaks] = mc::coord[k - 1][i - 1];
    }
    mz::nres = uni_res;
    mk::maxtxt = 26;
    // Add chain letters.
    reset_breaks();
    i = index1(mk::keywrd, " RESI");
    if (i > 0) {
      j = index1(mk::keywrd.substr(i), " ") + i;
      j = index1(mk::keywrd.substr(i, j - i), "0");
    }
    if (i > 0 && j > 0) {
      for (i = 1; i <= mk::numat; ++i) {
        setsub(mc::txtatm[i], 13, 16, fsub(mc::txtatm1[i], 13, 16));
      }
    }
    if (mk::moperr) { delete[] ioptl; return; }
    mbreaks = 1;
    mk::line = fsub(mc::txtatm[1], 23, 26);
    for (i = 1; i <= mk::numat; ++i) {
      setsub(mc::txtatm[i], 22, 22, std::string(1, mc::chains[mbreaks - 1]));
      if (i == mc::breaks[mbreaks]) {
        mbreaks = mbreaks + 1;
      } else {
        if (fsub(mc::txtatm[i], 23, 26) == "    ") continue;
        if (fsub(mc::txtatm[i], 1, 6) == "HETATM") {
          if (fsub(mc::txtatm[i], 23, 26) != mk::line.substr(0, 4)) {
            if (fsub(mc::txtatm[i - 1], 1, 4) != "ATOM") mbreaks = (26 < mbreaks + 1) ? 26 : (mbreaks + 1);
          }
          setsub(mc::txtatm[i], 22, 22, std::string(1, mc::chains[mbreaks - 1]));
        }
      }
      mk::line = fsub(mc::txtatm[i], 23, 26);
    }
    // Check for unknowns.
    for (i = 1; i <= mk::numat; ++i) {
      for (j = 1; j <= 20; ++j) {
        if (fsub(mc::txtatm[i], 18, 20) == mz::tyres[j]) break;
      }
      if (j == 21) {
        if (fsub(mc::txtatm[i], 15, 16) != "  ") continue;
        j = 1;
        for (k = i + 1; k <= mk::numat; ++k) {
          if (fsub(mc::txtatm[k], 14, 26) == fsub(mc::txtatm[i], 14, 26)) {
            j = j + 1;
            if (j < 10) {
              setsub(mc::txtatm[k], 15, 15, std::string(1, (char)('0' + j)));
            } else {
              char b[4];
              std::snprintf(b, sizeof b, "%2d", j);
              setsub(mc::txtatm[k], 15, 16, b);
            }
          }
        }
        if (j > 1) setsub(mc::txtatm[i], 15, 15, "1");
      }
    }
    if (n_new != 0) {
      // Re-name residues to use the XENO name.
      mbreaks = 1;
      for (i = 1; i <= n_new; ++i) old_name[i] = "---";
      for (i = 1; i <= max_sites; ++i) l_names[i] = false;
      for (i = 1; i <= mk::natoms; ++i) {
        ter = (i == mc::breaks[mbreaks]);
        if (ter) mbreaks = mbreaks + 1;
        j = nint(reada(mc::txtatm[i].substr(22), 1));
        if (fsub(mc::txtatm[i], 1, 6) == "HETATM") continue;
        for (k = 1; k <= n_new; ++k) {
          if (new_chain[k] == mc::chains[mbreaks - 1] && new_res[k] == j) {
            if (old_name[k] == "---") old_name[k] = fsub(mc::txtatm[i], 18, 20);
            setsub(mc::txtatm[i], 18, 20, new_name[k]);
            l_names[k] = true;
            break;
          }
        }
      }
      first = true;
      for (i = 1; i <= n_new; ++i) {
        if (l_names[i]) {
          if (old_name[i] == new_name[i]) {
            if (first) {
              std::fprintf(stdout, "\n%s\n%s\n",
                           "      Residue names that have been changed",
                           "      Residue No.  Calc'd name   XENO name");
              first = false;
            }
            std::fprintf(stdout, "%3d%10d%2s%c%7s%-3s%9s%-3s%s\n", i, new_res[i], "", new_chain[i],
                         "", old_name[i].c_str(), "", new_name[i].c_str(),
                         "  Name not changed!");
          } else {
            std::fprintf(stdout, "%3d%10d%2s%c%7s%-3s%9s%-3s\n", i, new_res[i], "", new_chain[i],
                         "", old_name[i].c_str(), "", new_name[i].c_str());
          }
        }
      }
    }
    bool l_use_old_labels =
        (index1(mk::keywrd, " SITE=") != 0 && index1(mk::keywrd, " ADD-H") == 0);
    l_use_old_labels = true;
    update_txtatm(l_use_old_labels, true);
    write_sequence();
    if (index1(mk::keywrd, " RAMA") != 0) {
      if (index1(mk::keywrd, " ADD-H") == 0 && uni_res > 1) {
        std::fprintf(stdout, "\n%10s\n%10s\n", "", "        Ramachandran Angles");
        std::fprintf(stdout, "%10s\n", "    Residue    phi    psi  omega");
      }
      for (i = 1; i <= uni_res; ++i) {
        if (std::abs(mz::angles[1][i]) + std::abs(mz::angles[3][i]) > 1.e-20 &&
            mz::res_start[i] > 0) {
          std::fprintf(stdout, "%14s%7.1f%7.1f%7.1f\n",
                       (fsub(mc::txtatm[mz::res_start[i]], 18, 20) +
                        fsub(mc::txtatm[mz::res_start[i]], 23, 26)).c_str(),
                       mz::angles[1][i], mz::angles[2][i], mz::angles[3][i]);
        }
      }
      std::fprintf(stdout, "\n");
    }
  }
  if (index1(mk::keywrd, " PDBOUT") != 0) {
    // Identify atoms where chain breaks occur.
    if (index1(mk::keywrd, " RESEQ") + index1(mk::keywrd, " ADD-H") != 0) reset_breaks();
    if (index1(mk::keywrd, " RESEQ") + index1(mk::keywrd, " ADD-H") == 0) {
      if (index1(mk::keywrd, " RESID") != 0) {
        if (!mc::txtatm1.empty() && !is_blank(mc::txtatm1[1])) compare_sequence(n_new);
      }
    }
  }

  // ---- Edit keywords to remove text not used in the next calculation ------
  for (i = 1; i <= 6; ++i) {
    mk::line = mk::refkey[i];
    upcase(mk::line, (int)len_trim(mk::line));
    j = index1(mk::line, " SITE=");
    if (j > 0) {
      k = index1(mk::line.substr(j - 1), ") ") + j;
      mk::refkey[i] = mk::refkey[i].substr(0, j) + mk::refkey[i].substr(k - 1);
    }
    j = index1(mk::line, " RESEQ");
    if (j != 0) mk::refkey[i] = mk::refkey[i].substr(0, j - 1) +
                                mk::refkey[i].substr(j + 5);
  }
  mk::line = " " + mk::refkey[1];
  upcase(mk::line, (int)len_trim(mk::line));
  i = index1(mk::line, " RESEQ");
  if (i != 0) mk::refkey[1] = mk::refkey[1].substr(0, i - 1) + mk::refkey[1].substr(i + 4);
  if (index1(mk::keywrd, " ADD-H") != 0) { delete[] ioptl; return; }

  // ---- Modify ions so that it refers to all atoms (real and dummy) --------
  j = 0;
  for (i = 1; i <= mk::maxatoms; ++i) mz::iz[i] = mz::ions[i];
  for (i = 1; i <= mk::maxatoms; ++i) mz::ions[i] = 0;
  for (i = 1; i <= mk::natoms; ++i) {
    if (mc::labels[i] == 99) {
      mz::ions[i] = 0;
    } else {
      j = j + 1;
      mz::ions[i] = mz::iz[j];
    }
  }
  ibad = 0;
  if (lres && !lreseq) {
    // Check all ions to see if any residue is a di-ion.
    for (i = 1; i <= mk::numat; ++i) {
      if (mz::ions[i] != 1 && mz::ions[i] != -1) continue;  // WARNING (F90 no-op loop)
    }
    if (ibad != 0) std::fprintf(stdout, "\n");
    jbad = ibad;
    ibad = 0;
    ibad = ibad + jbad;
  }

  charges = (lsite || index1(mk::keywrd, "CHARGES") != 0);
  ichrge = -mz::noccupied * 2;
  for (i = 1; i <= mk::numat; ++i) ichrge = ichrge + nint(pc::tore[mc::nat[i]]);
  mk::line = " ";
  if (mz::noccupied != 0 && (index1(mk::keywrd, " LEWIS") > 0 || mz::noccupied * 2 != mk::nelecs)) {
    mk::maxtxt = 0;
    for (i = 1; i <= mk::numat; ++i) {
      mk::maxtxt = (mk::maxtxt > (int)len_trim(mc::txtatm[i])) ? mk::maxtxt : (int)len_trim(mc::txtatm[i]);
    }
    if (mk::maxtxt == 0) {
      j = 1;
    } else {
      j = mk::maxtxt / 2 + 2;
    }
    update_txtatm(true, false);
    if (mk::prt_topo) {
      std::fprintf(stdout, "\n%s\n", "   TOPOGRAPHY OF SYSTEM");
      std::fprintf(stdout, "%s%s%s%s\n", "  ATOM No. ", mk::line.substr(0, j).c_str(),
                   "  LABEL  ", mk::line.substr(0, j).c_str());
      std::fprintf(stdout, "%s\n", "Atoms connected to this atom");
      if (j == 0) {
        for (i = 1; i <= mk::numat; ++i) {
          std::fprintf(stdout, "%7d%9s", i, (ec::elemnt[mc::nat[i]] + "  ").c_str());
          for (j = 1; j <= mc::nbonds[i]; ++j) std::fprintf(stdout, "%7d", mc::ibonds[j][i]);
          std::fprintf(stdout, "\n");
        }
      } else {
        if (mk::maxtxt > 2) {
          for (i = 1; i <= mk::numat; ++i) {
            std::fprintf(stdout, "%7d%9s", i, (ec::elemnt[mc::nat[i]] + " (" +
                          mc::txtatm[i].substr(0, mk::maxtxt) + ") ").c_str());
            for (j = 1; j <= mc::nbonds[i]; ++j) std::fprintf(stdout, "%7d", mc::ibonds[j][i]);
            std::fprintf(stdout, "\n");
          }
        } else {
          for (i = 1; i <= mk::numat; ++i) {
            std::fprintf(stdout, "%7d%9s", i, ec::elemnt[mc::nat[i]].c_str());
            for (j = 1; j <= mc::nbonds[i]; ++j) std::fprintf(stdout, "%7d", mc::ibonds[j][i]);
            std::fprintf(stdout, "\n");
          }
        }
      }
    }
  }
  if (mz::noccupied != 0 && index1(mk::keywrd, " LEWIS") > 0) {
    std::fprintf(stdout, "\n%37s\n", "   Lewis Structure");
    if (index1(mk::keywrd, " LARGE") != 0)
      std::fprintf(stdout, "\n%23s\n", "  ATOMS IN OCCUPIED LOCALIZED MOLECULAR ORBITALS");
    int l4 = 4;
    int lew_rows = mz::Lewis_tot / l4 + 5;
    Lewis_formatted.assign(lew_rows + 1, std::vector<std::string>(l4 + 1, " "));
    int kc = 0;
    int jc = 0;
    for (i = 1; i <= mz::Lewis_tot; ++i) {
      if (mz::Lewis_elem[1][i] > 0) jc = jc + 1;
    }
    int mm = jc / l4 + 1;
    std::fprintf(stdout, "%s\n", "     LMO  Atom  Atom        LMO  Atom  Atom        LMO  Atom  Atom        LMO  Atom  Atom    ");
    ii = 0;
    jj = 1;
    for (i = 1; i <= mz::Lewis_tot; ++i) {
      if (mz::Lewis_elem[1][i] > 0) {
        kc = kc + 1;
        ii = ii + 1;
        if (ii > mm) {
          ii = 1;
          jj = jj + 1;
        }
        char b[24];
        if (mz::Lewis_elem[2][i] > 0) {
          std::snprintf(b, sizeof b, "%8d%6d%6d", kc, mz::Lewis_elem[1][i], mz::Lewis_elem[2][i]);
        } else {
          std::snprintf(b, sizeof b, "%8d%6d", kc, mz::Lewis_elem[1][i]);
        }
        Lewis_formatted[ii][jj] = b;
      }
    }
    for (i = 1; i <= lew_rows; ++i) {
      if (is_blank(Lewis_formatted[i][1])) break;
      std::fprintf(stdout, "%s    %s    %s    %s\n",
                   Lewis_formatted[i][1].c_str(), Lewis_formatted[i][2].c_str(),
                   Lewis_formatted[i][3].c_str(), Lewis_formatted[i][4].c_str());
    }
    if (index1(mk::keywrd, " LARGE") != 0) {
      std::fprintf(stdout, "\n%23s\n", "  ATOMS IN UNOCCUPIED LOCALIZED MOLECULAR ORBITALS");
      Lewis_formatted.assign(lew_rows + 1, std::vector<std::string>(l4 + 1, " "));
      jc = 0;
      for (i = 1; i <= mz::Lewis_tot; ++i) {
        if (mz::Lewis_elem[2][i] > 0) jc = jc + 1;
      }
      mm = jc / l4 + 1;
      kc = 0;
      ii = 0;
      jj = 1;
      for (i = 1; i <= mz::Lewis_tot; ++i) {
        if (mz::Lewis_elem[2][i] > 0) {
          kc = kc + 1;
          ii = ii + 1;
          if (ii > mm) {
            ii = 1;
            jj = jj + 1;
          }
          char b[24];
          if (mz::Lewis_elem[1][i] > 0) {
            std::snprintf(b, sizeof b, "%8d%6d%6d", kc, mz::Lewis_elem[1][i], mz::Lewis_elem[2][i]);
          } else {
            std::snprintf(b, sizeof b, "%8d%6d", kc, mz::Lewis_elem[2][i]);
          }
          Lewis_formatted[ii][jj] = b;
        }
      }
      for (i = 1; i <= lew_rows; ++i) {
        if (is_blank(Lewis_formatted[i][1])) break;
        std::fprintf(stdout, "%s    %s    %s    %s\n",
                     Lewis_formatted[i][1].c_str(), Lewis_formatted[i][2].c_str(),
                     Lewis_formatted[i][3].c_str(), Lewis_formatted[i][4].c_str());
      }
    }
    std::fprintf(stdout, "\n");
    if (mz::noccupied > 0) {
      std::fprintf(stdout, "%s\n\n", "          Type          Number of Lewis structural elements identified");
    }
    if (numbon[1] != 0) std::fprintf(stdout, "%s%6d\n", "         SIGMA BONDS   ", numbon[1]);
    if (numbon[2] != 0) std::fprintf(stdout, "%s%6d\n", "         LONE PAIRS    ", numbon[2]);
    if (numbon[3] != 0) std::fprintf(stdout, "%s%6d\n", "         PI BONDS      ", numbon[3]);
    if (index1(mk::keywrd, " LEWIS") > 0 || (mz::noccupied * 2 != mk::nelecs && mz::noccupied != 0)) {
      std::fprintf(stdout, "\n%s%6d\n", " Number of filled levels from atoms and charge:", mk::nelecs / 2);
      std::fprintf(stdout, "%s%6d\n", " Number of filled levels from Lewis structure: ", mz::noccupied);
    }
    l = 0;
    m = 0;
    num = std::string(1, (char)('2' + (int)std::floor(std::log10(mk::numat + 1.0))));
    for (i = 1; i <= mk::numat; ++i) {
      if (!pc::main_group[mc::nat[i]]) {
        // Element is a transition metal.  Work out its formal oxidation state.
        k = mz::ions[i];
        for (j = 1; j <= mz::Lewis_tot; ++j) {
          if (mz::Lewis_elem[1][j] != 0 && mz::Lewis_elem[2][j] != 0) {
            if (mz::Lewis_elem[1][j] == i) k = k + 1;
            if (mz::Lewis_elem[2][j] == i) k = k + 1;
          }
        }
        if (m == 0) std::fprintf(stdout, "\n");
        m = 1;
        std::fprintf(stdout, "%10s%4d%s%3d\n", " Formal oxidation state of atom", i,
                     (", a " + ec::elemnt[mc::nat[i]] + ", is").c_str(), k);
        if (k < 0) l = 1;
        if (k > 3) l = 1;
      }
    }
    if (l == 1) web_message(ch::iw, "Lewis_structures.html");
  }

  // ---- Check for sulfate and phosphate ------------------------------------
  for (i = 1; i <= mk::numat; ++i) {
    if (mc::nat[i] == 16 && fsub(mc::txtatm[i], 18, 20) == "SO4") {
      k = 2;
      mz::ions[i] = 0;
      for (j = 1; j <= 4; ++j) {
        l = mc::ibonds[j][i];
        if (mz::ions[l] == -1) {
          mz::ions[l] = 0;
          k = k - 1;
          if (k == 0) break;
        }
      }
    }
    if (mc::nat[i] == 15 && fsub(mc::txtatm[i], 18, 20) == "PO4") {
      k = 1;
      mz::ions[i] = 0;
      for (j = 1; j <= 4; ++j) {
        l = mc::ibonds[j][i];
        if (mz::ions[l] == -1) {
          mz::ions[l] = 0;
          k = k - 1;
          if (k == 0) break;
        }
      }
    }
  }
  for (i = -6; i <= 6; ++i) num_ions[i + 6] = 0;
  for (i = 1; i <= mk::numat; ++i) {
    j = (mz::ions[i] < -6) ? -6 : ((mz::ions[i] > 6) ? 6 : mz::ions[i]);
    num_ions[j + 6] = num_ions[j + 6] + 1;
  }
  i = 0;
  for (j = 1; j <= 6; ++j) i = i + num_ions[j + 6] + num_ions[-j + 6];
  if (i > 0) {
    if (index1(mk::keywrd, " LEWIS") > 0) {
      std::fprintf(stdout, "\n%s\n", "          Type           Number of charged sites identified");
      for (i = 1; i <= 6; ++i) {
        if (num_ions[i + 6] > 0)
          std::fprintf(stdout, "%9s%2x%6d\n", ion_names[i + 6], "", num_ions[i + 6]);
        if (num_ions[-i + 6] > 0)
          std::fprintf(stdout, "%9s%2x%6d\n", ion_names[-i + 6], "", num_ions[-i + 6]);
      }
      i = 0;
      for (j = 1; j <= 6; ++j) i = i + num_ions[j + 6] * j;
      std::fprintf(stdout, "\n%s%5d\n", " SUM OF POSITIVE CHARGES", i);
      i = 0;
      for (j = 1; j <= 6; ++j) i = i + num_ions[-j + 6] * j;
      std::fprintf(stdout, "%s%5d\n", " SUM OF NEGATIVE CHARGES", -i);
    }
    padding = std::string(40, ' ');
    for (i = 1; i <= mk::numat; ++i) {
      if (mz::ions[i] != 0) {
        if (mc::nat[i] == 15 || mc::nat[i] == 16) {
          kk = 0;
          for (jj = 1; jj <= mc::nbonds[i]; ++jj) {
            ii = mc::ibonds[jj][i];
            if (mc::nat[ii] == 8) {
              if (mz::ions[ii] == -1) {
                // Found a PO4 or SO4.  Neutralize the P-O or S-O Zwitterion.
                mz::ions[i] = mz::ions[i] - 1;
                mz::ions[ii] = mz::ions[ii] + 1;
                kk = 1;
                break;
              }
            }
          }
          if (kk == 1) continue;
        }
      }
    }
    for (i = 1; i <= mk::numat; ++i) {
      if (mz::ions[i] != 0) break;
    }
    maxtxt_store = mk::maxtxt;
    if (mk::maxtxt < 0) mk::maxtxt = 14;
    if (i <= mk::numat) {
      if (mk::maxtxt > 1) {
        mk::line = "   Ion Atom No.           Label               Charge";
        std::fprintf(stdout, "\n%s\n", mk::line.substr(0, len_trim(mk::line)).c_str());
        l = (17 - mk::maxtxt / 2 > 1) ? (17 - mk::maxtxt / 2) : 1;
        residues = (index1(mk::keywrd, " RESID") != 0);
      } else {
        std::fprintf(stdout, "\n%s\n", "     Ion Atom No.  Type    Charge");
      }
    } else {
      std::fprintf(stdout, "\n%18s\n", "NO CHARGED ATOMS FOUND.");
    }
    for (j = 6; j >= -6; --j) {
      if (j == 0) continue;
      k = 0;
      for (i = 1; i <= mk::numat; ++i) {
        if (mz::ions[i] == j) {
          if (k == 0) std::fprintf(stdout, "\n");
          k = k + 1;
          mk::line = " ";
          jj = 0;
          if (j == 1) {
            for (ii = 1; ii <= mk::numat; ++ii) {
              if (mz::ions[ii] < 0) {
                sum = (mc::coord[0][i] - mc::coord[0][ii]) * (mc::coord[0][i] - mc::coord[0][ii]) +
                      (mc::coord[1][i] - mc::coord[1][ii]) * (mc::coord[1][i] - mc::coord[1][ii]) +
                      (mc::coord[2][i] - mc::coord[2][ii]) * (mc::coord[2][i] - mc::coord[2][ii]);
                if (sum < 99.0) {
                  jj = jj + 1;
                  if (jj > 100) break;
                  near_ions[jj] = ii;
                  r_ions[jj] = std::sqrt(sum);
                  mk::line = " Angstroms from anion";
                }
              }
            }
          } else if (j == -1) {
            for (ii = 1; ii <= mk::numat; ++ii) {
              if (mz::ions[ii] > 0) {
                sum = (mc::coord[0][i] - mc::coord[0][ii]) * (mc::coord[0][i] - mc::coord[0][ii]) +
                      (mc::coord[1][i] - mc::coord[1][ii]) * (mc::coord[1][i] - mc::coord[1][ii]) +
                      (mc::coord[2][i] - mc::coord[2][ii]) * (mc::coord[2][i] - mc::coord[2][ii]);
                if (sum < 99.0) {
                  jj = jj + 1;
                  if (jj > 100) break;
                  near_ions[jj] = ii;
                  r_ions[jj] = std::sqrt(sum);
                  mk::line = " Angstroms from cation";
                }
              }
            }
          }
          if (jj > 1) {
            // Sort near ions into increasing distance.
            jj = (jj < 100) ? jj : 100;
            for (ii = 1; ii <= jj; ++ii) {
              for (kk = ii + 1; kk <= jj; ++kk) {
                if (r_ions[kk] < r_ions[ii]) {
                  sum = r_ions[kk];
                  r_ions[kk] = r_ions[ii];
                  r_ions[ii] = sum;
                  m = near_ions[kk];
                  near_ions[kk] = near_ions[ii];
                  near_ions[ii] = m;
                }
              }
            }
            for (ii = 2; ii <= jj; ++ii) {
              if (r_ions[ii] > 5.0) {
                jj = ii - 1;
                break;
              }
            }
          }
          if (mk::maxtxt > 1) {
            if (jj > 0) {
              if (residues) {
                txtatm_1 = mc::txtatm[i];
                txtatm_2 = mc::txtatm[near_ions[1]];
              } else {
                txtatm_1 = mc::txtatm1[i];
                txtatm_2 = mc::txtatm1[near_ions[1]];
              }
              num = "5";
              if (fsub(ec::elemnt[mc::nat[near_ions[1]]], 1, 1) != " ") num = "5";
              m = (int)len_trim(mk::line) + 1;
              if (fsub(txtatm_2, mk::maxtxt, mk::maxtxt) == ")") {
                ii = mk::maxtxt + 1;
                txtatm_2 = "(" + txtatm_2.substr(0, len_trim(txtatm_2));
                kk = mk::maxtxt - 1;
              } else {
                ii = mk::maxtxt;
                kk = mk::maxtxt;
              }
              char b[256];
              std::snprintf(b, sizeof b, "%5d%2s%5d%3s%s%s%5d%2s%5.1f%s%5d%s%s%s%s\n",
                            k, ec::elemnt[mc::nat[i]].c_str(), i,
                            padding.substr(0, l - 5).c_str(),
                            ("(" + txtatm_1.substr(0, kk) + ")").c_str(),
                            padding.substr(0, l - 15).c_str(),
                            mz::ions[i], "  (", r_ions[1],
                            mk::line.substr(0, m).c_str(), ec::elemnt[mc::nat[near_ions[1]]].c_str(),
                            near_ions[1], ", Label: ", txtatm_2.substr(0, ii).c_str(), ")");
              std::fprintf(stdout, "%s", b);
              for (ii = 2; ii <= jj; ++ii) {
                if (residues) txtatm_2 = mc::txtatm[near_ions[ii]];
                else txtatm_2 = mc::txtatm1[near_ions[ii]];
                num = "5";
                if (fsub(ec::elemnt[mc::nat[near_ions[ii]]], 1, 1) != " ") num = "5";
                char b2[256];
                std::snprintf(b2, sizeof b2, "%49s%5.1f%s%5d%s%s\n",
                              "   (", r_ions[ii], mk::line.substr(0, m).c_str(),
                              near_ions[ii], ", Label: ",
                              txtatm_2.substr(0, mk::maxtxt).c_str());
                std::fprintf(stdout, "%s", b2);
              }
            } else {
              if (residues) txtatm_1 = mc::txtatm[i];
              else txtatm_1 = mc::txtatm1[i];
              kk = (int)len_trim(txtatm_1);
              if (fsub(txtatm_1, kk, kk) == ")") {
                kk = mk::maxtxt - 1;
                ii = l - 6;
                j2 = l - 1;
              } else {
                kk = mk::maxtxt;
                ii = l - 5;
                j2 = l - 5;
              }
              char b3[256];
              std::snprintf(b3, sizeof b3, "%5d%2s%5d%3s%s%s%5d\n",
                            k, ec::elemnt[mc::nat[i]].c_str(), i,
                            padding.substr(0, ii).c_str(),
                            ("(" + txtatm_1.substr(0, kk) + ")").c_str(),
                            padding.substr(0, j2).c_str(), mz::ions[i]);
              std::fprintf(stdout, "%s", b3);
            }
          } else {
            char b4[128];
            std::snprintf(b4, sizeof b4, "%5d%3s%5d%5s%2s%9d\n",
                          k, "", i, "", ec::elemnt[mc::nat[i]].c_str(), mz::ions[i]);
            std::fprintf(stdout, "%s", b4);
          }
        }
      }
    }
    mk::maxtxt = maxtxt_store;
  }
  if (index1(mk::keywrd, " LEWIS") != 0 && index1(mk::keywrd, " 0SCF") == 0)
    mopend("RUN STOPPED BECAUSE KEYWORD LEWIS WAS USED.");
  if (mz::noccupied != 0 &&
      (charges || !lreseq && (ichrge != 0 || irefq != 0))) {
    num = std::string(1, (char)('3' + std::max((int)std::floor(std::log10(std::abs(ichrge) + 0.05)), 0)));
    std::fprintf(stdout, "\n%17s%5d\n", " COMPUTED CHARGE ON SYSTEM:", ichrge);
  }

  for (i = 1; i <= mk::numat; ++i) mc::atmass[i] = pc::ams[mc::nat[i]];
  if (done && !lreseq) {
    // xyzint expects column-major contiguous xyz/geo buffers; pack them.
    std::vector<double> xyz_c(3 * mk::numat + 3, 0.0), geo_c(3 * mk::numat + 3, 0.0);
    for (i = 1; i <= mk::numat; ++i)
      for (k = 1; k <= 3; ++k) xyz_c[(i - 1) * 3 + (k - 1)] = mc::coord[k-1][i];
    xyzint(xyz_c.data(), mk::numat, mc::na.data(), mc::nb.data(), mc::nc.data(), 1.0,
           geo_c.data());
    for (i = 1; i <= mk::numat; ++i)
      for (k = 1; k <= 3; ++k) mc::geo[k][i] = geo_c[(i - 1) * 3 + (k - 1)];
  }
  if (lreseq) {
    for (i = 1; i <= mk::numat; ++i)
      for (k = 1; k <= 3; ++k) mc::geo[k][i] = mc::coord[k-1][i];
    for (i = 1; i <= mk::maxatoms; ++i) mc::na[i] = 0;
  }
  if (lreseq || lsite) {
    if (ch::archive_fn.size() >= 3) ch::archive_fn = ch::archive_fn.substr(0, ch::archive_fn.size() - 3) + "arc";
    opend = false;  // simplified: inquire(unit=iarc, opened=opend)
    // Remove "PDBOUT" from the reference keywords (F90: refkey(i)=refkey(i)(:j)//refkey(i)(j+7:)).
    for (i = 1; i <= 6; ++i) {
      mk::line = mk::refkey[i];
      upcase(mk::line, (int)len_trim(mk::line));
      j = index1(mk::line, " PDBOUT");
      if (j > 0) mk::refkey[i] = mk::refkey[i].substr(0, j) + mk::refkey[i].substr(j + 6);
    }
    geout(ch::iarc);
    if (index1(mk::keywrd, " PDBOUT") != 0) {
      if (ch::archive_fn.size() >= 3) ch::archive_fn = ch::archive_fn.substr(0, ch::archive_fn.size() - 3) + "pdb";
      for (i = 1; i <= mk::numat; ++i)
        for (k = 1; k <= 3; ++k) mc::coord[k-1][i] = mc::geo[k][i];
      update_txtatm(true, true);
      if (index1(mk::keywrd, " NORES") == 0) rectify_sequence();
      reset_breaks();
      pdbout(ch::iarc);
      mc::coorda.assign(4, std::vector<double>(mk::numat + 1, 0.0));
      for (i = 1; i <= mk::numat; ++i)
        for (k = 1; k <= 3; ++k) mc::coorda[k-1][i] = mc::coord[k-1][i];
    }
  }
  if (irefq != ichrge && !lreseq || charges) {
    // The calculated charge does not match CHARGE=n.
    if (mz::icharges != 0) std::fprintf(stdout, "\n");
    mk::line = " ";
    if (index1(mk::keywrd, " 0SCF") == 0) mk::line = "JOB STOPPED BECAUSE";
    i = (int)len_trim(mk::line);
    if (i > 0) i = i + 1;
    if (charges) {
      if (index1(mk::keywrd, " CHARGES") == 0) {
        mopend(mk::line.substr(0, i) + "CHARGES MODIFIED BY SITE COMMAND");
      } else {
        mopend(mk::line.substr(0, i) + "CHARGES KEYWORD USED");
      }
      mk::nelecs = mk::nelecs - ichrge;
    }
    if (lreseq) mopend(mk::line.substr(0, i) + "GEOMETRY RESEQUENCED");
    if (charges || lreseq) { delete[] ioptl; return; }
    if (index1(mk::keywrd, " LEWIS") != 0) {
      if (index1(mk::keywrd, " 0SCF") == 0) {
        if (!mk::moperr) mopend(mk::line.substr(0, i) + "KEYWORD LEWIS USED");
        delete[] ioptl;
        return;
      }
    }
    if (index1(mk::keywrd, " 0SCF") + index1(mk::keywrd, " RESEQ") +
            index1(mk::keywrd, " LEWIS") == 0) {
      std::fprintf(stdout, "\n");
      std::fprintf(stdout, "%10s%5d%s\n",
                   "In the data-set supplied, the charge specified (", irefq, ") is incorrect.");
      if (index1(mk::keywrd, " GEO-OK") != 0 && index1(mk::keywrd, " CHARGE=") != 0) {
        std::fprintf(stdout, "%10s\n", " KEYWORD 'GEO-OK' WAS PRESENT, SO THE CHARGE HAS BEEN RESET.");
        std::fprintf(stdout, "%10s\n",
                     " IF THE NEW CHARGE IS INCORRECT, EITHER MODIFY THE STRUCTURE OR "
                     "USE KEYWORD 'SETPI' TO CORRECT THE LEWIS STRUCTURE.");
      }
      mk::nelecs = mk::nelecs + irefq - ichrge;
      mk::nclose = mk::nelecs / 2;
      mk::nopen = mk::nclose;
      mk::nalpha = 0;
      mk::nbeta = 0;
      mk::uhf = false;
      if (irefq != ichrge) {
        if ((irefq - ichrge) % 2 == 0) {
          std::fprintf(stdout, "\n%10s\n",
                       "If the Lewis structure used by MOZYME is incorrect, "
                       "use keywords such as CVB or SETPI to correct it");
          web_message(ch::iw, "setpi.html");
        } else {
          std::fprintf(stdout, "\n%10s\n",
                       "The charge keyword, Lewis structure, or the chemical formula is faulty");
        }
        if (index1(mk::keywrd, " GEO-OK") == 0 && index1(mk::keywrd, " CHARGE=") != 0) {
          mopend("CHARGE SPECIFIED IS INCORRECT. CORRECT THE ERROR BEFORE CONTINUING");
          std::fprintf(stdout, "\n%10s\n", "If that is done, then the correct charge will be used.");
          delete[] ioptl;
          return;
        } else if (mk::id > 0) {
          std::fprintf(stdout, "%10s\n", " Infinite systems must have a zero charge on the unit cell.");
          mopend("Unit cell has a charge. Correct fault and re-submit ");
          delete[] ioptl;
          return;
        } else {
          fix_charges(ichrge);
        }
      }
    }
  }
  if (done) {
    if (mk::prt_coords) std::fprintf(stdout, "\n%10s\n", " GEOMETRY AFTER RE-SEQUENCING");
    update_txtatm(true, true);
    if (mk::prt_coords) geout(ch::iw);
    if (index1(mk::keywrd, "0SCF") + index1(mk::keywrd, " RESEQ") == 0 ||
        index1(mk::keywrd, " PDBOUT") == 0) {
      geout(ch::iarc);
    }
    if (index1(mk::keywrd, " PDBOUT") != 0) {
      if (mk::prt_coords) pdbout(1);
    }
    mopend("GEOMETRY RESEQUENCED");
    goto label_1100;
  }
  if (ibad != 0 && !let) {
    mopend("ERROR");
    goto label_1100;
  }
  if (index1(mk::keywrd, " NEWGEO") != 0) {
    newflg();
  }
  if (times) timer(" END OF GEOCHK");
  if (lreseq || lsite) {
    if (lsite) {
      mopend("Keyword SITE used");
      std::fprintf(stdout, "\n%s\n", " Run stopped because SITE used");
    } else {
      mopend("Keyword RESEQ used");
      std::fprintf(stdout, "\n%s\n", " Run stopped because RESEQ used");
    }
    delete[] ioptl;
    return;
  }

  // ---- Modify ions so that it refers to real atoms only --------------------
  j = 0;
  for (i = 1; i <= mk::natoms; ++i) {
    if (mc::labels[i] != 99) {
      j = j + 1;
      mz::iz[i] = j;
    }
  }
  // Restore charges, if present.
  if (mk::maxtxt != 26) {
    for (i = 1; i <= mk::natoms; ++i) {
      if (atom_charge[i] != ' ') setsub(mc::txtatm[i], 2, 2, std::string(1, atom_charge[i]));
    }
  }
  // Restore nbonds and ibonds in case they are modified within this subroutine.
  for (i = 1; i <= mk::numat; ++i) mc::nbonds[i] = nnbonds[i];
  for (i = 1; i <= mk::numat; ++i)
    for (j = 1; j <= 15; ++j) mc::ibonds[j][i] = iibonds[j][i];

label_1100:
  delete[] ioptl;
  return;
}
