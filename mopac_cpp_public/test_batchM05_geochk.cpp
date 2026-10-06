// test_batchM05_geochk.cpp — batch acceptance tests for the geochk.F90 group
// (M05d): geochk.cpp (main driver), geochk_small.cpp (fix_charges, add_sp*),
// geochk_txtatm.cpp (update_txtatm, rectify_sequence, compare_sequence,
// write_sequence), site.cpp, find_salt_bridges.cpp.
//
// All external Fortran routines referenced by the group are provided as
// minimal stubs below so the group can be linked and exercised standalone.
// Module globals (common_arrays_C, molkst_C, MOZYME_C, ...) are linked from
// their real *_C.cpp translation units.
//
// Exit code 0 = ALL PASS.  Any failed CHECK prints to stderr and exits 1.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "add_hydrogen_atoms.h"
#include "chanel_C.h"
#include "chklew.h"
#include "chkion.h"
#include "common_arrays_C.h"
#include "distance.h"
#include "elemts_C.h"
#include "extvdw_for_MOZYME.h"
#include "find_salt_bridges.h"
#include "findn1.h"
#include "geochk.h"
#include "geout.h"
#include "lewis.h"
#include "ligand.h"
#include "lyse.h"
#include "mod_atomradii.h"
#include "molkst_C.h"
#include "mopend.h"
#include "MOZYME_C.h"
#include "names.h"
#include "newflg.h"
#include "parameters_C.h"
#include "pdbout.h"
#include "reada.h"
#include "reseq.h"
#include "set_up_dentate.h"
#include "timer.h"
#include "upcase.h"
#include "web_message.h"
#include "xyzint.h"

namespace mc = common_arrays_C;
namespace mk = molkst_C;
namespace mz = MOZYME_C;
namespace ch = chanel_C;
namespace pc = parameters_C;
namespace el = elemts_C;
namespace mr = mod_atomradii;

static int failures = 0;
#define CHECK(cond, msg)                                                \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg); \
      ++failures;                                                       \
    }                                                                   \
  } while (0)

// ---------------------------------------------------------------------------
// Call counters
// ---------------------------------------------------------------------------
static int c_lewis = 0, c_chklew = 0, c_chkion = 0, c_findn1 = 0, c_lyse = 0,
    c_names = 0, c_ligand = 0, c_reseq = 0, c_geout = 0, c_pdbout = 0,
    c_extvdw = 0, c_web = 0, c_newflg = 0, c_setup = 0, c_reset = 0,
    c_addH = 0, c_xyzint = 0, c_mopend = 0, c_timer = 0, c_moiety = 0,
    c_cvs = 0, c_checkh = 0;
static bool lewis_done = false;

// ---------------------------------------------------------------------------
// 1-based Fortran helpers used by the stubs
// ---------------------------------------------------------------------------
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
static std::string fmt_pdb_atom(int serial, const std::string& name,
                                const std::string& res3, char chain,
                                int resnum) {
  std::string t(26, ' ');
  setsub(t, 1, 6, "ATOM  ");
  char buf[16];
  std::snprintf(buf, sizeof buf, "%5d", serial);
  setsub(t, 7, 11, buf);
  setsub(t, 13, 16, name);      // 4 chars
  setsub(t, 18, 20, res3);      // residue name
  setsub(t, 22, 22, std::string(1, chain));
  std::snprintf(buf, sizeof buf, "%4d", resnum);
  setsub(t, 23, 26, buf);
  return t;
}

// ---------------------------------------------------------------------------
// Stubs for external routines
// ---------------------------------------------------------------------------
void lewis(bool) {
  ++c_lewis;
  if (!lewis_done) {
    lewis_done = true;
    // H2O-like Lewis structure: 2 sigma bonds (O-H1, O-H2) + 2 lone pairs.
    mz::Lewis_tot = 4;
    mz::Lewis_elem.assign(3, std::vector<int>(5, 0));
    mz::Lewis_elem[1][1] = 1; mz::Lewis_elem[2][1] = 2;
    mz::Lewis_elem[1][2] = 1; mz::Lewis_elem[2][2] = 3;
    mz::Lewis_elem[1][3] = 1; mz::Lewis_elem[2][3] = 0;
    mz::Lewis_elem[1][4] = 1; mz::Lewis_elem[2][4] = 0;
  }
}
void chklew(std::vector<int>& mb, std::vector<int>& numbon, int& l, int, bool) {
  ++c_chklew;
  for (std::size_t i = 1; i < mb.size(); ++i) mb[i] = static_cast<int>(i);
  numbon[1] = 2; numbon[2] = 2; numbon[3] = 0;
  l = 1;
}
void chkion(std::vector<int>&, int&, const std::vector<char>&) { ++c_chkion; }
void findn1(int& n1, const bool* ioptl, int& io) {
  ++c_findn1;
  io = 0;
  n1 = 0;
  for (int i = 1; i <= mk::numat; ++i)
    if (!ioptl[i]) { n1 = i; break; }
}
void reseq(bool* iopt, int* lused, int n1, int& new_, int& io) {
  ++c_reseq;
  io = 0;
  for (int i = 1; i <= mk::numat; ++i)
    if (!iopt[i]) { iopt[i] = true; ++new_; lused[new_] = i; break; }
  (void)n1;
}
void lyse() { ++c_lyse; }
void names(bool* ioptl, int* lused, int n1, int& ires, int, int io,
           int& uni_res, int& mres) {
  ++c_names;
  for (int i = 1; i <= n1 && i <= mk::numat; ++i) {
    if (!ioptl[i]) {
      ioptl[i] = true;
      mc::txtatm[i] = fmt_pdb_atom(i, " N  ", "UNK", 'A', 11 + i);
    }
  }
  if (n1 >= 1 && n1 <= mk::numat) {
    ++uni_res;
    mres = n1;
  }
  ires = (n1 > 0) ? n1 : 0;
  (void)lused; (void)io;
}
void ligand(int&, const std::vector<int>&, int&) { ++c_ligand; }
void moiety(std::vector<bool>&, std::vector<int>&, int, int&) { ++c_moiety; }
void geout(int) { ++c_geout; }
void pdbout(int) { ++c_pdbout; }
void extvdw_for_MOZYME(std::vector<double>& radius, const std::vector<double>&) {
  ++c_extvdw;
  for (std::size_t i = 0; i < radius.size(); ++i) radius[i] = 1.7;
}
void web_message(int, const char*) { ++c_web; }
void timer(const std::string&) { ++c_timer; }
void newflg() { ++c_newflg; }
void set_up_dentate() { ++c_setup; }
void check_CVS(bool) { ++c_cvs; }
void check_h(int& ibad) { ++c_checkh; ibad = 0; }
void reset_breaks() { ++c_reset; }
void add_a_sp3_hydrogen_atom_ext(int, int, int, int, double,
                                 const std::vector<int>&, int) { ++c_addH; }
void mopend(const std::string& msg) {
  ++c_mopend;
  std::fprintf(stderr, "[mopend] %s\n", msg.c_str());
  mk::moperr = true;
}
double reada(const std::string& s, int istart) {
  if (istart < 1) istart = 1;
  const char* p = s.c_str() + (istart - 1);
  while (*p && !(std::isdigit((unsigned char)*p) || *p == '.' ||
                 *p == '+' || *p == '-'))
    ++p;
  if (!*p) return 0.0;
  char* end = nullptr;
  double v = std::strtod(p, &end);
  if (end == p) return 0.0;
  return v;
}
double distance(int a, int b) {
  double dx = mc::coord[1][a] - mc::coord[1][b];
  double dy = mc::coord[2][a] - mc::coord[2][b];
  double dz = mc::coord[3][a] - mc::coord[3][b];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}
void upcase(std::string& s, int n) {
  for (int i = 0; i < n && i < (int)s.size(); ++i)
    s[i] = (char)std::toupper((unsigned char)s[i]);
}
void xyzint(double*, int, int*, int*, int*, double, double*) { ++c_xyzint; }

// ---------------------------------------------------------------------------
// Test harness
// ---------------------------------------------------------------------------
static void reset_globals(int maxatoms) {
  mk::maxatoms = maxatoms;
  mk::numat = 0; mk::natoms = 0; mk::numat_old = 0; mk::numcal = 0;
  mk::nelecs = 0; mk::nalpha = 0; mk::nbeta = 0; mk::nclose = 0;
  mk::nopen = 0; mk::uhf = false; mk::id = 0; mk::maxtxt = 0;
  mk::moperr = false; mk::nvar = 0; mk::prt_coords = false;
  mk::prt_topo = false; mk::ncomments = 0;
  mk::keywrd = " "; mk::line = " "; mk::allkey = " ";
  mk::refkey.assign(7, "");

  mc::nat.assign(maxatoms + 1, 0);
  mc::labels.assign(maxatoms + 1, 0);
  mc::nfirst.assign(maxatoms + 1, 0);
  mc::nlast.assign(maxatoms + 1, 0);
  mc::na.assign(maxatoms + 1, 0);
  mc::nb.assign(maxatoms + 1, 0);
  mc::nc.assign(maxatoms + 1, 0);
  mc::nbonds.assign(maxatoms + 1, 0);
  mc::atmass.assign(maxatoms + 1, 0.0);
  mc::l_atom.assign(maxatoms + 1, ' ');
  mc::breaks.assign(401, 0);
  mc::break_coords.assign(4, std::vector<double>(401, 0.0));
  mc::txtatm.assign(maxatoms + 1, std::string(26, ' '));
  mc::txtatm1.assign(maxatoms + 1, std::string(26, ' '));
  mc::all_comments.assign(4, "");
  mc::chains = std::string(200, 'A');
  mc::geo.assign(4, std::vector<double>(maxatoms + 1, 0.0));
  mc::coord.assign(4, std::vector<double>(maxatoms + 1, 0.0));
  mc::coorda.assign(4, std::vector<double>(maxatoms + 1, 0.0));
  mc::ibonds.assign(16, std::vector<int>(maxatoms + 1, 0));
  mc::lopt.assign(4, std::vector<int>(maxatoms + 1, 0));
  mc::loc.assign(2, std::vector<int>(3 * maxatoms + 1, 0));

  mz::ions.assign(maxatoms + 1, 0);
  mz::iz.assign(maxatoms + 1, 0);
  mz::ib.assign(maxatoms + 1, 0);
  mz::at_res.assign(maxatoms + 1, 0);
  mz::allres.assign(maxatoms + 1, std::string(4, ' '));
  mz::allr.assign(maxatoms + 1, std::string(1, ' '));
  mz::angles.assign(maxatoms + 1, std::vector<double>(4, 0.0));
  mz::bbone.assign(4, std::vector<int>(maxatoms + 1, 0));
  mz::icharges = 0;
  mz::start_res.assign(101, -200);
  mz::res_start.assign(maxatoms + 1, 0);
  mz::Lewis_tot = 0;
  mz::Lewis_elem.assign(maxatoms + 1, std::vector<int>(5, 0));
  mz::noccupied = 0; mz::nvirtual = 0; mz::nres = 0; mz::maxres = 1000;
  mz::tyres.assign(21, "");
  mz::tyr.assign(21, "");
  const char* R3[21] = {"", "GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR",
                        "ASP", "ASN", "LYS", "GLU", "GLN", "ARG", "HIS",
                        "PHE", "TYR", "TRP", "CYS", "MET", "PRO"};
  const char* R1[21] = {"", "G", "A", "V", "L", "I", "S", "T",
                        "D", "N", "K", "E", "Q", "R", "H",
                        "F", "Y", "W", "C", "M", "P"};
  for (int i = 1; i <= 20; ++i) { mz::tyres[i] = R3[i]; mz::tyr[i] = R1[i]; }

  ch::iw = 6; ch::archive_fn = "test.arc"; ch::iarc = 1;

  mr::atom_radius_covalent.assign(108, 1.5);
  el::elemnt.assign(maxatoms + 1, "   ");

  // parameters_C: fill the missing natorb table and main_group flags.
  pc::natorb[1] = 1; pc::natorb[6] = 4; pc::natorb[7] = 4; pc::natorb[8] = 4;
  pc::natorb[15] = 5; pc::natorb[16] = 6; pc::natorb[17] = 7;
  pc::main_group[1] = true; pc::main_group[6] = true; pc::main_group[7] = true;
  pc::main_group[8] = true; pc::main_group[15] = true; pc::main_group[16] = true;
  pc::main_group[17] = true;
}

static void test_fix_charges() {
  std::string rk = " AM1 CHARGE=2";
  rk.resize(241, ' ');
  mk::refkey[1] = rk;
  mk::keywrd = rk;
  fix_charges(0);
  CHECK(mk::refkey[1].find("CHARGE") == std::string::npos,
        "fix_charges(0) should delete CHARGE from refkey(1)");

  rk = " AM1 CHARGE=2";
  rk.resize(241, ' ');
  mk::refkey[1] = rk;
  mk::keywrd = rk;
  fix_charges(-1);
  CHECK(mk::refkey[1].find(" CHARGE=-1") != std::string::npos,
        "fix_charges(-1) should write CHARGE= -1 into refkey(1)");
  CHECK(mk::keywrd.find(" CHARGE=-1") != std::string::npos,
        "fix_charges(-1) should write CHARGE= -1 into keywrd");
}

static void test_add_sp() {
  mk::maxatoms = 10; mk::numat = 3; mk::natoms = 3;
  mc::nat[1] = 8; mc::nat[2] = 1; mc::nat[3] = 1;
  for (int k = 1; k <= 3; ++k) {
    mc::geo[k][1] = 0.0; mc::geo[k][2] = 0.0; mc::geo[k][3] = 0.0;
  }
  mc::geo[1][2] = 1.0; mc::geo[2][3] = 1.0;

  add_sp_H(1, 2, 3);
  CHECK(mk::natoms == 4, "add_sp_H should add one atom");
  CHECK(mc::nat[4] == 1, "add_sp_H atom should be hydrogen");
  {
    double dx = mc::geo[1][4] - mc::geo[1][1];
    double dy = mc::geo[2][4] - mc::geo[2][1];
    double dz = mc::geo[3][4] - mc::geo[3][1];
    double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    CHECK(std::fabs(d - 1.1) < 1e-9,
          "add_sp_H geometry should place H 1.1 A from i1");
  }

  add_sp2_H(1, 2, 3);
  CHECK(mk::natoms == 5, "add_sp2_H should add one atom");
  {
    double dx = mc::geo[1][5] - mc::geo[1][2];
    double dy = mc::geo[2][5] - mc::geo[2][2];
    double dz = mc::geo[3][5] - mc::geo[3][2];
    double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    CHECK(std::fabs(d - 1.1) < 1e-9,
          "add_sp2_H geometry should place H 1.1 A from i");
  }

  add_sp3_H(1, 2, 3, 4);
  CHECK(mk::natoms == 6, "add_sp3_H should add one atom");
  {
    double dx = mc::geo[1][6] - mc::geo[1][2];
    double dy = mc::geo[2][6] - mc::geo[2][2];
    double dz = mc::geo[3][6] - mc::geo[3][2];
    double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    CHECK(std::fabs(d - 1.1) < 1e-9,
          "add_sp3_H geometry should place H 1.1 A from i");
  }
}

static void test_rectify_sequence() {
  // Ordered input: no chain break -> sequence must be preserved.
  reset_globals(50);
  mk::maxtxt = 26;
  mk::numat = 3; mk::natoms = 3;
  mc::nat[1] = 7; mc::nat[2] = 7; mc::nat[3] = 7;
  mc::labels[1] = 7; mc::labels[2] = 7; mc::labels[3] = 7;
  mc::txtatm1[1] = fmt_pdb_atom(1, " N  ", "GLY", 'A', 12);
  mc::txtatm1[2] = fmt_pdb_atom(2, " N  ", "GLY", 'A', 13);
  mc::txtatm1[3] = fmt_pdb_atom(3, " N  ", "GLY", 'A', 14);
  for (int i = 1; i <= 3; ++i) mc::txtatm[i] = mc::txtatm1[i];
  for (int i = 1; i <= 3; ++i) {
    mc::coord[1][i] = i * 1.0; mc::coord[2][i] = 0.0; mc::coord[3][i] = 0.0;
  }
  rectify_sequence();
  int r1 = (int)(reada(mc::txtatm[1], 23) + 0.5);
  int r2 = (int)(reada(mc::txtatm[2], 23) + 0.5);
  int r3 = (int)(reada(mc::txtatm[3], 23) + 0.5);
  CHECK(r1 == 12 && r2 == 13 && r3 == 14,
        "rectify_sequence: ordered input must be preserved");

  // Out-of-order input: fragment 2 (residues 15,16) followed by the
  // out-of-sequence residue 13, then fragment 3 (residue 17).  The Fortran
  // algorithm (cumulative slot counter, fragment-1 chain head excluded) moves
  // residue 17 up into slot 4 and pushes residue 13 into slot 5.
  reset_globals(50);
  mk::maxtxt = 26;
  mk::numat = 5; mk::natoms = 5;
  mc::nat.assign(51, 7);
  mc::labels.assign(51, 7);
  mc::txtatm1[1] = fmt_pdb_atom(1, " N  ", "GLY", 'A', 12);
  mc::txtatm1[2] = fmt_pdb_atom(2, " N  ", "GLY", 'A', 15);
  mc::txtatm1[3] = fmt_pdb_atom(3, " N  ", "GLY", 'A', 16);
  mc::txtatm1[4] = fmt_pdb_atom(4, " N  ", "GLY", 'A', 13);
  mc::txtatm1[5] = fmt_pdb_atom(5, " N  ", "GLY", 'A', 17);
  for (int i = 1; i <= 5; ++i) mc::txtatm[i] = mc::txtatm1[i];
  for (int i = 1; i <= 5; ++i) {
    mc::coord[1][i] = i * 1.0; mc::coord[2][i] = 0.0; mc::coord[3][i] = 0.0;
  }
  rectify_sequence();
  r1 = (int)(reada(mc::txtatm[1], 23) + 0.5);
  r2 = (int)(reada(mc::txtatm[2], 23) + 0.5);
  int r3b = (int)(reada(mc::txtatm[3], 23) + 0.5);
  int r4 = (int)(reada(mc::txtatm[4], 23) + 0.5);
  int r5 = (int)(reada(mc::txtatm[5], 23) + 0.5);
  CHECK(r1 == 12, "rectify_sequence: slot 1 keeps residue 12");
  CHECK(r2 == 15, "rectify_sequence: slot 2 keeps residue 15");
  CHECK(r3b == 16, "rectify_sequence: slot 3 keeps residue 16");
  CHECK(r4 == 17, "rectify_sequence: residue 17 moved up into slot 4");
  CHECK(r5 == 13, "rectify_sequence: out-of-order residue 13 moved to slot 5");
}

static void test_compare_sequence() {
  reset_globals(50);
  mk::maxtxt = 26;
  mk::numat = 2; mk::natoms = 2;
  mc::nat[1] = 7; mc::nat[2] = 6;
  mc::txtatm[1] = fmt_pdb_atom(1, " N  ", "UNK", 'A', 12);
  mc::txtatm[2] = fmt_pdb_atom(2, " CA ", "UNK", 'A', 12);
  mc::txtatm1[1] = fmt_pdb_atom(1, " N  ", "GLY", 'A', 12);
  mc::txtatm1[2] = fmt_pdb_atom(2, " CA ", "GLY", 'A', 12);
  mc::breaks.assign(401, 0);
  mk::line = " ";
  compare_sequence(0);
  CHECK(mk::line == "(Use the XENO keyword to re-define unrecognized residues.)",
        "compare_sequence should finish with the XENO advice line");
  CHECK(c_mopend == 0, "compare_sequence should not call mopend");
}

static void test_write_sequence() {
  reset_globals(50);
  mk::numat = 3; mk::natoms = 3;
  mc::nat[1] = 7; mc::nat[2] = 6; mc::nat[3] = 6;
  mc::txtatm[1] = fmt_pdb_atom(1, " N  ", "GLY", 'A', 5);
  mc::txtatm[2] = fmt_pdb_atom(2, " CA ", "GLY", 'A', 5);
  mc::txtatm[3] = fmt_pdb_atom(3, " C  ", "GLY", 'A', 5);
  mz::ions[1] = 0; mz::ions[2] = 0; mz::ions[3] = 0;
  write_sequence();
  CHECK(fsub(mz::allres[5], 1, 3) == "GLY",
        "write_sequence should record residue 5 as GLY");
  CHECK(mz::allr[5] == "G", "write_sequence should map GLY to one-letter G");
}

static void test_find_salt_bridges() {
  reset_globals(60);
  mk::numat = 14; mk::natoms = 14;
  // Cationic end (neutral Arg-like guanidinium with 4 H on the three N).
  mc::nat[1] = 6;   // guanidinium C, bonded to N2 N3 N4
  mc::nat[2] = 7;   // N, bonded to C1 H5 H6
  mc::nat[3] = 7;   // N, bonded to C1 H7 H8
  mc::nat[4] = 7;   // N, bonded to C1 C9 (no H)
  mc::nat[5] = 1; mc::nat[6] = 1; mc::nat[7] = 1; mc::nat[8] = 1;
  mc::nat[9] = 6;   // backbone C
  mc::nbonds[1] = 3; mc::nbonds[2] = 3; mc::nbonds[3] = 3;
  mc::nbonds[4] = 2; mc::nbonds[5] = 1; mc::nbonds[6] = 1;
  mc::nbonds[7] = 1; mc::nbonds[8] = 1; mc::nbonds[9] = 1;
  mc::ibonds[1][1] = 2; mc::ibonds[2][1] = 3; mc::ibonds[3][1] = 4;
  mc::ibonds[1][2] = 1; mc::ibonds[2][2] = 5; mc::ibonds[3][2] = 6;
  mc::ibonds[1][3] = 1; mc::ibonds[2][3] = 7; mc::ibonds[3][3] = 8;
  mc::ibonds[1][4] = 1; mc::ibonds[2][4] = 9;
  mc::ibonds[1][5] = 2; mc::ibonds[1][6] = 2; mc::ibonds[1][7] = 3;
  mc::ibonds[1][8] = 3; mc::ibonds[1][9] = 4;
  // Anionic end (protonated Asp-like -COOH: 1 H on the OH oxygen).
  mc::nat[10] = 6;  // carboxyl C, bonded to C11 O12 O13
  mc::nat[11] = 6;  // backbone C
  mc::nat[12] = 8;  // O, bonded to C10
  mc::nat[13] = 8;  // O, bonded to C10 H14
  mc::nat[14] = 1;
  mc::nbonds[10] = 3; mc::nbonds[11] = 1; mc::nbonds[12] = 1;
  mc::nbonds[13] = 2; mc::nbonds[14] = 1;
  mc::ibonds[1][10] = 11; mc::ibonds[2][10] = 12; mc::ibonds[3][10] = 13;
  mc::ibonds[1][11] = 10; mc::ibonds[1][12] = 10;
  mc::ibonds[1][13] = 10; mc::ibonds[2][13] = 14;
  mc::ibonds[1][14] = 13;
  // Labels: Arg residues (chain A, residue 12), Asp residues (chain A, 30).
  for (int i = 1; i <= 9; ++i)
    mc::txtatm[i] = fmt_pdb_atom(i, " N  ", "ARG", 'A', 12);
  for (int i = 10; i <= 14; ++i)
    mc::txtatm[i] = fmt_pdb_atom(i, " C  ", "ASP", 'A', 30);
  // Coordinates: place a guanidinium N 3.0 A from a carboxyl O.
  for (int i = 1; i <= 14; ++i)
    for (int k = 1; k <= 3; ++k) mc::coord[k][i] = 0.0;
  mc::coord[1][2] = 0.0; mc::coord[1][3] = 1.2; mc::coord[1][4] = -1.2;
  mc::coord[1][13] = 3.0; mc::coord[1][12] = 3.6;
  mc::coord[1][10] = 3.3; mc::coord[1][11] = 4.2;
  // n_cat = 0 : auto-detect and build the expanded SITE keyword.
  mk::keywrd = " SITE=(SALT)";
  mk::maxtxt = 26;
  static int empty_cat[1] = {0}, empty_ani[1] = {0};
  int mopend_before = c_mopend;
  find_salt_bridges(empty_cat, empty_ani, 0, 0);
  CHECK(c_mopend == mopend_before,
        "find_salt_bridges should not call mopend when a bridge is found");
  CHECK(mk::keywrd.find("SITE=(") != std::string::npos,
        "find_salt_bridges should rebuild the SITE keyword");
  CHECK(mk::keywrd.find("SALT") == std::string::npos,
        "find_salt_bridges should remove SALT from the SITE keyword");
  CHECK(mk::keywrd.find("(+)") != std::string::npos,
        "find_salt_bridges should tag the cationic site with (+)");
  CHECK(mk::keywrd.find("(-)") != std::string::npos,
        "find_salt_bridges should tag the anionic site with (-)");
}

static void test_site_minimal() {
  reset_globals(50);
  mk::numat = 2; mk::natoms = 2;
  mc::nat[1] = 8; mc::nat[2] = 1;
  mc::nbonds[1] = 2; mc::nbonds[2] = 1;
  mc::ibonds[1][1] = 2; mc::ibonds[2][1] = 0; mc::ibonds[1][2] = 1;
  mc::coord[1][1] = 0.0; mc::coord[2][1] = 0.0; mc::coord[3][1] = 0.0;
  mc::coord[1][2] = 0.96; mc::coord[2][2] = 0.0; mc::coord[3][2] = 0.0;
  mc::txtatm[1] = fmt_pdb_atom(1, " O  ", "HOH", 'A', 1);
  mc::txtatm[2] = fmt_pdb_atom(2, " H  ", "HOH", 'A', 1);
  mk::keywrd = " SITE=(COO) ";
  std::string allkey = " ";
  std::vector<bool> neutral(11, false);
  std::vector<char> chain(2, ' ');
  std::vector<int> res(2, 0);
  std::vector<std::vector<char>> charge(2, std::vector<char>(4, ' '));
  int mopend_before = c_mopend;
  site(neutral, chain, res, charge, 0, 0, allkey);
  CHECK(c_mopend == mopend_before,
        "site with no matching residues should be a no-op");
  CHECK(mk::numat == 2, "site no-op should not change the atom count");
}

static void test_geochk_geo_ok() {
  // Scenario A: H2O, keyword GEO-OK.  lewis() stub provides the Lewis
  // structure; geochk should compute charges, restore connectivity, and stop.
  reset_globals(100);
  mk::keywrd = " GEO-OK";
  mk::numat = 3; mk::natoms = 3; mk::numat_old = 3;
  mk::nelecs = 8; mk::nalpha = 0; mk::nbeta = 0;
  mk::nclose = 4; mk::nopen = 4; mk::uhf = false; mk::id = 0;
  mk::maxtxt = 0;
  mc::nat[1] = 8; mc::nat[2] = 1; mc::nat[3] = 1;
  mc::labels[1] = 8; mc::labels[2] = 1; mc::labels[3] = 1;
  mc::nbonds[1] = 2; mc::nbonds[2] = 1; mc::nbonds[3] = 1;
  mc::ibonds[1][1] = 2; mc::ibonds[2][1] = 3;
  mc::ibonds[1][2] = 1; mc::ibonds[1][3] = 1;
  for (int k = 1; k <= 3; ++k) {
    mc::coord[k][1] = 0.0; mc::geo[k][1] = 0.0;
    mc::coord[k][2] = 0.0; mc::geo[k][2] = 0.0;
    mc::coord[k][3] = 0.0; mc::geo[k][3] = 0.0;
  }
  mc::coord[1][2] = 0.96; mc::geo[1][2] = 0.96;
  mc::coord[1][3] = -0.24; mc::coord[2][3] = 0.93;
  mc::geo[1][3] = -0.24; mc::geo[2][3] = 0.93;
  mc::coorda = mc::coord;

  int lewis_before = c_lewis, chklew_before = c_chklew;
  int mopend_before = c_mopend;
  geochk();
  CHECK(c_lewis > lewis_before, "geochk should call lewis()");
  CHECK(c_chklew > chklew_before, "geochk should call chklew()");
  CHECK(c_mopend == mopend_before,
        "GEO-OK neutral H2O should not hit mopend");
  CHECK(mz::noccupied == 4, "geochk should set noccupied from Lewis data");
  CHECK(mz::ions[1] == 0, "oxygen of neutral H2O should have 0 ions");
  CHECK(mz::ions[2] == 0, "hydrogen 1 of H2O should have 0 ions");
  CHECK(mz::ions[3] == 0, "hydrogen 2 of H2O should have 0 ions");
  CHECK(mk::numat == 3, "geochk should keep numat unchanged");
  CHECK(mc::nbonds[1] == 2 && mc::nbonds[2] == 1 && mc::nbonds[3] == 1,
        "geochk should restore nbonds");
  CHECK(mc::ibonds[1][1] == 2 && mc::ibonds[2][1] == 3 &&
        mc::ibonds[1][2] == 1 && mc::ibonds[1][3] == 1,
        "geochk should restore ibonds");
}

static void test_geochk_resi_xeno() {
  // Scenario B: three-atom chain, " RESI 0SCF" + XENO renames UNK residues.
  // lres path: lyse -> names -> chain letters -> XENO rename -> write_sequence.
  reset_globals(100);
  mk::keywrd = " RESI 0SCF XENO=(A12=GLY,A13=ASP) ";
  mk::numat = 3; mk::natoms = 3; mk::numat_old = 3;
  mk::nelecs = 15; mk::nalpha = 0; mk::nbeta = 0;
  mk::nclose = 7; mk::nopen = 1; mk::uhf = true; mk::id = 0;
  mk::maxtxt = 26;
  mc::nat[1] = 8; mc::nat[2] = 7; mc::nat[3] = 6;
  mc::labels[1] = 8; mc::labels[2] = 7; mc::labels[3] = 6;
  mc::nbonds[1] = 1; mc::nbonds[2] = 1; mc::nbonds[3] = 1;
  mc::ibonds[1][1] = 2; mc::ibonds[1][2] = 1; mc::ibonds[1][3] = 0;
  mc::ibonds[1][2] = 1;
  for (int k = 1; k <= 3; ++k) {
    mc::coord[k][1] = 0.0; mc::geo[k][1] = 0.0;
    mc::coord[k][2] = 0.0; mc::geo[k][2] = 0.0;
    mc::coord[k][3] = 0.0; mc::geo[k][3] = 0.0;
  }
  mc::coord[1][1] = 0.0; mc::coord[1][2] = 1.0; mc::coord[1][3] = 2.0;
  mc::geo = mc::coord;
  mc::coorda = mc::coord;
  mc::breaks.assign(401, 0);

  int lyse_before = c_lyse, names_before = c_names;
  int mopend_before = c_mopend;
  geochk();
  CHECK(c_lyse > lyse_before, "geochk RESI path should call lyse()");
  CHECK(c_names > names_before, "geochk RESI path should call names()");
  CHECK(c_mopend == mopend_before,
        "valid RESI + XENO path should not hit mopend");
  CHECK(fsub(mc::txtatm[1], 18, 20) == "GLY",
        "XENO should rename residue A12 UNK -> GLY");
  CHECK(fsub(mc::txtatm[2], 18, 20) == "ASP",
        "XENO should rename residue A13 UNK -> ASP");
  CHECK(mc::nbonds[1] == 1 && mc::nbonds[2] == 1 && mc::nbonds[3] == 1,
        "RESI path should restore nbonds");
  CHECK(mc::ibonds[1][1] == 2 && mc::ibonds[1][2] == 1,
        "RESI path should restore ibonds");
}

int main() {
  reset_globals(100);

  test_fix_charges();
  test_add_sp();
  test_rectify_sequence();
  test_compare_sequence();
  test_write_sequence();
  test_find_salt_bridges();
  test_site_minimal();
  test_geochk_geo_ok();
  test_geochk_resi_xeno();

  if (failures == 0) {
    std::printf("ALL PASS: M05d geochk group\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) FAILED\n", failures);
  return 1;
}
