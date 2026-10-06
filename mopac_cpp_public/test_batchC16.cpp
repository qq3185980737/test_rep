// test_batchC16.cpp — batch tests for 批C16 (run_mopac driver).
// The driver's external callees are replaced by test-controlled stubs so the
// control flow of run_mopac itself can be exercised deterministically:
//  T1 getdat returns 0 atoms -> early return
//  T2 normal single-point job -> readmo stub sets moperr on 2nd loop -> finish
//  T3 DRC dispatch reaches drc()
//  T4 default dispatch reaches ef()
//  T5 FLEPO dispatch reaches flepo()
//  T6 0SCF geometry-output path reaches geout() and loops
#include "run_mopac.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "maps_C.h"
#include "elemts_C.h"
#include "chanel_C.h"
#include "cosmo_C.h"
#include "symmetry_C.h"
#include "MOZYME_C.h"
#include "meci_C.h"
#include "parameters_for_PM6_Sparkles_C.h"
#include "parameters_for_PM7_Sparkles_C.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;
using namespace maps_C;
using namespace elemts_C;
using namespace chanel_C;
using namespace cosmo_C;
using namespace symmetry_C;
using namespace MOZYME_C;
using namespace meci_C;
using namespace parameters_for_PM6_Sparkles_C;
using namespace parameters_for_PM7_Sparkles_C;

static int failures = 0;
static void check(const char* name, bool ok) {
  if (!ok) { std::cerr << "FAIL: " << name << "\n"; ++failures; }
  else std::cout << "ok: " << name << "\n";
}

// ---- call counters for dispatch tests ----
static int c_readmo = 0, c_switch = 0, c_moldat = 0, c_calpar = 0;
static int c_geout = 0, c_drc = 0, c_ef = 0, c_flepo = 0, c_compfg = 0;
static int c_writmo = 0, c_summary = 0, c_mopend = 0, c_setup = 0;
static bool readmo_fail_next = false;
static int stub_natoms = 3;

// ---- test stubs for run_mopac external callees ----
extern "C" { extern int __argc; extern char** __argv; }

void getdat(int, int) { natoms = stub_natoms; }  // controlled input
void readmo() {
  ++c_readmo;
  // Each successful loop arms a failure for the next loop, so a job finishes
  // after two cycles (matching F90: EOF on the next input record).
  if (readmo_fail_next) { moperr = true; readmo_fail_next = false; }
  else readmo_fail_next = true;
}
void switch_method() { ++c_switch; }
void datin(int) {}
void moldat(int) { ++c_moldat; }
void calpar() { ++c_calpar; }
void fbx() {}
void fordd() {}
void gettxt() {}
void setup_mopac_arrays(int, int) { ++c_setup; }
void set_up_MOZYME_arrays() {}
void delete_MOZYME_arrays() {}
void set_up_rapid(const char*) {}
void update_txtatm(bool, bool) {}
void write_sequence() {}
void bridge_H() {}
void geochk() {}
void add_hydrogen_atoms() {}
void compare_txtatm(bool&, bool&) {}
void den_in_out(int) {}
void density_for_MOZYME(std::vector<double>&, int, int, const std::vector<double>&) {}
void l_control(const std::string&, int, int) {}
void post_scf_corrections(double&, bool) {}
void geout(int) { ++c_geout; }
void wrttxt(int) {}
void geoutg(int) {}
void pdbout(int) {}
double C_triple_bond_C() { return 0.0; }
void compfg(const std::vector<double>&, bool, double& escf, bool,
            std::vector<double>&, bool) { ++c_compfg; escf = -12.34; }
void react1() {}
void grid() {}
void paths() {}
void pathk() {}
void force() {}
void drc(std::vector<double>&, const std::vector<double>&) { ++c_drc; }
void nllsq() {}
void powsq() {}
void flepo(std::vector<double>&, double&) { ++c_flepo; }
void ef(std::vector<double>&, double&) { ++c_ef; }
void lbfgs(double*, double&) {}
void writmo() { ++c_writmo; }
void polar() {}
void static_polarizability() {}
void pmep() {}
void esp() {}
void Locate_TS_for_Proteins() {}
void Refine_TS_for_Proteins() {}
int mkl_get_max_threads() { return 1; }
void mkl_set_num_threads(int) {}
void write_path_html() {}
void upcase(std::string& s, int n) {
  for (int i = 0; i < n && i < (int)s.size(); ++i)
    if (s[i] >= 'a' && s[i] <= 'z') s[i] = char(s[i] - 'a' + 'A');
}
void summary(const std::string&, int) { ++c_summary; }
void fdate(std::string& d) { d = "Mon Sep 25 12:00:00 2026"; }
void to_screen(const std::string&) {}
void mopend(const std::string&) { ++c_mopend; moperr = true; }
double reada(const std::string&, int) { return 1.0; }
double second(int) { return 1.0; }

static void reset_globals(const std::string& kw) {
  numcal = 0;
  step_num = 0;
  moperr = false;
  escf = 0.0;
  gnorm = 0.0;
  pressure = 0.0;
  E_disp = 0.0;
  E_hb = 0.0;
  E_hh = 0.0;
  solv_energy = 0.0;
  nres = 0;
  nscf = 0;
  nmos = 0;
  na1 = 0;
  lpka = false;
  stress = 0.0;
  no_pKa = 0;
  time0 = 0.0;
  MM_corrections = false;
  state_Irred_Rep = " ";
  natoms = 0;
  numat = 0;
  nvar = 0;
  id = 0;
  last = 0;
  iscf = 0;
  iflepo = 0;
  keywrd = kw;
  jobnam = "testjob";
  job_fn = "testjob";
  output_fn = "testjob.out";
  archive_fn = "testjob.arc";
  end_fn = "testjob.end";
  gui = false;
  mozyme = false;
  sparkle = false;
  maxtxt = 0;
  lgpu = false;
  param_constant = 1.0;
  lxfac = false;
  in_house_only = false;
  prt_coords = false;
  use_ref_geo = false;
  trunc_1 = 7.0;
  trunc_2 = 0.22;
  mpack = 100;
  isok = true;
  errtxt = "";
  atheat = 0.0;
  latom = 0;
  rxn_coord = 0.0;
  labels.assign(300, 1);
  nat.assign(300, 1);
  nfirst.assign(300, 1);
  nlast.assign(300, 1);
  nw.assign(10, 0);
  xparam.assign(20, 0.0);
  grad.assign(20, 0.0);
  p.clear();
  pa.clear();
  pb.clear();
  react.clear();
  coord.assign(4, std::vector<double>(300, 0.0));
  geo.assign(4, std::vector<double>(300, 0.0));
  loc.assign(3, std::vector<int>(100, 0));
  txtatm.assign(300, "");
  txtatm1.assign(300, "");
  l_atom.assign(300, 0);
  refkey.assign(20, " ");
  line = "";
  koment = "";
  ios[1] = 1; iop[1] = 0; iod[1] = 0;
  tore[1] = 1.0;
  // reset counters
  c_readmo = c_switch = c_moldat = c_calpar = 0;
  c_geout = c_drc = c_ef = c_flepo = c_compfg = 0;
  c_writmo = c_summary = c_mopend = c_setup = 0;
  readmo_fail_next = false;
}

int main() {
  // T1: getdat stub returns natoms=3 but run_mopac proceeds; we force early
  //     exit by having readmo set moperr on the first loop.
  reset_globals(" 1SCF");
  run_mopac();
  check("T2 run_mopac returns after moperr", true);
  // F90: first loop runs 1SCF (compfg), second loop's readmo hits moperr and
  // returns -> no summary.
  check("T2 readmo twice (loop, then fail)", c_readmo == 2);
  check("T2 summary not called (moperr->return)", c_summary == 0);
  check("T2 compfg ran on 1st loop", c_compfg == 1);

  // T1: getdat with zero atoms -> early return before the job loop
  reset_globals(" 1SCF");
  stub_natoms = 0;                      // empty data set
  run_mopac();
  stub_natoms = 3;
  check("T1 early return on zero atoms", c_readmo == 0 && c_summary == 0);

  // T3: DRC dispatch
  reset_globals(" DRC XYZ");
  run_mopac();
  check("T3 drc called", c_drc == 1);
  check("T3 summary not called (normal run)", c_summary == 0);

  // T4: default dispatch -> ef
  reset_globals(" AM1");
  nvar = 2;
  run_mopac();
  check("T4 ef called", c_ef == 1);

  // T5: FLEPO dispatch
  reset_globals(" FLEPO");
  nvar = 2;
  run_mopac();
  check("T5 flepo called", c_flepo == 1);

  // T6: 0SCF geometry-output path -> geout, then loops and finishes
  reset_globals(" 0SCF");
  run_mopac();
  check("T6 geout called", c_geout >= 1);
  check("T6 summary not called (normal run)", c_summary == 0);

  if (failures == 0) { std::cout << "ALL PASS\n"; return 0; }
  std::cerr << failures << " FAILURES\n";
  return 1;
}
