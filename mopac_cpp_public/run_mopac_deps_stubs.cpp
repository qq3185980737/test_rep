// run_mopac_deps_stubs.cpp — legacy linkage data + platform stand-ins for
// run_mopac and old-style ports.
//
// The routines declared in run_mopac_deps_stubs.h all resolve to real
// translations (2026-09-28):
//   update_txtatm / write_sequence   -> geochk_txtatm.cpp
//   delete_MOZYME_arrays             -> set_up_MOZYME_arrays.cpp
//   write_path_html / add_path       -> pathk.cpp / pdbout.cpp
//   Locate/Refine_TS_for_Proteins    -> big_swap.cpp
//   check_h / check_CVS              -> lewis.cpp
//   set_up_rapid / set_up_dentate    -> set_up_RAPID.cpp / set_up_dentate.cpp
//   hcore_for_MOZYME / picopt / pinout / symtrz -> their own translation files
//   setcup / setup_nhco              -> moldat_helpers.cpp
// CPU_0 / CPU_1 / wall_clock_0 / wall_clock_1 are defined by timer.cpp.
//
// This file keeps only:
//   * MKL thread-control stand-ins (the real MOPAC links against MKL),
//   * legacy bare globals used by old-style ports (symtrz/symtry/tidy/thermo/
//     timout), which coexist with the namespaced module globals.
#include "run_mopac_deps_stubs.h"

#include <string>
#include <vector>

// Real MOPAC uses MKL; minimal stand-in returns 1 thread.
int mkl_get_max_threads() { return 1; }
void mkl_set_num_threads(int /*n*/) {}

int msdel = 0;

// --- legacy bare globals used by symtrz.cpp / symtry.cpp / tidy.cpp /
//     thermo.cpp / timout.cpp (old-style ports) ---
int numat = 0, norbs = 0, maxci = 0, lab = 0, ndep = 0, nclass = 0, iw = 0;
bool moperr = false;
std::vector<int> nfirst, nlast, jndex;
std::vector<double> coord[4];
int nat[1024] = {0}, atmass[1024] = {0};
int idepfn[1024] = {0}, locpar[1024] = {0}, locdep[1024] = {0}, na[1024] = {0};
double depmul[1024] = {0}, geo[1024] = {0};
double elem[4][4][21] = {{{0}}};
std::vector<std::vector<int>> jelem;
std::vector<std::string> namo;
const char* keywrd = "";
