// run_mopac.cpp — C++ translation of "run_mopac.F90" (868 lines).
//
// Main driver: read the data set (getdat/readmo), load the method (switch,
// datin, moldat, calpar), dispatch to the requested calculation (react1,
// grid, paths, force, drc, nllsq, powsq, compfg, flepo, ef, lbfgs, ...),
// print results (writmo, polar, static_polarizability, pmep, esp), then loop
// back for the next job (label 10) or finish (label 101).
//
// Notes on faithful adaptation:
//  - F90 open/write/close on unit numbers (iw/ir/iarc/ilog) have no direct
//    C++ mapping; output goes to stdout/stderr, open/close/rewind are kept as
//    no-op placeholders with the F90 intent in comments.
//  - The #if GPU block is not translated (no GPU implementation in this build).
//  - MKL thread calls map to minimal stubs (mkl_get_max_threads/mkl_set_num_threads).

#include "run_mopac.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

#include "add_hydrogen_atoms.h"
#include "big_swap.h"
#include "calpar.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "compfg.h"
#include "cosmo_C.h"
#include "datin.h"
#include "density_for_MOZYME.h"
#include "den_in_out.h"
#include "drc.h"
#include "ef.h"
#include "elemts_C.h"
#include "esp.h"
#include "fdate.h"
#include "flepo.h"
#include "force.h"
#include "funcon_C.h"
#include "geochk.h"
#include "geout.h"
#include "geoutg.h"
#include "getdat_lines_C.h"
#include "getdat.h"
#include "gettxt.h"
#include "grid.h"
#include "lbfgs.h"
#include "maps_C.h"
#include "meci_C.h"
#include "mndod.h"
#include "moldat.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "mopend.h"
#include "nllsq.h"
#include "parameters_C.h"
#include "parameters_for_PM6_Sparkles_C.h"
#include "parameters_for_PM7_Sparkles_C.h"
#include "pathk.h"
#include "paths.h"
#include "pdbout.h"
#include "pmep.h"
#include "polar.h"
#include "post_scf_corrections.h"
#include "powsq.h"
#include "react1.h"
#include "reada.h"
#include "readmo.h"
#include "run_mopac_deps_stubs.h"
#include "second.h"
#include "set_up_dentate.h"
#include "set_up_MOZYME_arrays.h"
#include "set_up_RAPID.h"
#include "setup_mopac_arrays.h"
#include "static_polarizability.h"
#include "switch.h"
#include "symmetry_C.h"
#include "to_screen.h"
#include "upcase.h"
#include "writmo.h"
#include "wrttxt.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace cosmo_C;
using namespace elemts_C;
using namespace funcon_C;
using namespace maps_C;
using namespace meci_C;
using namespace molkst_C;
using namespace MOZYME_C;
using namespace parameters_C;
using namespace parameters_for_PM7_Sparkles_C;
using namespace parameters_for_PM6_Sparkles_C;
using namespace symmetry_C;

extern void geoutg(int iprt);
extern void compare_txtatm(bool& a, bool& b);
extern void summary(const std::string& txt, int ntxt);
// Timers consumed by timout(): defined in timer.cpp, updated here per job.
extern double wall_clock_0, wall_clock_1, CPU_0, CPU_1;

void run_mopac() {
  // ---- one-time initialization (before the job loop) ----
  for (int i = 1; i <= N_PARAM; ++i) tore[i] = ios[i] + iop[i] + iod[i];
  fbx();                // factorials and Pascal's triangle (pure constants)
  fordd();              // more constants for MNDO-d
  param_constant = 1.0;
  lgpu = false;
  if (param_constant < -1.e5) return;
  trunc_1 = 7.0;        // beyond 7.0 A, use exact point-charge
  trunc_2 = 0.22;       // multiplier in Gaussian: exp(-trunc_2*(trunc_1 - Rab)^2)
  //
  // Read in all data; put it into a scratch file "ir"
  //
  getdat(ir, iw);       // ir: temp file; natoms: number of lines incl. blanks/comments
  to_screen("To_file: Start of reading in data");
  if (natoms == 0 || moperr) return;
  //
  // Time budget (wrtkey.F90 equivalent): default two days, overridable by T=
  //
  tleft = 172800.0;
  {
    const std::size_t ti = keywrd.find(" T=");
    if (ti != std::string::npos) {
      tleft = reada(keywrd, static_cast<int>(ti));
      tleft = std::min(1e7 - 1.0, tleft);
    }
  }
  //
  // CLOSE UNIT IW IN CASE IT WAS ALREADY PRE-ASSIGNED
  // (F90 open(unit=iw,file=output_fn,...): retry up to 20 times if unavailable)
  //
  int l = 0;
  {
    // open(unit=iw, file=output_fn, status='UNKNOWN', position='asis', iostat=i)
    // C++ mapping: redirect stdout to output_fn so all write(iw) go to the .out file
    // (official MOPAC behavior: the .out file is the main product).
    FILE* out_fp = std::freopen(output_fn.c_str(), "w", stdout);
    if (out_fp == nullptr) {
      ++l;
      std::fprintf(stderr, "%3d File \"%s\" is unavailable for use\n", 21 - l,
                   output_fn.c_str());
      std::fprintf(stderr, "    Correct fault (probably the output file is in use elsewhere)\n");
      mopend("Job abandoned due to output file being unavailable");
      goto l101;
    }
  }
  if (l > 0) std::fprintf(stderr, " Fault successfully corrected.  Job continuing\n");
  // rewind iw — no-op in C++ mapping (F90: move file pointer to start).
  isok = true;
  errtxt = "Job stopped by operator";
  double tim = second(1);
  {
    std::time_t now = std::time(nullptr);
    std::tm tmv{};
    localtime_s(&tmv, &now);
    time_start[0] = tmv.tm_year + 1900;
    time_start[1] = tmv.tm_mon + 1;
    time_start[2] = tmv.tm_mday;
    time_start[3] = 0;                     // GMT offset (dummy)
    time_start[4] = tmv.tm_hour;
    time_start[5] = tmv.tm_min;
    time_start[6] = tmv.tm_sec;
    time_start[7] = 0;                     // milliseconds
  }
  if (moperr) goto l101;
  //
  // Set up essential arrays needed for reading in the data
  //
  natoms = natoms + 200;                    // 200+ total lines, e.g. hf=200+5
  setup_mopac_arrays(natoms, 1);            // mode=1: basic arrays (coord etc.)
  maxatoms = natoms;
  // Unique definition of in_house_only: check for a file.
  in_house_only = false;
  {
    std::FILE* fp = std::fopen("M:\\Survey_Of_Solids\\bits.txt", "r");
    if (fp) { in_house_only = true; std::fclose(fp); }
  }
  if (!in_house_only) {
    std::FILE* fp = std::fopen("/Users/jstewart/fa3207.cif", "r");
    if (fp) { in_house_only = true; std::fclose(fp); }
  }

  // ===================== main job loop (F90 label 10) =====================
  // Loop-body locals are hoisted to function scope so forward jumps (goto
  // LendBody / L101) do not cross initialized declarations.
  int i = 0, j = 0;
  double eat = 0.0, rmx = 0.0;
  bool exists = false, opend = false, sparkles = false;
  int st = 0;            // 0: normal; 1: goto label 100; 2: goto label 101
  while (true) {
    // ---- label 10 body: initialize a new calculation ----
    numcal = numcal + 1;          // a new calculation -> numcal=1
    step_num = step_num + 1;      // new electronic structure
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
    time0 = second(1);
    MM_corrections = false;
    state_Irred_Rep = " ";

    if (numcal > 1) to_screen("To_file: Leaving MOPAC");
    if (numcal > 1 && numcal < 4 && keywrd.find(" GEO_DAT") != std::string::npos) {
      // Quickly jump over the first three lines
      // close(ir) — no-op in C++ mapping.
      i = natoms;
      getdat(ir, iw);
      natoms = i;
      gettxt();                       // record the first three lines
    }
    //
    // Read in all the data for the current job
    //
    i = numcal;
    readmo();                         // atoms: 3 header lines + geometry
    if (moperr && numcal == 1 && keywrd.find(" GEO_DAT") == std::string::npos)
      goto LendBody;                  // F90: goto 100
    if (moperr) { st = 2; goto LendBody; }  // F90: goto 101 (finish with TOTAL JOB TIME)
    if (numcal == 1) {
      // Reference run_mopac.F90 sets num_threads = mkl_get_max_threads(); the
      // reference exe returns 4 under OMP_NUM_THREADS=4.  Our static-MKL build
      // currently resolves mkl_get_max_threads() to 1 even with the env set
      // (OpenMP runtime init quirk), which would pin everything to 1 thread and
      // change every MKL reduction order.  Read the env directly instead.
      num_threads = 1;
      const char* en = std::getenv("OMP_NUM_THREADS");
      if (en) { int v = std::atoi(en); if (v >= 1) num_threads = v; }
      const char* em = std::getenv("MKL_NUM_THREADS");
      if (em) { int v = std::atoi(em); if (v >= 1) num_threads = v; }
      if (in_house_only) num_threads = 1;
      i = static_cast<int>(keywrd.find(" THREADS"));
      if (i > 0) {
        i = static_cast<int>(std::lround(reada(keywrd, i)));
        num_threads = std::min(std::max(1, i), num_threads);
      }
      mkl_set_num_threads(num_threads);
      // (F90 #if GPU block not translated.)
    }
    if (!gui && numcal == 1 && natoms == 0) {
      line = " Data set exists, but does not contain any atoms.";
      std::fprintf(stderr, "\n\n%20s\n\n", line.c_str());
      mopend(line);
      std::fprintf(stderr, " (Check the first few lines of the data-set for an extra blank line.\n");
      std::fprintf(stderr, "  If there is an extra line, delete it and re-submit.)\n");
      // inquire(unit=iw, opened=opend); if (opend) close(ir,status='delete')
      return;
    }
    if (numat > 46000) {
      line = " Data set '" + jobnam + "' exists, but is much too large to run.";
      std::fprintf(stderr, "\n\n%20s\n\n", line.c_str());
      std::printf("\n\n%20s\n\n", line.c_str());
      mopend(line);
      numat = 0;
      natoms = 0;
      goto LendBody;                  // F90: goto 100
    }
    if ((numcal == 1 && moperr) || natoms == 0) {
      // Check for spurious "extra" data
      // (F90: read up to 15 lines from unit ir; kept as a comment stub)
      // i > 14 -> WARNING: extra data at end of input
      st = 2;
      goto LendBody;                  // F90: goto 101
    }
    if (moperr) {
      if (keywrd.rfind(" MODEL", 0) == 0) { st = 2; goto LendBody; }  // goto 101
      // goto 10 — loop again
      continue;
    }
    lxfac = (keywrd.find(" XFAC") != std::string::npos);
    //
    // Load in parameters for the method to be used
    //
    switch_method();                  // F90 call switch
    if (keywrd.find(" EXTERNAL") != std::string::npos) datin(iw);
    if (moperr) goto LendBody;        // F90: goto 100
    sparkle = (keywrd.find(" SPARKL") != std::string::npos);
    sparkles = false;
    for (i = 1; i <= natoms; ++i) {
      if (labels[i] > 57 && labels[i] < 71) {
        sparkles = true;
        break;
      }
    }
    if ((method_pm7 || method_PM6 || method_rm1) && sparkles) {
      if (!sparkle) {
        line = " ";
        for (j = 0; j < 10; ++j) {
          if (atom_names[labels[i]][j] != ' ') break;
        }
        line = " Data are not available for " + atom_names[labels[i]].substr(j) + ".";
        std::printf("%s\n", line.c_str());
        if (keywrd.find(" 0SCF") == std::string::npos) {
          if ((method_pm7 && std::fabs(gss7sp[labels[i]]) > 0.1) ||
              (method_PM6 && std::fabs(parameters_for_PM6_Sparkles_C::gss6sp[labels[i]]) > 0.1))
            std::printf(" (Parameters are available if SPARKLE is used)\n");
          mopend(line);
          goto LendBody;              // F90: goto 100
        }
      }
    } else {
      sparkle = false;
    }
    for (i = 57; i <= 71; ++i) {
      if (zs[i] < 0.1) tore[i] = 3.0;
    }
    //
    // Set up all the data for the molecule
    //
    moldat(0);                        // orbitals, spin, two-electron count
    calpar();                         // derived parameters
    if (moperr) goto LendBody;        // F90: goto 100
    to_screen("To_file: Data read in");
    //
    // If no SCF calculations are needed, output geometry and quit
    // (F90 dead-code block "if (.false.) then" is omitted; see source lines
    //  362-394 for the archive/geout variant that is unreachable.)
    //
    if (keywrd.find(" 0SCF") != std::string::npos ||
        keywrd.find(" RESEQ") != std::string::npos) {
      // open(unit=iarc, file=archive_fn...); rewind iarc — no-op mapping.
      if (keywrd.find(" PDBOUT") != std::string::npos) {
        line = archive_fn.substr(0, archive_fn.size() - 3) + "pdb";
      } else {
        line = archive_fn;
      }
      if (keywrd.find(" HTML") != std::string::npos) {
        for (i = static_cast<int>(line.size()); i >= 1; --i) {
          if (line[i - 1] == '/' || line[i - 1] == '\\') break;
        }
      }
      if (keywrd.find(" 0SCF") != std::string::npos &&
          keywrd.find(" MINI") != std::string::npos) {
        // open(unit=l, file=xyz_fn); write header; coordinates
        std::printf("%6d \n", nl_atoms);
        std::printf("POINT \n");
        for (i = 1; i <= numat; ++i) {
          if (l_atom[i])
            std::printf("%5s%15.5f%15.5f%15.5f\n", elemnt[nat[i]].c_str(),
                        coord[0][i], coord[1][i], coord[2][i]);
        }
      }
    }
    if (keywrd.find(" PDBOUT") != std::string::npos && maxtxt != 26 &&
        keywrd.find(" RESID") == std::string::npos) {
      if (maxtxt == 0) {
        maxtxt = 26;
        for (i = 1; i <= natoms - id; ++i) {
          char buf[64];
          std::snprintf(buf, sizeof(buf), "HETATM%5d %s   HET A   1", i,
                        elemnt[nat[i]].c_str());
          txtatm[i] = buf;
        }
        for (i = 1; i <= natoms - id; ++i) txtatm1[i] = txtatm[i];
      } else {
        line = "PDBOUT only works when the atom labels are in PDB format or keyword RESIDUES is also present";
        mopend(line);
        std::printf("\n%18s\n", "(Before using PDBOUT, either add keyword RESIDUES or run a job using keyword RESIDUES to add PDB atom labels.)");
        std::printf("%18s\n", "(Keyword RESIDUES can only be used when one of MOZYME, LEWIS, CHARGES, or RESEQ is also present)");
        st = 2;                       // F90: return
        goto LendBody;
      }
    }
    if (keywrd.find(" ADD-H") != std::string::npos ||
        keywrd.find(" SITE=") != std::string::npos) nelecs = 0;
    if (keywrd.find(" 0SCF") != std::string::npos ||
        keywrd.find(" RESEQ") != std::string::npos ||
        keywrd.find(" ADD-H") != std::string::npos ||
        keywrd.find(" SITE=") != std::string::npos) {
      if (keywrd.find(" DISP") != std::string::npos) {
        l_control("0SCF", static_cast<int>(std::string("0SCF").size()), 1);
        l_control("PRT", static_cast<int>(std::string("PRT").size()), 1);
        post_scf_corrections(eat, false);
      }
      if (keywrd.find(" XYZ") != std::string::npos)
        line = " GEOMETRY IN CARTESIAN COORDINATES";
      else if (keywrd.find(" INT") != std::string::npos)
        line = " GEOMETRY IN MOPAC Z-MATRIX FORMAT";
      else
        line = " GEOMETRY OF SYSTEM SUPPLIED";
      if (prt_coords) std::printf("%s\n", line.c_str());
      xparam[1] = -1.0;
      if (keywrd.find(" OLDEN") != std::string::npos &&
          keywrd.find(" 0SCF") == std::string::npos) {
        // read in density so that charges can be calculated
        if (mozyme) {
          set_up_MOZYME_arrays();
        } else {
          p.resize(mpack);
          pa.resize(mpack);
          pb.resize(mpack);
        }
        den_in_out(0);
        if (moperr) { st = 2; goto LendBody; }   // F90: return
        if (mozyme) density_for_MOZYME(p, 0, nelecs / 2, pa);
      }
      if (keywrd.find(" ADD-H") != std::string::npos ||
          keywrd.find(" SITE=") != std::string::npos) {
        // Add RESEQ by default
        if (keywrd.find(" NORES") == std::string::npos &&
            (maxtxt == 26 || keywrd.find(" RESID") != std::string::npos) &&
            keywrd.find(" RESEQ") == std::string::npos)
          l_control("RESEQ", static_cast<int>(std::string("RESEQ").size()), 1);
      }
      if (prt_coords) geout(iw);
      if (keywrd.find(" AIGOUT") != std::string::npos) {
        std::printf("\n\n  GEOMETRY IN GAUSSIAN Z-MATRIX FORMAT\n");
        wrttxt(iw);
        geoutg(iw);
        std::printf("\n\n  GEOMETRY IN GAUSSIAN Z-MATRIX FORMAT\n");
        wrttxt(iw);
        geoutg(iw);
      } else if (mozyme ||
                 (keywrd.find(" PDBOUT") != std::string::npos ||
                  keywrd.find(" RESEQ") != std::string::npos ||
                  keywrd.find(" RESID") != std::string::npos)) {
        i = static_cast<int>(coorda.size());
        j = static_cast<int>(coord.size());
        if (i != j) {
          coorda.assign(coord.begin(), coord.end());
        }
        if (keywrd.find(" ADD-H") == std::string::npos &&
            keywrd.find(" SITE=") == std::string::npos) geochk();
        if (keywrd.find(" ADD-H") == std::string::npos &&
            keywrd.find(" SITE=") == std::string::npos &&
            keywrd.find(" RESEQ") != std::string::npos) {
          moperr = false;
          // goto 10 — restart job loop
          goto L10_next;
        }
        if (keywrd.find(" SITE=") != std::string::npos ||
            keywrd.find(" ADD-H") != std::string::npos) {
          moperr = false;
          if (keywrd.find(" ADD-H") != std::string::npos) {
            add_hydrogen_atoms();
            if (moperr) {
              // close(iarc, status="DELETE") — no-op
              st = 2;                 // F90: goto 101
              goto LendBody;
            }
            if (keywrd.find(" NORESEQ") != std::string::npos)
              update_txtatm(true, false);
            for (j = 1; j <= 6; ++j) {
              line = refkey[j];
              upcase(line, static_cast<int>(line.size()));
              i = static_cast<int>(line.find(" ADD-H"));
              if (i > 0) {
                refkey[j] = line.substr(0, i) + line.substr(i + 6);
                break;
              }
            }
            l_control("ADD-H", static_cast<int>(std::string("ADD-H").size()), -1);
            numat = natoms;
            numcal = numcal + 1;
            mopend("HYDROGEN ATOMS ADDED");
            moperr = false;
          }
          l_control("0SCF", static_cast<int>(std::string("0SCF").size()), 1);
          geochk();
          if (moperr) { st = 2; goto LendBody; }   // F90: goto 101
          if (coorda.size() >= coord.size()) {
            for (i = 1; i <= numat; ++i) {
              coorda[0][i] = coord[0][i];
              coorda[1][i] = coord[1][i];
              coorda[2][i] = coord[2][i];
              txtatm1[i] = txtatm[i];
            }
          }
          moperr = false;
          i = static_cast<int>(refkey[1].find("ADD-H"));
          if (i != 0) refkey[1] = refkey[1].substr(0, i) + refkey[1].substr(i + 5);
        }
        if (keywrd.find(" SITE=") != std::string::npos ||
            keywrd.find(" ADD-H") != std::string::npos ||
            keywrd.find(" RESEQ") != std::string::npos ||
            keywrd.find(" RESID") != std::string::npos)
          update_txtatm(true, false);  // switch to input labels after checks
        write_sequence();
        if (chanel_C::log) bridge_H();
        if (moperr) { st = 2; goto LendBody; }     // F90: goto 101
        // close(iarc); archive_fn = ...//"arc"; open(iarc); rewind — no-op.
        archive_fn = archive_fn.substr(0, archive_fn.size() - 3) + "arc";
        if (keywrd.find(" NOOPT") == std::string::npos &&
            keywrd.find(" OPT") == std::string::npos) {
          // Set all parameters to "opt"
          l = 0;
          for (i = 1; i <= numat; ++i) {
            for (j = 1; j <= 3; ++j) {
              l = l + 1;
              loc[1][l] = i;
              loc[2][l] = j;
            }
          }
        }
        geout(iw);
        if (keywrd.find(" PDBOUT") != std::string::npos) {
          // close(iarc); line = ...//"pdb"; open(iarc) — no-op mapping.
          line = archive_fn.substr(0, archive_fn.size() - 3) + "pdb";
          if (keywrd.find(" HTML") != std::string::npos) {
            for (i = static_cast<int>(line.size()); i >= 1; --i) {
              if (line[i - 1] == '/' || line[i - 1] == '\\') break;
            }
          }
          pdbout(iw);
        }
      } else {
        if (keywrd.find(" ADD-H") != std::string::npos) add_hydrogen_atoms();
        // open(iarc, archive_fn); rewind — no-op mapping.
        if (maxtxt == 26) {
          // PDB format, renumber atoms in atom label
          j = 1;
          for (i = 1; i <= numat; ++i) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%6s%5d%s", txtatm[i].substr(0, 6).c_str(),
                          i + j - 1, txtatm[i].substr(11).c_str());
            txtatm[i] = buf;
          }
        }
        geout(iw);
      }
      goto LendBody;                  // F90: go to 100
    }
    if (maxtxt == 26) compare_txtatm(moperr, moperr);
    if (moperr) {
      if (keywrd.find(" GEO-OK") != std::string::npos) {
        moperr = false;
      } else {
        goto LendBody;                // F90: goto 100
      }
    }
    //
    // Everything is ready - now set up the arrays used by the SCF, etc.
    //
    useps = false;
    iseps = (keywrd.find(" EPS=") != std::string::npos ||
             keywrd.find(" PKA") != std::string::npos);
    setup_mopac_arrays(1, 2);
    iseps = (keywrd.find(" EPS=") != std::string::npos);
    if (moperr) goto LendBody;        // F90: goto 100
    nw.resize(numat + 1);   // 1-based indices nw[1..numat]
    l = 1;
    for (i = 1; i <= numat; ++i) {
      nw[i] = l;
      l = l + ((nlast[i] - nfirst[i] + 1) * (nlast[i] - nfirst[i] + 2)) / 2;
    }
    //
    // CALCULATE THE ATOMIC ENERGY
    //
    if (sparkle) {
      atheat = 0.0;
      for (i = 1; i <= numat; ++i) {
        if (nat[i] > 56 && nat[i] < 72 && zs[nat[i]] < 0.1)
          atheat = atheat + eheat_sparkles[nat[i]];
        else
          atheat = atheat + eheat[nat[i]];
      }
    } else {
      atheat = 0.0;
      for (i = 1; i <= numat; ++i) atheat += eheat[nat[i]];
    }
    eat = 0.0;
    for (i = 1; i <= numat; ++i) eat += eisol[nat[i]];
    atheat = atheat - eat * fpc_9;
    atheat = atheat + C_triple_bond_C();
    rxn_coord = 1.e9;
    //
    // All data for the current job are now read in, and all parameters are
    // available in the arrays.  Now decide what type of calculation to do.
    //
    if (mozyme) {
      geochk();
      if (moperr) {
        if (keywrd.find(" PDBOUT") != std::string::npos) {
          archive_fn = archive_fn.substr(0, archive_fn.size() - 3) + "pdb";
          pdbout(iw);
        }
        // If the current calculation cannot continue, increment numcal to
        // indicate that a partial calculation has already been done.
        numcal = numcal + 1;
        // goto 10 — restart job loop
        goto L10_next;
      }
      if (keywrd.find(" RESID") != std::string::npos ||
          keywrd.find(" RESEQ") != std::string::npos ||
          keywrd.find(" ADD-H") != std::string::npos ||
          keywrd.find(" SITE=") != std::string::npos)
        update_txtatm(true, false);
      if (moperr) {
        numcal = numcal + 1;
        goto L10_next;                // goto 10
      }
      rapid = (keywrd.find(" RAPID") != std::string::npos ||
               keywrd.find(" LOCATE-TS") != std::string::npos);
      set_up_MOZYME_arrays();
      if (moperr) { st = 2; goto LendBody; }       // F90: goto 101
      if (keywrd.find(" RAPID") != std::string::npos) set_up_rapid("ON");
    }
    //
    // Dispatch to the requested calculation type
    //
    if (keywrd.find(" SADDLE") != std::string::npos) {
      to_screen(" Transition state geometry calculated using SADDLE");
      react1();
    } else if (keywrd.find(" STEP1") != std::string::npos) {
      to_screen(" Grid calculation");
      grid();
      iflepo = -1;                    // prevent printing of results
    } else if (latom != 0) {
      to_screen(" Path calculation");
      if (keywrd.find(" STEP") == std::string::npos ||
          keywrd.find(" POINT") == std::string::npos)
        paths();
      else
        pathk();
      iflepo = -1;
    } else if (keywrd.find(" FORCE") != std::string::npos ||
               keywrd.find(" IRC=") != std::string::npos ||
               keywrd.find(" THERM") != std::string::npos ||
               keywrd.find(" DFORCE") != std::string::npos) {
      to_screen(" Force constant calculation");
      last = 1;
      force();
      iflepo = -1;
    } else if (keywrd.find(" DRC") != std::string::npos ||
               keywrd.find(" IRC") != std::string::npos) {
      to_screen(" Reaction coordinate calculation");
      std::fill(na.begin(), na.end(), 0);
      if (react.empty()) react.assign(3 * numat + 1, 0.0);
      if (keywrd.find(" HTML") != std::string::npos) write_path_html();
      drc(react, react);
      iflepo = -1;
    } else if (keywrd.find(" NLLSQ") != std::string::npos) {
      std::printf(" Transition state refinement using NLLSQ\n");
      to_screen(" Transition state refinement using NLLSQ");
      nllsq();
    } else if (keywrd.find(" SIGMA") != std::string::npos) {
      std::printf(" Transition state refinement using SIGMA\n");
      to_screen(" Transition state refinement using SIGMA");
      powsq();
    } else if (keywrd.find(" 1SCF") != std::string::npos || nvar == 0) {
      iflepo = 1;
      iscf = 1;
      last = 1;
      i = static_cast<int>(keywrd.find(" GRAD"));
      grad.assign(nvar + 1, 0.0);
      numcal = numcal + 1;
      tim = second(1);
      wall_clock_0 = tim;   // timout() start markers (CPU ≈ wall-clock here)
      CPU_0 = tim;
      to_screen(" Single point calculation");
      compfg(xparam, true, escf, true, grad, i >= 0);
    } else if (keywrd.find(" LOCATE-TS") != std::string::npos) {
      if (!use_ref_geo) {
        if (keywrd.find(" LOCATE-TS(SET)") == std::string::npos) {
          line = " LOCATE-TS requires GEO_REF to be used";
          mopend(line);
          st = 2;                     // F90: return
          goto LendBody;
        }
      }
      Locate_TS_for_Proteins();
    } else if (keywrd.find(" REFINE-TS") != std::string::npos) {
      Refine_TS_for_Proteins();
    } else if (keywrd.find(" DFP") != std::string::npos ||
               keywrd.find(" FLEPO") != std::string::npos ||
               keywrd.find(" BFGS") != std::string::npos || nvar == 1) {
      std::printf(" Geometry optimization using BFGS\n");
      to_screen(" Geometry optimization using BFGS");
      flepo(xparam, nvar, escf);
    } else if (keywrd.find(" TS") != std::string::npos) {
      std::printf(" Transition state refinement using EF\n");
      to_screen(" Transition state refinement using EF");
      ef(xparam, escf);
    } else if (keywrd.find(" LBFGS") != std::string::npos ||
               (id == 0 && keywrd.find(" GEO_REF") != std::string::npos) ||
               (nvar > 1000 && keywrd.find(" EF") == std::string::npos)) {
      std::printf("\n%9s\n", "Geometry optimization using L-BFGS");
      to_screen(" Geometry optimization using L-BFGS");
      lbfgs(xparam.data(), escf);
    } else {
      std::printf("\n%9s\n", " 运行ef");
      ef(xparam, escf);
    }
    //
    // Calculation done, now print results
    //
    if (moperr) goto LendBody;        // F90: go to 100
    last = 1;
    if (iflepo >= 0) {
      wall_clock_1 = second(1);   // timout() end markers
      CPU_1 = wall_clock_1;
      writmo();
      if (moperr) goto LendBody;
      if (keywrd.find(" POLAR") != std::string::npos) {
        polar();
        if (moperr) goto LendBody;
      }
      if (keywrd.find(" STATIC") != std::string::npos) {
        numcal = numcal + 1;          // in case POLAR was also used
        static_polarizability();
        if (moperr) goto LendBody;
      }
      if (keywrd.find("PMEP") != std::string::npos) pmep();
      if (moperr) goto LendBody;
      if (keywrd.find(" ESP") != std::string::npos) {
        esp();
        if (moperr) goto LendBody;
      }
    }

  LendBody:
    // ===================== F90 label 100 tail =====================
    if (st == 2) break;               // goto 101: finish
    tim = tim + second(2);
    if (tim > 1.e7) tim = tim - 1.e7;
    // inquire(unit=ilog, opened=opend) — no-op mapping; LOG handling below.
    if (chanel_C::log && keywrd.find("LOG") == std::string::npos) {
      std::printf("\n == MOPAC DONE ==\n");
      // (F90: scan ilog for " COVALENT"/"charged"/"Type    Charge"/... and
      //  delete or keep the log file; kept as a comment stub.)
    }
    if (!gui) {
      p.clear();
      react.clear();
      {
        std::FILE* fp = std::fopen(end_fn.c_str(), "r");
        exists = (fp != nullptr);
        if (fp) std::fclose(fp);
      }
      if (exists) {
        std::remove(end_fn.c_str());  // open+close(status='delete')
      }
      itemp_1 = ncomments;
      if (keywrd.find(" ADD-H PDBOUT") == std::string::npos ||
          koment.find(" From PDB file") == std::string::npos) ncomments = 0;
      if (keywrd.find(" GEO_DAT") != std::string::npos) {
        i = static_cast<int>(keywrd.find(" GEO_DAT")) + 9;
        j = static_cast<int>(keywrd.find("\" ", i + 10)) + i + 9;
        line = "GEO_DAT=" + keywrd.substr(i, j - i);
        l_control(line, static_cast<int>(line.size()), 1);
      }
      delete_MOZYME_arrays();
      // goto 10 — next job
      continue;
    }
    break;                            // gui: fall through to label 101

  L10_next:
    // restart the job loop without the label-100 tail
    continue;
  }

  // ===================== F90 label 101 =====================
l101:
  if (!gui) {
    setup_mopac_arrays(0, 0);
    delete_MOZYME_arrays();
  }
  summary(" ", 1);
  if (tim > 1.e7) tim = tim - 1.e7;
  std::printf("\n\n\n TOTAL JOB TIME: %16.2f SECONDS\n", tim);
  std::printf("\n == MOPAC DONE ==\n");
  {
    std::string datestr;
    fdate(datestr);
    std::fprintf(stderr, "\n\n%20s\n", ("MOPAC Job: \"" + job_fn +
                                        "\" ended normally on " + datestr).c_str());
  }
  // Delete files that are definitely not wanted
  {
    std::FILE* fp = std::fopen(end_fn.c_str(), "r");
    exists = (fp != nullptr);
    if (fp) std::fclose(fp);
  }
  if (exists) std::remove(end_fn.c_str());
  // close(ir, status='delete') — no-op.
  jobnam = " ";
  return;
}
