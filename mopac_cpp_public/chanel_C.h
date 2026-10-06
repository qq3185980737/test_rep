// chanel_C.h — C++ mapping of Fortran module "chanel_C".
#pragma once
#include <string>
namespace chanel_C {
extern int iw; extern int iw0; extern int ires; extern int iarc, ixyz;     // main output unit
extern int iump;  // grid surface data unit (grid.F90)
extern int ir;
extern int imep;    // input data unit (getsym.F90)
extern int iden;  // density file unit (pinout.F90)
extern int igpt;   // GRAPH file unit (mullik.F90)
extern bool log;   // true if a log file is open
extern int ilog;   // log output unit
extern std::string input_fn, output_fn, restart_fn, density_fn, log_fn, end_fn;
extern int lbfgs_it;  // iteration-summary unit (mainlb)
extern std::string archive_fn, brillouin_fn, esp_fn, ump_fn, mep_fn, pol_fn;
extern std::string gpt_fn, esr_fn, xyz_fn, syb_fn, cosmo_fn;
extern std::string job_fn;   // current job / dataset filename
extern int iscr;   // scratch file unit (geoutg.F90)
}
