// chanel_C.cpp — storage.
#include "chanel_C.h"
#include <vector>
namespace chanel_C {
int iw = 6;
int ir = 5;
int imep = 0; int iden = 4;
int igpt = 0; int iw0 = 0; int ires = 9; int iarc=8, ixyz=10;
bool log = false;
int ilog = 7;
std::string input_fn, output_fn, restart_fn, density_fn, log_fn, end_fn;
std::string archive_fn, brillouin_fn, esp_fn, ump_fn, mep_fn, pol_fn;
std::string gpt_fn, esr_fn, xyz_fn, syb_fn, cosmo_fn;
std::string job_fn;
int iscr = 12;
int lbfgs_it = 13;
int iump = 0;
std::vector<int> ioda(2001,0);
std::vector<int> ifilen(2001,0);
int irecst = 0;
int idaf = 17;
int irecln = 1023;
}
