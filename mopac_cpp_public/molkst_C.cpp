// molkst_C.cpp — C++ translation of the Fortran module "molkst_C" (MOPAC 2016).
// Definitions for every member declared in molkst_C.h.
#include "molkst_C.h"

namespace molkst_C {

int numat = 0;
int numat_old = 0;
int norbs = 0;

int natoms = 0;
int lm61 = 0;
double enuclr = 0.0;
std::string keywrd = " ";
double mol_weight = 0.0;
int nopen = 0;
double rjkab1 = 0.0;
int msdel = 0;
double fract = 0.0;
int nclose = 0;
int nelecs = 0;
int nalpha = 0;
int nbeta = 0;
bool mozyme = false;
bool method_pm7 = false;
bool method_pm7_ts = false;
bool N_3_present = false;
bool Si_O_H_present = false;
bool method_pm6 = false;
bool method_am1 = false;
bool method_mndod = false;
bool is_PARAM = false;
bool sparkle = false;
bool method_rm1 = false;
bool method_mndo = false;
bool method_pm3 = false;
bool in_house_only = false;
int maxtxt = 0;
int nl_atoms = 0;
bool isok = false;
int maxatoms = 0;
bool moperr = false;
bool units = false;
bool Angstroms = true;
std::vector<std::string> refkey(6, " ");
std::string koment = " ";
std::string title = " ";
std::string line = " ";
int nbreaks = 0;
int ncomments = 0;
int id = 0;
int l11 = 0, l21 = 0, l31 = 0;

int mpack = 0;
bool lgpu = false;
int numcal = 0;
double clower = 0.0, cupper = 0.0, cutofp = 0.0;
int l1u = 0, l2u = 0, l3u = 0, l123 = 1;
int nvar = 0;
double pressure = 0.0, cosine = 0.0, escf = 0.0;
double tleft = 0.0, time0 = 0.0;
int nscf = 0, tdump = 0, iflepo = 0, last = 0;
int ndep = 0;
double efield[4] = {0.0, 0.0, 0.0, 0.0};
double E_disp = 0.0, E_hb = 0.0, E_hh = 0.0;
bool method_PM6_DH2X = false;
int N_Hbonds = 0;
bool method_PM6 = false;
bool method_pm6_d3_not_h4 = false;
int n2elec = 0, ispd = 0, itemp_2 = 0;
int itemp_1 = 0;
int step_num = 0;
std::string allkey = " ";
bool rhf = false, uhf = false;
std::string geo_dat_name = " ", geo_ref_name = " ";
double step = 0.0;
double Rab = 0.0;
int P_Hbonds = 0;
bool method_pm6_d3h4x = false, method_pm6_d3h4 = false, method_pm6_d3 = false;
bool method_pm6_dh_plus = false, method_pm6_dh2 = false, method_pm6_dh2x = false;
bool method_pm7_hh = false, method_pm7_minus = false;
double ux = 0.0, uy = 0.0, uz = 0.0;
int na1 = 0;
std::string jobnam = " ";
bool gui = false;
std::string errtxt = " ";
bool use_ref_geo = false;
double stress = 0.0, density = 0.0, elect = 0.0, atheat = 0.0;
std::string verson = " ";
int ijulian = 0, iscf = 0;
double hpress = 0.0, nsp2_corr = 0.0, Si_O_H_corr = 0.0, sum_dihed = 0.0;
double emin = 0.0;              // lowest energy found (compfg.F90)
int nalpha_open = 0, nbeta_open = 0;
bool prt_gradients = false, prt_coords = false, prt_cart = false,
     prt_pops = false, prt_charges = false;
bool prt_topo = false;
double gnorm = 0.0;
double zpe = 0.0;
std::string formula = " ";
double sz = 0.0, ss2 = 0.0;
int ltxt = 0;
double temp_1 = 0.0, temp_2 = 0.0;
int ilim = 0;
int num_threads = 0;
int run = 0;
double arc_hof_1 = 0.0, arc_hof_2 = 0.0;
int site_no = 0;
bool academic = false;
int num_bits = 64;

double param_constant = 1.0;
double trunc_1 = 7.0;
double trunc_2 = 0.22;
bool lxfac = false;
bool MM_corrections = false;
bool limscf = false;
int no_pKa = 0;
std::vector<int> mers(4, 0);
std::string program_name = "GUI             ";
}  // namespace molkst_C
