// MOZYME_C.cpp — storage for MOZYME_C (subset).
#include "MOZYME_C.h"
namespace MOZYME_C {
std::vector<std::string> afn;
std::vector<std::vector<std::string>> atomname;
std::vector<int> n_add;
bool lijbo = false;
std::vector<std::vector<int>> nijbo;
bool direct = false;
bool semidr = false;
bool rapid = false;
double thresh = 0.0;
std::vector<double> partp;
std::vector<double> parth;
std::vector<double> partf;
std::vector<int> nbackb(5, 0);
int jatom = 0;
int iatom = 0;
int loop = 1;
int mxeno = 0;
std::vector<std::string> txeno;
std::vector<std::vector<int>> nxeno;
std::vector<std::string> allres;
int maxres = 0;
std::vector<std::vector<int>> bbone;
std::vector<std::vector<double>> angles;
std::vector<std::string> allr;
std::vector<int> ions;
std::vector<int> res_start;
std::vector<int> start_res;
std::vector<int> at_res;
std::vector<std::string> tyres;
std::vector<std::string> tyr;
int size_mres = 23;
bool odd_h = false;
std::vector<int> iorbs;
std::vector<double> ws;
int Lewis_tot=0;
std::vector<std::vector<int>> Lewis_elem;
std::vector<int> iz;
std::vector<int> ib;
int icharges=0;
std::vector<int> iij, numij, ijall, iijj;
int nres = 0;
int k = 0;
int nvirtual=0, noccupied=0;

double pmax = 0.0;
bool use_three_point_extrap = false;

std::vector<int> isort, nce, ncf, ncocc, ncvir, nnce, nncf, icocc, icvir;
std::vector<double> cocc, cvir;
int cocc_dim=0, cvir_dim=0, icocc_dim=0, icvir_dim=0;
double tiny=0.0, sumt=0.0, ovmax=0.0, energy_diff=0.0; int ijc=0;
double sumb=0.0; std::vector<double> fmo;
std::vector<std::vector<int>> ifmo;
std::vector<int> jopt;
int numred=0;
int norred=0, nelred=0;
int mode=0;
std::vector<int> kopt;
std::vector<double> part_dxyz, p1, p2, p3;
std::vector<int> idiag, nfmo;
int fmo_dim=0;
double refnuc=0.0;   // reference nuclear energy (hcore_for_MOZYME); global per MOZYME_C.h


double shift=0.0; double cutofs=49.0; int ij_dim=0; int morb=4; std::vector<double> tvec; int ipad2=1000000; int ipad4=1000000;
std::vector<std::vector<double>> geo_1, geo_2;
std::vector<double> dxyz_1, dxyz_2;
std::vector<double> xparam_1, xparam_2;
std::vector<double> partf_1, partp_1, parth_1;
std::vector<double> partf_2, partp_2, parth_2;
std::vector<double> f_1, p_1, f_2, p_2;
std::vector<int> nbonds_1, nbonds_2;
std::vector<std::vector<int>> ibonds_1, ibonds_2;
std::vector<int> icocc_1, icvir_1, ncocc_1, ncvir_1;
std::vector<int> nncf_1, nnce_1, ncf_1, nce_1;
std::vector<int> iij_1, iijj_1, ijall_1, numij_1, iorbs_1;
std::vector<std::vector<int>> nijbo_1;
std::vector<double> cocc_1, cvir_1;
std::vector<int> icocc_2, icvir_2, ncocc_2, ncvir_2;
std::vector<int> nncf_2, nnce_2, ncf_2, nce_2;
std::vector<int> iij_2, iijj_2, ijall_2, numij_2, iorbs_2;
std::vector<std::vector<int>> nijbo_2;
std::vector<double> cocc_2, cvir_2;
}  // namespace MOZYME_C
