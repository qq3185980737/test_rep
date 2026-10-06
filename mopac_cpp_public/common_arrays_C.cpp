// common_arrays_C.cpp — storage for the module-level data of Common_arrays_C.
#include "common_arrays_C.h"

namespace common_arrays_C {
std::vector<int> na_store;  // stored NA connectivity (getgeo IRC/DRC)
std::vector<double> fmatrx; // force matrix, packed lower triangle (force.F90)

std::vector<double> f;
std::vector<double> eigs; std::vector<double> eigb;

std::vector<int> nat;
std::vector<int> nbonds;
std::vector<int> labels;
std::vector<int> na;
std::vector<int> nb, nc;
std::vector<int> nfirst;
std::vector<double> uspd; // one-electron integrals (hcore.F90)
std::vector<int> nlast;
std::vector<int> breaks;
std::vector<std::vector<int>> ibonds;
std::vector<std::vector<int>> lopt;
std::vector<std::vector<int>> loc;
std::vector<double> atmass;
std::vector<double> xparam;
std::vector<std::vector<double>> coord;
std::vector<std::vector<double>> geo;
std::vector<std::vector<double>> coorda;
std::vector<std::vector<double>> geoa;
std::vector<std::vector<double>> fcint;
std::vector<std::vector<double>> break_coords;
std::vector<char> l_atom;
std::vector<double> hesinv;
std::vector<double> profil;
int time_start[8] = {};
std::vector<double> ch;
std::vector<std::string> txtatm;
std::vector<std::string> txtatm1;
std::vector<std::string> all_comments;
std::vector<std::vector<double>> tvec(4, std::vector<double>(4, 0.0));
std::string chains;
std::vector<double> Vab;
std::vector<int> cell_ijk;
std::vector<int> acceptor_a, acceptor_b;
std::vector<int> bonding_a_h, bonding_b_h;
std::vector<std::string> H_txt;      // hydrogen-bond diagnostic text (prt_hbonds)
std::vector<double> H_energy;        // hydrogen-bond energies (prt_hbonds)

std::vector<double> p;
std::vector<double> pdiag;
std::vector<std::vector<double>> c;
std::vector<double> pa;
std::vector<double> pb;
std::vector<std::vector<double>> cb;
std::vector<double> bondab;
std::vector<double> h;
std::vector<double> w;
std::vector<std::vector<double>> po(11, std::vector<double>(110, 0.0));
std::vector<double> wk;
std::vector<std::vector<int>> pibonds;
std::vector<double> dxyz;
std::vector<int> hblist;
std::vector<double> grad;
std::vector<double> gnext1, gmin1, aicorr;
std::vector<std::string> simbol;
std::vector<double> errfn, fb;
std::vector<double> T_range, HOF_tot, H_tot, Cp_tot, S_tot;
std::vector<std::string> pibonds_txt;
std::vector<double> q;
std::vector<int> ifact;
std::vector<int> i1fact;
std::vector<double> xparef;
std::vector<double> ptot2;
std::vector<int> nw;
}  // namespace common_arrays_C
