// meci_C.cpp — storage for the module-level data of the meci_C module.
#include "meci_C.h"

namespace meci_C {

std::vector<double> occa;

int is = 0;
int iiloop = 0;
int jloop = 0;

std::vector<std::vector<int>> ispqr;
int lab=0;
int nmeci=22; int nmos=0; int nstate=0; int nelec=-20; int nbo[4]={0,0,0,0};
std::vector<int> nalmat;
std::vector<std::vector<int>> microa;
std::vector<std::vector<int>> microb;
std::vector<double> conf;
std::vector<double> dijkl, xy;
std::vector<double> vectci;
std::vector<std::vector<double>> deltap;
std::vector<std::vector<double>> rjkaa;
std::vector<std::vector<double>> rjkab;
std::vector<double> eig;
std::vector<int> ispin;
int maxci = 20000;  // official readmo.F90: BITS64 -> maxci=20000 (BITS32 -> 5000)
std::vector<double> cdiag;
double cdiagi = 0.0;
double cif1 = 0.0;
double cif2 = 0.0;
int k = 0;
int dummy = 0;
int labsiz = 0;
int root_requested = 0;
int msdel = 0;
std::vector<double> spin;
std::vector<int> nfa;
std::vector<double> eiga;
std::vector<double> cimat;
}  // namespace meci_C
