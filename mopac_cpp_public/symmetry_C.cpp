// symmetry_C.cpp — storage.
#include "symmetry_C.h"
namespace symmetry_C {
std::vector<int> jndex;
std::vector<std::string> namo;
double cub[4][4] = {};
double elem[4][4][21] = {};
std::vector<std::vector<int>> jelem;
// group is defined (initialized to 20*5) in symmetry_tables.cpp.
int ielem[21] = {};
std::vector<int> jy(7, 0);          // Fortran: integer jy(6), 1-based (index 1..6)
std::string name;
int nclass=1;
int nirred=0;
std::vector<std::string> jx(21, " ");  // Fortran: character jx(20)*4, 1-based (index 1..20)
std::vector<int> locdep;
std::vector<int> locpar;
std::vector<int> idepfn;
std::vector<double> depmul;
int nsym=0;
std::vector<std::vector<int>> ipo;
double r[10][121] = {};
int nent=0;
std::string state_spin="SINGLET "; std::string state_Irred_Rep="    "; int state_QN=0;
int igroup = 0;  // point-group index (force.F90)
}

