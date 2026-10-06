#include "ef_C.h"
namespace ef_C {
int nstep=0, negreq=0, iprnt=0, ef_mode=0;
double ddx=1.0, xlamd=0, xlamd0=0, skal=1, rmin=0, rmax=0, omin=0;
double x0=0, x1=0, x2=0; int iloop=0;
std::vector<std::vector<double>> alparm;
std::vector<std::vector<double>> hess;
std::vector<std::vector<double>> bmat;
std::vector<double> u, oldhss, oldu;
std::vector<std::vector<double>> pmat, uc, oldf, d, vmode;
std::vector<double> hessc;
}
