// superd.cpp
#include "superd.h"
#include "parameters_C.h"
#include <cstdio>
using namespace parameters_C;
extern const char* elemnt(int);
extern int iw;
void superd(double* c, double* eigs, int norbs, int nelecs, int numat, int* nat) {
    int ihomo = nelecs/2;
    int isize = norbs;
    double alpha = (eigs[ihomo-1]+eigs[ihomo])/2.0;
    std::printf(" Mulliken electronegativity: %12.6f\n", -alpha);
    std::printf(" Parr & Pople absolute hardness: %12.6f\n", (eigs[ihomo]-eigs[ihomo-1])/2.0);
    std::printf(" Schuurmann MO shift alpha: %12.6f\n", alpha);
    std::printf(" Ehomo: %12.6f\n Elumo: %12.6f\n", eigs[ihomo-1], eigs[ihomo]);
    double totloc=0, totnuc=0;
    int i=1;
    for (int ii=1;ii<=numat;++ii) {
        if (nat[ii-1]==1) { ++i; }
        else {
            double cdens=0, deloc=0, denuc=0;
            int nb = natorb[nat[ii-1]];
            for (int k=i;k<=i+nb-1;++k)
                for (int j=1;j<=ihomo;++j){
                    double ck=c[(j-1)*norbs+(k-1)];
                    cdens += ck*ck;
                    deloc += ck*ck/(eigs[j-1]-alpha);
                }
            for (int j=ihomo+1;j<=isize;++j){
                for (int k=i;k<=i+nb-1;++k){
                    double ck=c[(j-1)*norbs+(k-1)];
                    denuc -= ck*ck/(eigs[j-1]-alpha);
                }
            }
            std::printf(" %s %d %13.6f%13.6f%13.6f\n", elemnt(nat[ii-1]), ii, 2*denuc, 2*deloc, -2*cdens);
            i += 4;
            totnuc += 2*denuc;
            totloc += 2*deloc;
        }
    }
    std::printf(" Total: %13.6f%13.6f\n", totnuc, totloc);
    double seltot=0; i=1;
    for (int ii=1;ii<=numat;++ii) {
        if (nat[ii-1]==1) { ++i; }
        else {
            double selpol=0; int nb=natorb[nat[ii-1]];
            for (int l=i;l<=i+nb-1;++l)
                for (int j=1;j<=ihomo;++j)
                    for (int k=ihomo+1;k<=isize;++k){
                        double clj=c[(j-1)*norbs+(l-1)], clk=c[(k-1)*norbs+(l-1)];
                        selpol += clj*clj*clk*clk/(eigs[k-1]-eigs[j-1]);
                    }
            selpol = -4*selpol;
            std::printf(" %s %d %13.6f\n", elemnt(nat[ii-1]), ii, selpol);
            seltot += selpol;
            i += nb;
        }
    }
    std::printf(" Total: %13.6f\n", seltot);
}
