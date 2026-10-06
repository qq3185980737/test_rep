// scfcri.cpp — C++ translation of MOPAC 2016 "scfcri.F90".
#include "scfcri.h"
#include <cstring>
#include <cmath>
#include <string>
#include "molkst_C.h"
#include "chanel_C.h"
using namespace molkst_C;
using namespace chanel_C;
static int icalcn=0;
static double scfcrt=0, scfref=0;
static int findkw(const std::string& s, const char* kw){
    size_t n=std::strlen(kw);
    for (size_t i=0; i+n<=s.size(); ++i) if (s.compare(i,n,kw)==0) return (int)(i+1);
    return 0;
}
extern double reada(const std::string&, int);
void scfcri(double& selcon) {
    if (icalcn != numcal) {
        icalcn = numcal;
        scfcrt = 1e-2;
        int i = findkw(keywrd, " TS") + findkw(keywrd, " FORCETS") + findkw(keywrd, " IRC=");
        if (i != 0) scfcrt = 1e-3;
        bool precis = (findkw(keywrd, " PRECIS") != 0);
        i = findkw(keywrd, " RELSCF");
        if (i != 0) { scfcrt = reada(keywrd, i) * scfcrt; scfref = scfcrt; }
        i = findkw(keywrd, " SCFCRT");
        if (i != 0) { scfcrt = reada(keywrd, i); scfref = scfcrt; }
        if (precis) { scfcrt = scfcrt * 0.01; scfref = scfcrt; }
        if (findkw(keywrd, " POLAR") != 0 && scfref == 0.0) scfcrt = 1e-4;
        selcon = scfcrt;
    } else {
        if (std::fabs(efield[1]) + std::fabs(efield[2]) + std::fabs(efield[3]) > 1e-6)
            selcon = 1e-4;
        if (scfref != 0.0) selcon = scfref;
    }
}
