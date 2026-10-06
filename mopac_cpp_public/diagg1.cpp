// diagg1.cpp — C++ translation of MOPAC 2016 "diagg1.F90".
// MOZYME significant occupied-virtual matrix-element construction; block-sparse
// extraction loops kept as skeleton (ijbo/nijbo external).

#include "diagg1.h"

#include <algorithm>
#include <cmath>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace MOZYME_C;
using namespace molkst_C;

namespace { int ijbo(int, int) { return -1; } void timer(const char*) {} }

void diagg1(const std::vector<double>&, int nocc, int,
            std::vector<double>&, std::vector<double>&,
            std::vector<char>& latoms,
            std::vector<std::vector<int>>&, std::vector<double>&,
            int, int& nij, int idiagg,
            std::vector<double>&, std::vector<double>& aocc,
            std::vector<double>&) {
    static int icalcn = 0;
    static double fref = 10.0, oldlim = 0.0, safety = 1.0;
    bool times = (keywrd.find(" TIMES") != std::string::npos);
    if (numcal != icalcn) {
        icalcn = numcal;
        fref = 10.0; safety = 1.0; oldlim = 0.0;
        if (keywrd.find(" OLDENS") != std::string::npos) fref = 0.0;
    }
    for (int k = 0; k < (int)aocc.size(); ++k) aocc[k] = 0.0;
    double cutlim = 1e-8;
    double cutoff = std::max(cutlim, tiny * 10.0 * cutlim);
    double flim = std::min(3.0, fref * 0.5);
    fref = 0.0;
    if (idiagg <= 5) cutoff = cutlim;
    sumt = 0.0; ijc = 0; tiny = 0.0;
    for (int i = 1; i <= numat; ++i) latoms[i] = 0;
    ovmax = tiny;
    (void)nocc; (void)nij; (void)cutoff; (void)flim; (void)times;
}
