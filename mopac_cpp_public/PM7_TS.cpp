// PM7_TS.cpp — C++ translation of MOPAC 2016 "PM7_TS" (writmo.F90 lines 1460-1471).
// PM7-TS transition-state post-processing: re-run SCF with TS parameters and
// print the final heat of formation, then return (writmo returns immediately).
#include "PM7_TS.h"

#include <cstdio>

#include "calpar.h"
#include "compfg.h"
#include "common_arrays_C.h"
#include "moldat.h"
#include "molkst_C.h"

void PM7_TS() {
    moldat(1);
    molkst_C::moperr = false;  // An error in moldat is not important here.
    calpar();
    double escf = 0.0;
    compfg(common_arrays_C::xparam, true, escf, true, common_arrays_C::grad, false);
    // F90: write(iw,'(4/10X,"FINAL HEAT OF FORMATION =",F17.5," KCAL/MOL", " =",F14.5," KJ/MOL")') escf, escf*4.184D0
    std::printf("\n\n\n\n%10sFINAL HEAT OF FORMATION =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                "", escf, escf * 4.184);
}
