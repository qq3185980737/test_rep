// alpb_and_xfac_mndod.cpp — C++ translation of MOPAC 2016 "alpb_and_xfac_mndod.F90" (alpb_and_xfac_mndod).
#include "alpb_and_xfac.h"
#include "parameters_C.h"

void alpb_and_xfac_mndod() {
    using parameters_C::alpb;
    using parameters_C::xfac;
    alpb[11][1] = 1.05225212e0;
    xfac[11][1] = 1.00000000e0;
    alpb[11][6] = 1.05225212e0;
    xfac[11][6] = 1.00000000e0;
    alpb[12][1] = 1.35052992e0;
    xfac[12][1] = 1.00000000e0;
    alpb[12][6] = 1.48172071e0;
    xfac[12][6] = 1.00000000e0;
    alpb[16][12] = 1.48172071e0;
    xfac[16][12] = 1.00000000e0;
    alpb[13][1] = 1.38788000e0;
    xfac[13][1] = 1.00000000e0;
    alpb[13][6] = 1.38788000e0;
    xfac[13][6] = 1.00000000e0;
    alpb[13][13] = 1.38788000e0;
    xfac[13][13] = 1.00000000e0;
}
