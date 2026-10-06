// alpb_and_xfac.h — C++ translation of the alpb_and_xfac_* subroutines
// (MOPAC 2016: alpb_and_xfac_am1/pm3/mndo/mndod.F90 and the pm6/pm7/pm7_TS
// copies embedded in parameters_for_PMx_C.F90). Each fills parameters_C::alpb
// and parameters_C::xfac (1..100).
#pragma once

void alpb_and_xfac_am1();
void alpb_and_xfac_pm3();
void alpb_and_xfac_mndo();
void alpb_and_xfac_mndod();
void alpb_and_xfac_pm6();
void alpb_and_xfac_pm7();
void alpb_and_xfac_pm7_TS();
