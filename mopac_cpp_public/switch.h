// switch.h — C++ translation of MOPAC 2016 "switch.F90".
#pragma once

// Copies reference parameter sets (parameters_for_*_C) into the runtime
// arrays of parameters_C according to the active method flags, then runs
// the alpb_and_xfac_<model> fill, symmetrization and fractional-metal fixes.
void switch_method();
