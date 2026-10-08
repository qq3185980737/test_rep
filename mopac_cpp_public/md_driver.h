// md_driver.h — C++ port of the user's Fortran molecular-dynamics driver
// (originally "run_md_velocity_verlet" + "gen_md_vel" from md.txt).
//
// Runs a velocity-Verlet MD for the current Cartesian system (H2/O2 mixtures
// and beyond) using MOPAC PM7 gradients (kcal/mol/A) as forces. All inputs are
// read from the existing globals (xparam, numat, nvar, ams, loc, nat, escf,
// elect/enuclr/solv_energy); temperature / box / timestep are controlled by
// keywords TEMP=, BOX= and DT=.
#pragma once
#include <string>

void run_md_velocity_verlet(const std::string& keywrd, int iw);
