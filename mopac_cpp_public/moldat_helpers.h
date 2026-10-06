// moldat_helpers.h — C++ translation of moldat.F90 internal subprograms.
#pragma once

// setcup: determine CUTOFP and unit-cell counts (moldat.F90 lines 897-1035).
void setcup();

// setup_nhco: MM correction to -(C=O)-(NH)- linkages (moldat.F90 1250-1323).
void setup_nhco(int& ii);
