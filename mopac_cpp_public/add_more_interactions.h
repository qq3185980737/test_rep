// add_more_interactions.h — C++ translation of MOPAC 2016
// "add_more_interactions.F90".
#pragma once

// Re-evaluate the interaction matrix during a MOZYME direct-SCF geometry
// optimisation and, if new close contacts appear, grow the packed arrays
// (pold, p, h, f and, when rapid, partp/parth/partf) by 20%.
void add_more_interactions();
