// extvdw_for_MOZYME.h — C++ translation of "extvdw_for_MOZYME" (from
// geochk.F90, lines 1490-1579).  Fills per-atom van der Waals radii,
// applying keyword "METAL" and "VDWM(:sym=value:...)" overrides.
#pragma once
#include <vector>

void extvdw_for_MOZYME(std::vector<double>& radius,
                       const std::vector<double>& refvdw);
