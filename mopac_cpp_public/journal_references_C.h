// journal_references_C.h — C++ translation of MOPAC 2016 "journal_references_C.F90".
// Per-element journal reference strings for each Hamiltonian. Fortran uses
// equivalence (refmn..refrm1) <-> allref(107,7); here the 7 arrays are kept
// independent (allref view omitted; semantic content identical).
#pragma once

#include <string>

namespace journal_references_C {
// 1-based [108]; entries without a data statement stay "" (e.g. refmn/refam/
// refpm3 are filled at runtime by the parameter loader in Fortran).
extern std::string refmn[108];
extern std::string refpm6[108];
extern std::string refam[108];
extern std::string refpm3[108];
extern std::string refmd[108];
extern std::string refrm1[108];
extern std::string refpm7[108];
}  // namespace journal_references_C
