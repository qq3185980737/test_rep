// elemts_C.h — C++ mapping of Fortran module "elemts_C".
#pragma once
#include <string>
#include <vector>
namespace elemts_C {
// elemnt(z): chemical symbol (e.g. " C ", " N ") of element z (1-based).
extern std::vector<std::string> elemnt;
// cap_elemnt(z): capitalized element symbol used in PDB labels.
extern std::vector<std::string> cap_elemnt;
// atom_names(z): full atom name per element (len 12).
extern std::vector<std::string> atom_names;
}
