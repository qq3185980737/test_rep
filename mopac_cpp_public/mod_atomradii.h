// mod_atomradii.h — C++ mapping of Fortran module "mod_atomradii".
#pragma once
#include <vector>
namespace mod_atomradii {
// is_metal(z): true for elements treated as metals (negative VDW radius).
extern std::vector<bool> is_metal;
// atom_radius_covalent(z): covalent radius of element z; -1 marks metals.
extern std::vector<double> atom_radius_covalent;
// radius(z): per-atom covalent radius (1-based, filled by setup_mopac_arrays).
extern std::vector<double> radius;
}
