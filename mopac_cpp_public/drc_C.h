// drc_C.h — C++ mapping of Fortran module "drc_C" (subset used).
#pragma once
#include <vector>
namespace drc_C {
extern std::vector<double> vref;     // reference velocities (DRC)
extern std::vector<double> vref0;    // initial reference velocities
extern std::vector<double> allxyz;   // all cartesian coords along DRC
extern std::vector<double> allvel;   // all velocities along DRC
extern std::vector<double> xyz3;     // current xyz (3 x natoms)
extern std::vector<double> vel3;     // current velocities
extern std::vector<double> allgeo;   // all internal coords along DRC
extern std::vector<double> geo3;     // current internal coords
extern std::vector<double> parref;   // reference parameters
extern double time;          // DRC elapsed time (femtoseconds)
}
