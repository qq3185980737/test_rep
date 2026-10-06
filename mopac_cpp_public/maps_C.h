#pragma once
#include <vector>
namespace maps_C {
extern double rxn_coord, rc_escf, ekin; extern int lparam, latom; extern std::vector<double> react; extern int kloop;
extern int lpara1, latom1, lpara2, latom2;   // reaction-path marker atoms (geout.F90)
// grid.F90 state.
extern double rxn_coord1, rxn_coord2;
extern int ione, ijlp, ilp, jlp, jlp1;
extern std::vector<double> surf;   // grid energies per point
}
