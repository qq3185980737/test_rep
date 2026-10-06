// funcon_C.cpp — storage.
#include "funcon_C.h"
namespace funcon_C {
double pi = 3.14159265358979323846;
double a0 = 0.5291772083;
double ev = 27.2113834;
double fpc_9 = 23.060529;
double fpc_10 = 6.0221367e23;
double fpc_8 = 2.99792458e10;
double fpc_6 = 6.6260755e-27;
double fpc_2=14.399643;   // a0 * (one au in eV); au_to_kcal = fpc_9*fpc_2/a0 = 627.51
double fpc_1 = 1.602176462;  // Elementary charge, in Coulombs*10**19 (1998 Codata)
double fpc_5 = 1.9872065;    // R, gas constant (cal/mol/K)
double fpc_7 = 1.3806503e-16; // Boltzmann constant (erg/K)
}  // namespace funcon_C