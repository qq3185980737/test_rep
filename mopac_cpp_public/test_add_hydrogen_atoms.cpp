// test_add_hydrogen_atoms.cpp — numeric check of the pure-geometry routine
// add_a_generic_hydrogen_atom by including the translation unit.
//
// External routines not yet ported are given trivial stubs here.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "mod_atomradii.h"
#include "molkst_C.h"
#include "parameters_C.h"

using common_arrays_C::coord;
namespace mk = molkst_C;

// ---- stubs for external routines ---------------------------------------
double distance(int i, int j) {
    double dx = coord[1][i] - coord[1][j];
    double dy = coord[2][i] - coord[2][j];
    double dz = coord[3][i] - coord[3][j];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
double reada(const std::string& s, int) { return std::atof(s.c_str()); }
void dihed(const std::vector<std::vector<double>>&, int, int, int, int, double& sum) { sum = 0.0; }
void bangle(const std::vector<std::vector<double>>&, int, int, int, double& sum) { sum = 2.0; }
void upcase(std::string&, int) {}
void l_control(const std::string&, int, int) {}
void geochk() {}
void lewis(bool) {}
void set_up_dentate() {}
void mopend(const std::string&) {}

// Include the translation unit under test so its anonymous-namespace helpers
// (add_a_generic_hydrogen_atom, near_a_metal, ...) are reachable here.
#include "add_hydrogen_atoms.cpp"

int main() {
    using common_arrays_C::coord;
    using common_arrays_C::nat;
    using common_arrays_C::nbonds;
    using common_arrays_C::ibonds;

    // System: na=1 at origin, nb=2 along +x at 1.4 A, nc=0.
    mk::numat = 2;
    coord.assign(4, {0.0, 0.0, 0.0, 0.0});
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.4; coord[2][2] = 0.0; coord[3][2] = 0.0;
    nat.assign(5, 0);
    nat[1] = 6; nat[2] = 6;
    nbonds.assign(5, 0);
    ibonds.assign(16, std::vector<int>(5, 0));
    nbonds[1] = 1; ibonds[1][1] = 2;
    nbonds[2] = 1; ibonds[1][2] = 1;

    std::vector<int> metals(1, 0);
    int nmetals = 0;

    // Case A: angle = 180 deg -> H placed opposite nb, at bond_length=1.0.
    double ang = 3.14159265358979323846;
    double dih = 3.14159265358979323846;
    add_a_generic_hydrogen_atom(1, 2, 0, 1.0, ang, dih, metals, nmetals);
    double hx = coord[1][mk::numat], hy = coord[2][mk::numat], hz = coord[3][mk::numat];
    bool ok1 = std::fabs(hx + 1.0) < 1e-9 && std::fabs(hy) < 1e-9 && std::fabs(hz) < 1e-9;
    std::printf("case A: H=(%.6f, %.6f, %.6f)  expect (-1,0,0)  %s\n",
                hx, hy, hz, ok1 ? "PASS" : "FAIL");

    // Case B: angle=90 deg, dihedral=0, bond_length=1.0 -> H perpendicular
    // to nb, norm 1.
    mk::numat = 2;
    nbonds[1] = 1; nbonds[2] = 1;
    double ang2 = 1.5707963267948966;
    double dih2 = 0.0;
    add_a_generic_hydrogen_atom(1, 2, 0, 1.0, ang2, dih2, metals, nmetals);
    hx = coord[1][mk::numat]; hy = coord[2][mk::numat]; hz = coord[3][mk::numat];
    double norm = std::sqrt(hx * hx + hy * hy + hz * hz);
    bool ok2 = std::fabs(hx) < 1e-9 && std::fabs(norm - 1.0) < 1e-9;
    std::printf("case B: H=(%.6f, %.6f, %.6f) norm=%.6f  expect (0,~0.555,~0.832), norm 1  %s\n",
                hx, hy, hz, norm, ok2 ? "PASS" : "FAIL");

    return (ok1 && ok2) ? 0 : 1;
}
