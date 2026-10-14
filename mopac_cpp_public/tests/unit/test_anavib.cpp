// test_anavib.cpp
#pragma warning(disable: 4459)
#include <cmath>
#include <cstdio>
#include <vector>

#include "anavib.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "to_screen_C.h"

using namespace common_arrays_C;

int main() {
    // 2-atom H-H geometry: atom1 at origin, atom2 at 0.5 along x.
    molkst_C::numat = 2;
    nat.assign(3, 0); nat[1] = 1; nat[2] = 1;
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[1][2] = 0.5;  // x coords
    // elemnt[1] = H
    elemts_C::elemnt.assign(3, "  ");
    elemts_C::elemnt[1] = " H"; elemts_C::elemnt[2] = " H";

    const int n3 = 6, nv = 1;
    std::vector<double> eigs(n3 + 1, 0.0), dipt(n3 + 1, 0.0);
    eigs[1] = 1000.0; dipt[1] = 0.0;
    // vibs: mode 1 moves only dof 1 (atom 1 x).
    std::vector<std::vector<double>> vibs(n3 + 1, std::vector<double>(nv + 1, 0.0));
    vibs[1][1] = 1.0;
    // hess packed: hess[1] = H(1,1) = 5.0
    std::vector<double> hess(3 * molkst_C::numat * (3 * molkst_C::numat + 1) / 2 + 1, 0.0);
    hess[1] = 5.0;

    int npair = molkst_C::numat * (molkst_C::numat + 1) / 2;
    std::vector<double> rij(npair + 1, 0.0), f(npair + 1, 0.0);

    // symmetry & output arrays
    symmetry_C::jndex.assign(3, 1); symmetry_C::namo.assign(3, " A ");
    to_screen_C::travel.assign(nv + 1, 0.0);
    to_screen_C::redmas.assign(nv + 1, std::vector<double>(3, 0.0));
    to_screen_C::redmas[1][1] = 1.0; to_screen_C::redmas[1][2] = 1.0;
    chanel_C::iw = 6;

    anavib(eigs, dipt, n3, vibs, rij, nv, hess, f);

    bool ok = true;
    bool r1 = std::fabs(rij[1] - (0.5 + 1.0e-10)) < 1.0e-9;
    std::printf("rij[1]=%.6f (expect ~0.500000) %s\n", rij[1], r1 ? "PASS" : "FAIL");
    // After the selection sort, the picked pair's f(l) is zeroed to -1e-9.
    bool r2 = std::fabs(f[1] + 1.0e-9) < 1.0e-9;
    std::printf("f[1]=%.2e (after sort, expect ~-1e-9) %s\n", f[1], r2 ? "PASS" : "FAIL");
    ok = r1 && r2;
    return ok ? 0 : 1;
}
