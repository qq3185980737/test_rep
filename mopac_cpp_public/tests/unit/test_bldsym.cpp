// test_bldsym.cpp
#include <cmath>
#include <cstdio>

#include "bldsym.h"
#include "symmetry_C.h"

int main() {
    // cub is a Fortran 3x3 column-major array; set it to identity in flat layout.
    double* cb = &symmetry_C::cub[0][0];
    for (int i = 0; i < 9; ++i) cb[i] = 0.0;
    cb[0] = 1.0; cb[4] = 1.0; cb[8] = 1.0;
    // ioper=1: C2 about x -> diag(1,-1,-1)
    bldsym(1, 1);
    bool r1 = std::fabs(symmetry_C::elem[1][1][1] - 1.0) < 1e-9 &&
              std::fabs(symmetry_C::elem[2][2][1] + 1.0) < 1e-9 &&
              std::fabs(symmetry_C::elem[3][3][1] + 1.0) < 1e-9;
    std::printf("C2(x) diag=(%.1f,%.1f,%.1f) %s\n",
                symmetry_C::elem[1][1][1], symmetry_C::elem[2][2][1],
                symmetry_C::elem[3][3][1], r1 ? "PASS" : "FAIL");

    // ioper=8: C3 about z -> cos120=-0.5, sin120=0.866
    bldsym(8, 2);
    bool r2 = std::fabs(symmetry_C::elem[1][1][2] + 0.5) < 1e-9 &&
              std::fabs(symmetry_C::elem[2][1][2] - 0.8660254) < 1e-6 &&
              std::fabs(symmetry_C::elem[3][3][2] - 1.0) < 1e-9;
    std::printf("C3 top-left=(%.3f,%.3f) %.3f %s\n",
                symmetry_C::elem[1][1][2], symmetry_C::elem[2][1][2],
                symmetry_C::elem[3][3][2], r2 ? "PASS" : "FAIL");

    // ioper=20: infinite axis -> swap off-diagonal
    bldsym(20, 3);
    bool r3 = std::fabs(symmetry_C::elem[1][2][3] - 1.0) < 1e-9 &&
              std::fabs(symmetry_C::elem[2][1][3] - 1.0) < 1e-9;
    std::printf("infinite offdiag swap %s\n", r3 ? "PASS" : "FAIL");

    return (r1 && r2 && r3) ? 0 : 1;
}
