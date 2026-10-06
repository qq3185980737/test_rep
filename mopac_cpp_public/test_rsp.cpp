// test_rsp.cpp — verify rsp (Jacobi) on 2x2 and 3x3 symmetric matrices.
#include "rsp.h"
#include <cstdio>
#include <cmath>

static bool near(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; }

int main() {
    bool ok = true;
    // 2x2: [[2,1],[1,2]] eigenvalues 3,1
    double a2[3] = {2.0, 1.0, 2.0};  // packed lower: (1,1),(2,1),(2,2)
    double root[2], vect[4];
    rsp(a2, 2, root, vect);
    if (!(near(root[0], 1.0) && near(root[1], 3.0))) {
        std::printf("FAIL 2x2 roots: %g %g\n", root[0], root[1]);
        ok = false;
    }
    // 3x3 diagonal: [[3,0,0],[0,2,0],[0,0,1]]
    double a3[6] = {3.0, 0.0, 2.0, 0.0, 0.0, 1.0};
    double root3[3], vect3[9];
    rsp(a3, 3, root3, vect3);
    if (!(near(root3[0], 1.0) && near(root3[1], 2.0) && near(root3[2], 3.0))) {
        std::printf("FAIL 3x3 roots: %g %g %g\n", root3[0], root3[1], root3[2]);
        ok = false;
    }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
