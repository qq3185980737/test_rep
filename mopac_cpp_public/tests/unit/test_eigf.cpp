// test_eigf.cpp — verify eigenvectors_LAPACK against official water it1 Fock
#include "eigenvectors_LAPACK.h"
#include <cstdio>
#include <vector>
int main() {
    // official it1 f (lower-triangle packed, 1-based f[1..21])
    double fv[22] = {0,
        -30.219998, 0, -11.072378, 0, 0, -11.072378, 0, 0, 0, -11.072378,
        -6.165701, -5.885069, 0, 0, -6.304245, -6.160871, 1.469857,
        -5.695696, 0, -1.691744, -6.303056};
    std::vector<double> f(fv, fv + 22);
    std::vector<double> c(37, 0.0), eigs(7, 0.0);
    eigenvectors_LAPACK(c.data(), f.data(), eigs.data(), 6);
    std::printf("eigs: ");
    for (int i = 1; i <= 6; ++i) std::printf("%.4f ", eigs[i]);
    std::printf("\n");
    std::printf("c[1..6,1]= ");
    for (int i = 1; i <= 6; ++i) std::printf("%.4f ", c[i]);
    std::printf("\n");
    return 0;
}
