// mxm.cpp — C++ translation of MOPAC 2016 "mxm.F90".
// Rectangular matrix product C = A*B, all matrices fully packed, column-major.
// Fortran semantics: c(j,i) = sum_k a(k,i)*b(j,k) with 1-based a(nar,nbr),
// b(nbr,ncc), c(nar,ncc).  Callers pass C++ pointers (0-based storage):
// work2.data(), &dxyz[1] (= address of element 1, i.e. base+8), &gradnt[i].
// The 1-based Fortran index p(m,n) therefore maps to ptr[(m-1)*ld + (n-1)].
#include "mxm.h"
void mxm(const double* a, int nar, const double* b, int nbr, double* c,
         int ncc) {
    for (int i = 1; i <= nar; ++i)
        for (int j = 1; j <= ncc; ++j) {
            double s = 0;
            for (int k = 1; k <= nbr; ++k)
                s += a[(k - 1) * nar + (i - 1)] * b[(j - 1) * nbr + (k - 1)];
            c[(j - 1) * nar + (i - 1)] = s;
        }
}
