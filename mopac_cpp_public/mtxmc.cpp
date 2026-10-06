// mtxmc.cpp
#include "mtxmc.h"
#include "mxm.h"
// C(NAR,NAR) = (A(NBR,NAR))' * B(NBR,NAR), packed lower triangle
void mtxmc(const double* a, int nar, const double* b, int nbr, double* c) {
    int l = 0;
    for (int i = 1; i <= nar; ++i) {
        mxm(a + (i-1)*nbr, 1, b, nbr, c + l, i);
        l += i;
    }
}
