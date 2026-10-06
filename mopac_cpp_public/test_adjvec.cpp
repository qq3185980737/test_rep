// test_adjvec.cpp — small CSR-style check of adjvec orthogonalisation.

#include <cmath>
#include <cstdio>
#include <vector>

#include "MOZYME_C.h"
#include "adjvec.h"

int main() {
    MOZYME_C::thresh = 0.001;  // cutoff = 0.01

    // One orbital per atom.
    std::vector<int> iorbs(3, 0);
    iorbs[1] = 1; iorbs[2] = 1;

    // LMO-a (lmoa=1): atoms 1,2.
    int lmoa = 1;
    std::vector<int> nnca(3, 0), nca_loc(3, 0), ncveca(3, 0), icveca(5, 0);
    nnca[1] = 0; nca_loc[1] = 2; ncveca[1] = 0;
    icveca[1] = 1; icveca[2] = 2;
    std::vector<double> cveca(11, 0.0);
    cveca[1] = 1.0; cveca[2] = 2.0;

    // LMO-b (lmob=1, only atom 1 initially).
    int lmob = 1, nnb = 1, nib = 100;
    std::vector<int> nncb(3, 0), ncb_loc(3, 0), ncvecb(3, 0), icvecb(10, 0);
    nncb[1] = 0; ncb_loc[1] = 1; ncvecb[1] = 0;
    icvecb[1] = 1;
    std::vector<double> cvecb(11, 0.0);

    std::vector<int> iused(3, 0);
    double beta = 0.5;
    double sumtot = 0.0;

    adjvec(cvecb, 10, icvecb, nib, nncb, ncb_loc, nnb, ncvecb, lmob,
           iorbs, cveca, 10, icveca, 4, nnca, nca_loc, 2, ncveca, lmoa,
           beta, iused, sumtot);

    // Expected: common atom 1 -> cvecb[1] = -0.5*1.0 = -0.5;
    // new atom 2 appended -> cvecb[2] = -0.5*2.0 = -1.0; ncb_loc=2; icvecb[2]=2.
    bool ok = true;
    ok &= std::fabs(cvecb[1] + 0.5) < 1e-12;
    ok &= std::fabs(cvecb[2] + 1.0) < 1e-12;
    ok &= (ncb_loc[1] == 2) && (icvecb[2] == 2);
    ok &= std::fabs(sumtot - 0.5) < 1e-12;
    std::printf("cvecb[1]=%.3f (expect -0.5)\n", cvecb[1]);
    std::printf("cvecb[2]=%.3f (expect -1.0)\n", cvecb[2]);
    std::printf("ncb_loc[1]=%d icvecb[2]=%d sumtot=%.3f\n",
                ncb_loc[1], icvecb[2], sumtot);
    std::printf("%s\n", ok ? "PASS" : "FAIL");

    // Small-beta early-return path.
    std::vector<double> cb2(11, 9.0);
    double st2 = 1.0;
    adjvec(cb2, 10, icvecb, nib, nncb, ncb_loc, nnb, ncvecb, lmob,
           iorbs, cveca, 10, icveca, 4, nnca, nca_loc, 2, ncveca, lmoa,
           1.0e-6, iused, st2);
    bool ok2 = (std::fabs(st2 - 1.0) < 1e-15);  // unchanged by early return
    std::printf("small-beta early return: %s\n", ok2 ? "PASS" : "FAIL");

    return (ok && ok2) ? 0 : 1;
}
