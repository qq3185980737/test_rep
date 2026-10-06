// dernvo.cpp — C++ translation of MOPAC 2016 "dernvo.F90" (complete).
// Analytical open-shell/CI Cartesian gradient driver. For each coordinate
// (one per ilast block) deri1 computes the non-relaxed (frozen cloud)
// contribution into dxyz and the non-relaxed Fock matrices fmooff/fmoon;
// deri2 then adds the electronic relaxation contribution to dxyzr.
// Layout: fmooff/fmoon are column-major supervectors (minear/ninear per
// coordinate, 1-based-padding); per-coordinate columns are copied into the
// single-column containers that deri1/deri2 expect.
#include "dernvo.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "deri0.h"
#include "deri1.h"
#include "deri2.h"
#include "meci_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;

namespace {
bool debug = false, dcar = false, large = false;
int minear = 0, ninear = 0, nvax = 0, icalcn = 0;
}  // namespace

void dernvo() {
    nbo[1] = nclose;
    nbo[2] = nopen - nclose;
    nbo[3] = norbs - nopen;
    minear = nbo[2] * nbo[1] + nbo[3] * nopen;
    ninear = nmos;
    int k = 0, n2 = 0, l = 0;
    for (int i = 1; i <= 3; ++i) {
        l = k + 1;
        k = k + nbo[i];
        int n1 = std::max(0, nelec - l);
        n2 = std::min(k, nelec + nmos) - l;
        if (n2 > n1) ninear += (n2 * (n2 + 1) - n1 * (n1 + 1)) / 2;
    }
    int j = n2 + 1;
    n2 = n2 + l;
    if (j > 0 && n2 < norbs)
        ninear += j * (norbs - n2);
    else if (ninear <= 0)
        ninear = 1;
    int i = std::max(norbs * norbs + 45 * lm61, (lab * (lab + 1)) / 2);
    i = std::max(i, n2elec + mpack);
    std::vector<double> work(i + 1, 0.0);
    std::vector<double> scalar(minear + 1, 0.0), diag(minear + 1, 0.0);
    std::vector<double> fmooff(2 * minear + 1, 0.0);
    std::vector<double> fmoon(nmos * norbs + 1, 0.0);

    if (icalcn != numcal) {
        icalcn = numcal;
        debug = keywrd.find("DERNVO") != std::string::npos;
        large = keywrd.find("LARGE") != std::string::npos;
        dcar = (keywrd.find("FORC") != std::string::npos) +
                   (keywrd.find("PREC") != std::string::npos) !=
               0;
        nvax = 3 * numat;
    }
    std::vector<double> dxyzr(nvax + 1, 0.0);

    // Scaling row factors to speed convergence of the relaxation procedure.
    std::vector<int> nbovec(4);
    for (int ii = 1; ii <= 3; ++ii) nbovec[ii] = nbo[ii];
    deri0(eigs, norbs, scalar, diag, fract, nbovec);
    fprintf(stderr, "[DSCAL] %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n",
        scalar[1], scalar[2], scalar[3], scalar[4], scalar[5], scalar[6], scalar[7], scalar[8]); fflush(stderr);
    fprintf(stderr, "[DDIAG] %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n",
        diag[1], diag[2], diag[3], diag[4], diag[5], diag[6], diag[7], diag[8]); fflush(stderr);

    // eigbb: scratch area for deri2 (F90 eigbb(nmos*norbs*20)).
    std::vector<std::vector<double>> eigbb(
        std::max(1, nmos * norbs * 20),
        std::vector<double>(std::max(1, ninear) + 1, 0.0));
    double sum = 0.0;
    bool relaxd = false;
    if (gnorm < 1.0 || dcar) {
        std::fill(dxyzr.begin(), dxyzr.end(), 0.0);
        relaxd = false;
    }
    for (i = 1; i <= nvax; ++i) sum += std::fabs(dxyzr[i]);
    relaxd = sum > 1e-7;

    // Work matrix for deri1 (2D), separate flat buffer for deri2.
    std::vector<std::vector<double>> work2d(
        (size_t)norbs + 1, std::vector<double>((size_t)norbs + 1, 0.0));

    int ilast = 0;
    for (;;) {
        int ifirst = ilast + 1;
        ilast = std::min(nvax, ilast + 1);
        j = 1 - minear;
        k = 1 - ninear;
        for (i = ifirst; i <= ilast; ++i) {
            k += ninear;
            j += minear;
            // Non-relaxed contribution (frozen cloud) into dxyz, Fock
            // matrices into fmooff(j) / fmoon(k).
            std::vector<double> fv(minear + 1, 0.0), fdv(nmos * norbs + 1, 0.0);
            for (int kk = 1; kk <= minear; ++kk) fv[kk] = fmooff[j + kk - 1];
            for (int kk = 1; kk <= ninear; ++kk) fdv[kk] = fmoon[k + kk - 1];
            double dgrad = 0.0;
            deri1(i, dgrad, fv, minear, fdv, scalar, work2d);
            dxyz[i] = dgrad;
            fprintf(stderr, "[DFMO] i=%d fv=%+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f fdv=%+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n",
                i, fv[1], fv[2], fv[3], fv[4], fv[5], fv[6], fv[7], fv[8], fdv[1], fdv[2], fdv[3], fdv[4], fdv[5], fdv[6]); fflush(stderr);
            for (int kk = 1; kk <= minear; ++kk) fmooff[j + kk - 1] = fv[kk];
            for (int kk = 1; kk <= ninear; ++kk) fmoon[k + kk - 1] = fdv[kk];
        }
        (void)chanel_C::iw;
        // Electronic relaxation contribution.
        if (!relaxd) {
            const int ncol = ilast - ifirst + 1;
            const int ncolmax =
                nmeci * (nmeci + 1) * 10 / std::max(1, ninear);
            std::vector<std::vector<double>> f2d(
                ncol + 1, std::vector<double>(minear + 1, 0.0));
            std::vector<std::vector<double>> fd2d(
                ncol + 1, std::vector<double>(ninear + 1, 0.0));
            std::vector<std::vector<double>> fci2d(
                std::max(ncol, ncolmax) + 1,
                std::vector<double>(ninear + 1, 0.0));
            for (int iv = 1; iv <= ncol; ++iv) {
                for (int kk = 1; kk <= minear; ++kk) f2d[iv][kk] = fmooff[kk];
                for (int kk = 1; kk <= ninear; ++kk)
                    fd2d[iv][kk] = fmoon[kk];
            }
            std::vector<double> dseg(ncol, 0.0);
            for (int iv = 1; iv <= ncol; ++iv)
                dseg[iv - 1] = dxyzr[ifirst + iv - 1];
            fprintf(stderr, "[DV] before deri2 ncol=%d minear=%d ninear=%d ncolmax=%d\n", ncol, minear, ninear, ncolmax); fflush(stderr);
            deri2(minear, f2d, fd2d, fci2d, ninear, ncol, dseg, 1e-5, diag,
                  scalar, work);
            fprintf(stderr, "[DV] after deri2\n"); fflush(stderr);
            for (int iv = 1; iv <= ncol; ++iv) {
                dxyzr[ifirst + iv - 1] = dseg[iv - 1];
                for (int kk = 1; kk <= minear; ++kk)
                    fmooff[kk] = f2d[iv][kk];
                for (int kk = 1; kk <= ninear; ++kk)
                    fmoon[kk] = fd2d[iv][kk];
            }
        }
        if (moperr) break;
        if (ilast < nvax) continue;
        break;
    }

    dxyz.resize(nvax + 1);
    for (i = 1; i <= nvax; ++i) dxyz[i] = dxyz[i] + dxyzr[i];
    fprintf(stderr, "[DVOUT] dxyz=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
        dxyz[1], dxyz[2], dxyz[3], dxyz[4], dxyz[5], dxyz[6], dxyz[7], dxyz[8], dxyz[9]); fflush(stderr);
    fprintf(stderr, "[DVOUT] dxyzr=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
        dxyzr[1], dxyzr[2], dxyzr[3], dxyzr[4], dxyzr[5], dxyzr[6], dxyzr[7], dxyzr[8], dxyzr[9]); fflush(stderr);
    if (relaxd) std::fill(dxyzr.begin(), dxyzr.end(), 0.0);
    (void)eigbb;
}
