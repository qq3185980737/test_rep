// deri1.cpp — C++ translation of MOPAC 2016 "deri1.F90" (complete).
// Computes the non-relaxed derivative of the non-variationally optimized
// wavefunction energy with respect to one Cartesian coordinate, and the
// non-relaxed Fock matrix derivative in MO basis (F/FD) as required by the
// relaxation section (deri2).
//
// Packed arrays (f/fd/scalar/fmat/hmat/wmat) follow the Fortran 1-based
// convention: element 0 is padding. xy follows the F90 column-major index
// xy(i,j,k,l) = ((l-1)*nmos + (k-1))*nmos*nmos + (j-1)*nmos + i (1-based
// semantic; physical element 0 padding). work comes in as
// vector<vector<double>> sized (norbs+1) x (norbs+1); the internal flat
// buffer wf is 1-based semantic (physical offset i == Fortran element i).

#include "deri1.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "derivs_C.h"
#include "dhcore.h"
#include "dfock2.h"
#include "dijkl1.h"
#include "funcon_C.h"
#include "helect.h"
#include "mecid.h"
#include "mecih.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mxm.h"
#include "mtxm.h"
#include "supdot.h"

using namespace common_arrays_C;
using namespace derivs_C;
using namespace funcon_C;
using namespace meci_C;
using namespace molkst_C;

static int icalcn = 0;
static double const_step = 0.0;  // F90 save
static bool debug = false;

void deri1(int number, double& grad,
           std::vector<double>& fv, int minear,
           std::vector<double>& fd,
           const std::vector<double>& scalar,
           std::vector<std::vector<double>>& work) {
    if (icalcn != numcal) {
        wmat.clear();
        wmat.resize(std::max(n2elec + mpack,
                     std::max((lab * (lab + 1)) / 2,
                     std::max(20 * mpack, maxci * (maxci + 1) / 2))));
        hmat.clear();
        hmat.resize(std::max((lab * (lab + 1)) / 2,
                     std::max(20 * minear, mpack)));
        fmat.clear();
        fmat.resize(std::max((lab * (lab + 1)) / 2,
                     std::max(20 * minear, mpack)));
        const_step = fpc_9;
        debug = (keywrd.find("DERI1") != std::string::npos);
        icalcn = numcal;
    }
    double step = 1.0e-3;

    // 2-point finite difference for the integral derivatives: stored in
    // hmat (1-electron) and wmat (2-electron), without dividing by step.
    int nati = (number - 1) / 3 + 1;
    int natx = number - 3 * (nati - 1);
    std::vector<double> coord_flat(3 * numat + 1);
    for (int a = 1; a <= numat; ++a)
        for (int x = 1; x <= 3; ++x)
            coord_flat[(a - 1) * 3 + (x - 1)] = coord[x-1][a];
    double enucl2 = 0.0;
    dhcore(coord_flat.data(), hmat.data(), wmat.data(), enucl2, nati, natx,
           step);
    step = 0.5 / step;

    // Non-relaxed Fock matrix derivative in AO basis, divided by step.
    for (int i = 1; i <= mpack; ++i) fmat[i] = hmat[i];
    dfock2(fmat.data(), p.data(), pa.data(), wmat.data(), numat,
           nfirst.data(), nlast.data(), nati);

    // Derivative of the SCF-only energy (before CI correction).
    grad = (helect(norbs, p.data(), hmat.data(), fmat.data()) + enucl2) *
           step;
    {
        double he = helect(norbs, p.data(), hmat.data(), fmat.data());
    }
    for (int i = 1; i <= mpack; ++i) fmat[i] *= step;

    // Column-major flat MO coefficients: c(ip,i) -> cflat[(i-1)*norbs+ip].
    std::vector<double> cflat(norbs * norbs + 1, 0.0);
    for (int i = 1; i <= norbs; ++i)
        for (int ip = 1; ip <= norbs; ++ip)
            cflat[(i - 1) * norbs + ip] = c[ip][i];

    // PART 1: work(n,n) = fmat * c(n,n). fmat is a 1-based-padding buffer
    // (Fortran element n sits at physical n), so hoff=1.
    std::vector<double> wf(norbs * norbs + 1, 0.0);
    for (int i = 1; i <= norbs; ++i)
        supdot(&wf[(i - 1) * norbs], fmat.data(), &cflat[(i - 1) * norbs],
               norbs, 1);

    // PART 2: f(ij) = (C'*WORK)(i,j), off-diagonal blocks only.
    // 1) open-closed, 2) virtual-closed, 3) virtual-open.
    // cflat/wf are 1-based-padding buffers (physical offset == Fortran
    // element number); matrix-segment starts therefore carry a +1 so the
    // 0-based access inside mtxm/mxm hits the Fortran 1-based element.
    int l = 1;
    if (nbo[2] != 0 && nbo[1] != 0) {
        mtxm(&cflat[nbo[1] * norbs + 1], nbo[2], &wf[1], norbs, &fv[l],
             nbo[1]);
        l += nbo[2] * nbo[1];
    }
    if (nbo[3] != 0 && nbo[1] != 0) {
        mtxm(&cflat[nopen * norbs + 1], nbo[3], &wf[1], norbs, &fv[l],
             nbo[1]);
        l += nbo[3] * nbo[1];
    }
    if (nbo[3] != 0 && nbo[2] != 0)
        mtxm(&cflat[nopen * norbs + 1], nbo[3],
             &wf[nbo[1] * norbs + 1], norbs, &fv[l], nbo[2]);
    for (int i = 1; i <= minear; ++i) fv[i] *= scalar[i];

    // PART 3: super-vector FD, CI-active diagonal blocks, unscaled.
    l = 1;
    int nend = 0, ninit = 0, n2 = 0;
    for (int loop = 1; loop <= 3; ++loop) {
        ninit = nend + 1;
        nend = nend + nbo[loop];
        int n1 = std::max(ninit, nelec + 1);
        n2 = std::min(nend, nelec + nmos);
        if (n2 < n1) continue;
        for (int i = n1; i <= n2; ++i) {
            if (i <= ninit) continue;
            mxm(&cflat[(i - 1) * norbs + 1], 1,
                &wf[(ninit - 1) * norbs + 1], norbs, &fd[l], i - ninit);
            l += i - ninit;
        }
    }
    int ncol = n2 - ninit + 1;
    if (ncol > 0 && n2 < norbs) {
        mtxm(&cflat[n2 * norbs + 1], norbs - n2,
             &wf[(ninit - 1) * norbs + 1], norbs, &fd[l], ncol);
        l += ncol * (norbs - n2);
    }

    // CI-active Fock eigenvalue derivatives, stored in FD (continued).
    int lcut = l;
    if (number == 3) {
    }
    for (int i = nelec + 1; i <= nelec + nmos; ++i) {
        double s = 0.0;
        for (int k = 1; k <= norbs; ++k)
            s += cflat[(i - 1) * norbs + k] * wf[(i - 1) * norbs + k];
        fd[l] = s;
        ++l;
    }

    // CI-active 2-electron integral derivatives, stored in xy. F90 passes
    // (wmat, fmat, hmat, fmat) as (w, cij, wcij, ckl); fmat is reused as
    // scratch inside dijkl1. cv is the active-MO slice c(ip, nelec+i).
    std::vector<std::vector<double>> cv(norbs, std::vector<double>(nmos, 0.0));
    for (int ip = 1; ip <= norbs; ++ip)
        for (int i = 1; i <= nmos; ++i)
            cv[ip - 1][i - 1] = c[ip][nelec + i];
    dijkl1(cv, norbs, nati, wmat, fmat, hmat, fmat, xy);
    // xy is indexed by the F90 column-major linear index minus one
    // (physical element 0 == xy(1,1,1,1)); scale the whole active block.
    for (int i = 0; i < nmos * nmos * nmos * nmos; ++i) xy[i] *= step;

    // Build the CI matrix derivative, stored in wmat. eigs starts at the
    // Fortran element lcut (fd is 1-based-padding), so pass offset lcut-1.
    double gse = 0.0;
    mecid(&fd[lcut], gse, eigb.data(), wf.data(), xy.data());
    mecih(wf.data(), wmat.data(), nmos, &lab, xy.data());
    {
    }

    // Non-relaxed CI contribution to the energy derivative. vectci is
    // 1-based-padding: state i starts at Fortran element j=(i-1)*lab+1.
    double sum = 0.0;
    for (int i = 1; i <= nstate; ++i) {
        int j = (i - 1) * lab + 1;
        // vectci is 1-based-padding (meci writes vectci[1..lab]); F90 reads
        // vectci(j) with j=1 => physical element 1. supdot's g[] is
        // physically 0-based (g[1] reads physical element 1), so pass
        // &vectci[j-1] to make g[1] == vectci(j) == physical j.
        supdot(wf.data(), wmat.data(), &vectci[j - 1], lab);
        for (int k = 1; k <= lab; ++k)
            sum += vectci[j + k - 1] * wf[k];
    }
    grad = (grad + sum / nstate) * const_step;

    // Copy the flat work buffer back into the out-parameter (1-based rows).
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= norbs; ++j)
            work[i][j] = wf[(j - 1) * norbs + i];
}
