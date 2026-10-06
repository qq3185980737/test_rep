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
    fprintf(stderr, "[D1] enter number=%d lab=%d nmos=%d nelec=%d mpack=%d wmat=%zu\n", number, lab, nmos, nelec, mpack, wmat.size()); fflush(stderr);
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
        fprintf(stderr, "[D1] number=%d he=%.8f enucl2=%.8f step=%.6f grad_pre=%.8f\n",
                number, he, enucl2, step, (he+enucl2)*step); fflush(stderr);
        fprintf(stderr, "[D1] hmat1..8=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f fmat1..8=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f p1..8=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f he=%.4f enucl2=%.4f step=%.4f\n",
            hmat[1], hmat[2], hmat[3], hmat[4], hmat[5], hmat[6], hmat[7], hmat[8],
            fmat[1], fmat[2], fmat[3], fmat[4], fmat[5], fmat[6], fmat[7], fmat[8],
            p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], he, enucl2, step);
        fprintf(stderr, "[D1] fmat9..21=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f\n",
            fmat[9], fmat[10], fmat[11], fmat[12], fmat[13], fmat[14], fmat[15],
            fmat[16], fmat[17], fmat[18], fmat[19], fmat[20], fmat[21]);
        fprintf(stderr, "[D1] hmat9..21=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f\n",
            hmat[9], hmat[10], hmat[11], hmat[12], hmat[13], hmat[14], hmat[15],
            hmat[16], hmat[17], hmat[18], hmat[19], hmat[20], hmat[21]);
        fflush(stderr);
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
    if (number == 1) {
        fprintf(stderr, "[D1] wf_col4=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f col5=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f\n",
            wf[18+1], wf[18+2], wf[18+3], wf[18+4], wf[18+5], wf[18+6],
            wf[24+1], wf[24+2], wf[24+3], wf[24+4], wf[24+5], wf[24+6]); fflush(stderr);
        fprintf(stderr, "[D1] c_col4=%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f s4=",
            cflat[18+1], cflat[18+2], cflat[18+3], cflat[18+4], cflat[18+5], cflat[18+6]); fflush(stderr);
        double s4 = 0.0;
        for (int k = 1; k <= norbs; ++k) s4 += cflat[18+k]*wf[18+k];
        fprintf(stderr, "%.6f\n", s4); fflush(stderr);
    }

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
    // Scale F by the diagonal metric tensor 'scalar'.
    if (number==1) {
        fprintf(stderr, "[D1] nbo1=%d nbo2=%d nbo3=%d minear=%d nopen=%d nelec=%d nmos=%d\n",
                nbo[1], nbo[2], nbo[3], minear, nopen, nelec, nmos); fflush(stderr);
        fprintf(stderr, "[D1] f_raw:"); for (int ii=1; ii<=minear; ++ii) fprintf(stderr, " %+.6f", fv[ii]); fprintf(stderr, "\n"); fflush(stderr);
        fprintf(stderr, "[D1] fmat_all:"); for (int ii=1; ii<=mpack; ++ii) fprintf(stderr, " %+.6f", fmat[ii]); fprintf(stderr, "\n"); fflush(stderr);
        // expand fv[7] = sum_pq c(p,4) fmat(p,q) c(q,5)
        double f7 = 0.0;
        for (int p=1; p<=norbs; ++p) for (int q=1; q<=norbs; ++q) {
            int ij = (p>q) ? (p*(p-1))/2+q : (q*(q-1))/2+p;
            double fm = fmat[ij];
            if (fm == 0.0) continue;
            double c4 = cflat[(4-1)*norbs+p], c5 = cflat[(5-1)*norbs+q];
            double term = c4*fm*c5;
            f7 += term;
            fprintf(stderr, "[D1] fv7term p=%d q=%d ij=%d fm=%+.6f c4=%+.4f c5=%+.4f term=%+.6f\n",
                    p, q, ij, fm, c4, c5, term); fflush(stderr);
        }
        fprintf(stderr, "[D1] fv7sum=%+.6f\n", f7); fflush(stderr);
        // hmat-only and dfock2-only contributions to fv7 (number==1)
        double f7h = 0.0, f7d = 0.0;
        for (int p = 1; p <= norbs; ++p) for (int q = 1; q <= norbs; ++q) {
            int ij = (p > q) ? (p * (p - 1)) / 2 + q : (q * (q - 1)) / 2 + p;
            double c4 = cflat[(4 - 1) * norbs + p], c5 = cflat[(5 - 1) * norbs + q];
            f7h += c4 * hmat[ij] * c5;
            f7d += c4 * (fmat[ij] - hmat[ij]) * c5;
        }
        fprintf(stderr, "[D1] fv7_hmat=%+.6f fv7_dfock2=%+.6f (pre-step500)\n", f7h, f7d); fflush(stderr);
        // element-wise hmat*500 vs dfock2*500 for number==1 (O-x)
        if (number == 1) {
            int qs[6] = {1, 2, 3, 4, 5, 6};
            for (int qi = 0; qi < 6; ++qi) {
                int q = qs[qi];
                int ij = (2 > q) ? (2 * (2 - 1)) / 2 + q : (q * (q - 1)) / 2 + 2;
                fprintf(stderr, "[D1] fxy p=2 q=%d ij=%d hmat500=%+.6f dfock500=%+.6f fmat500=%+.6f\n",
                        q, ij, hmat[ij] * 500.0, (fmat[ij] - hmat[ij]) * 500.0, fmat[ij]);
            }
            // full non-zero map of fmat*500 + MO-basis (4,4),(5,5)
            fprintf(stderr, "[D1] hmat_raw:");
            for (int ii = 1; ii <= 21; ++ii)
                if (std::fabs(hmat[ii]) > 1e-6) fprintf(stderr, " %d:%+.6f", ii, hmat[ii]);
            fprintf(stderr, "\n");
            fprintf(stderr, "[D1] fmat_raw:");
            for (int ii = 1; ii <= 21; ++ii)
                if (std::fabs(fmat[ii]) > 1e-6) fprintf(stderr, " %d:%+.6f", ii, fmat[ii]);
            fprintf(stderr, "\n");
            for (int ii = 1; ii <= 21; ++ii)
                if (std::fabs(fmat[ii]) > 1e-6) fprintf(stderr, " %d:%+.6f", ii, fmat[ii]);
            fprintf(stderr, "\n");
            double f44 = 0, f55 = 0, f45 = 0;
            for (int p = 1; p <= norbs; ++p) for (int q = 1; q <= norbs; ++q) {
                int ij = (p > q) ? (p * (p - 1)) / 2 + q : (q * (q - 1)) / 2 + p;
                f44 += cflat[(4 - 1) * norbs + p] * fmat[ij] * cflat[(4 - 1) * norbs + q];
                f55 += cflat[(5 - 1) * norbs + p] * fmat[ij] * cflat[(5 - 1) * norbs + q];
                f45 += cflat[(4 - 1) * norbs + p] * fmat[ij] * cflat[(5 - 1) * norbs + q];
            }
            fprintf(stderr, "[D1] MO44=%+.6f MO55=%+.6f MO45=%+.6f (fmat500)\n", f44, f55, f45);
            fflush(stderr);
        }
        fprintf(stderr, "[D1] cmo4:"); for (int p=1; p<=norbs; ++p) fprintf(stderr, " %+.4f", cflat[(4-1)*norbs+p]); fprintf(stderr, "\n"); fflush(stderr);
        fprintf(stderr, "[D1] cmo5:"); for (int p=1; p<=norbs; ++p) fprintf(stderr, " %+.4f", cflat[(5-1)*norbs+p]); fprintf(stderr, "\n"); fflush(stderr);
    }
    for (int i = 1; i <= minear; ++i) fv[i] *= scalar[i];
    {
        fprintf(stderr, "[D1] fv number=%d:", number); fflush(stderr);
        for (int ii = 1; ii <= minear; ++ii) fprintf(stderr, " %+.6f", fv[ii]);
        fprintf(stderr, " scalar1..8:"); fflush(stderr);
        for (int ii = 1; ii <= minear; ++ii) fprintf(stderr, " %+.6f", scalar[ii]);
        fprintf(stderr, "\n"); fflush(stderr);
    }

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
        fprintf(stderr, "[D1] fd_pre l=%d lcut=%d\n", l, lcut), fflush(stderr);
    }
    for (int i = nelec + 1; i <= nelec + nmos; ++i) {
        double s = 0.0;
        for (int k = 1; k <= norbs; ++k)
            s += cflat[(i - 1) * norbs + k] * wf[(i - 1) * norbs + k];
        fd[l] = s;
        ++l;
    }
    if (number <= 2) {
        fprintf(stderr, "[D1] fd_all n=%d lcut=%d:", number, lcut); fflush(stderr);
        for (int ii = 1; ii <= lcut + nmos; ++ii) fprintf(stderr, " %+.6f", fd[ii]);
        fprintf(stderr, "\n"); fflush(stderr);
    }
    if (number == 3) {
        fprintf(stderr, "[D1] fd_pre2 l=%d lcut=%d\n", l, lcut); fflush(stderr);
        fprintf(stderr, "[D1] fd_lcut_terms i=5 cw:"); fflush(stderr);
        double s5 = 0.0;
        for (int k = 1; k <= norbs; ++k) {
            double cv = cflat[(5 - 1) * norbs + k], wv = wf[(5 - 1) * norbs + k];
            s5 += cv * wv;
            fprintf(stderr, " %d(%+.4f*%+.4f=%+.4f)", k, cv, wv, cv * wv);
        }
        fprintf(stderr, " sum5=%.6f\n", s5); fflush(stderr);
        fprintf(stderr, "[D1] fd_all lcut=%d:", lcut); fflush(stderr);
        for (int ii = 1; ii <= lcut + nmos; ++ii) fprintf(stderr, " %+.6f", fd[ii]);
        fprintf(stderr, "\n"); fflush(stderr);
    }

    // CI-active 2-electron integral derivatives, stored in xy. F90 passes
    // (wmat, fmat, hmat, fmat) as (w, cij, wcij, ckl); fmat is reused as
    // scratch inside dijkl1. cv is the active-MO slice c(ip, nelec+i).
    std::vector<std::vector<double>> cv(norbs, std::vector<double>(nmos, 0.0));
    for (int ip = 1; ip <= norbs; ++ip)
        for (int i = 1; i <= nmos; ++i)
            cv[ip - 1][i - 1] = c[ip][nelec + i];
    dijkl1(cv, norbs, nati, wmat, fmat, hmat, fmat, xy);
    fprintf(stderr, "[D1] after dijkl1 nati=%d\n", nati); fflush(stderr);
    // xy is indexed by the F90 column-major linear index minus one
    // (physical element 0 == xy(1,1,1,1)); scale the whole active block.
    for (int i = 0; i < nmos * nmos * nmos * nmos; ++i) xy[i] *= step;

    // Build the CI matrix derivative, stored in wmat. eigs starts at the
    // Fortran element lcut (fd is 1-based-padding), so pass offset lcut-1.
    double gse = 0.0;
    fprintf(stderr, "[D1] occa1..nmos=%+.6f %+.6f eigs1..nmos=%+.6f %+.6f\n",
            occa[1], occa[2], fd[lcut], fd[lcut+1]); fflush(stderr);
    mecid(&fd[lcut], gse, eigb.data(), wf.data(), xy.data());
    fprintf(stderr, "[D1] gse=%.6f eigb=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f wf1..10=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
        gse, eigb[1], eigb[2], eigb[3], eigb[4], eigb[5], eigb[6], eigb[7], eigb[8], eigb[9], eigb[10],
        wf[1], wf[2], wf[3], wf[4], wf[5], wf[6], wf[7], wf[8], wf[9], wf[10]); fflush(stderr);
    fprintf(stderr, "[D1] after mecid gse=%.6f\n", gse); fflush(stderr);
    mecih(wf.data(), wmat.data(), nmos, &lab, xy.data());
    {
        fprintf(stderr, "[D1] mecih-dbg number=%d fd_lcut..+3=%+.6f %+.6f %+.6f %+.6f wmat1..10=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f lcut=%d fd_sz=%zu xy0..15=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
                number, fd[lcut], fd[lcut+1], fd[lcut+2], fd[lcut+3],
                wmat[1], wmat[2], wmat[3], wmat[4], wmat[5], wmat[6], wmat[7], wmat[8], wmat[9], wmat[10],
                lcut, fd.size(),
                xy[0], xy[1], xy[2], xy[3], xy[4], xy[5], xy[6], xy[7],
                xy[8], xy[9], xy[10], xy[11], xy[12], xy[13], xy[14], xy[15]); fflush(stderr);
    }
    fprintf(stderr, "[D1] after mecih\n"); fflush(stderr);

    // Non-relaxed CI contribution to the energy derivative. vectci is
    // 1-based-padding: state i starts at Fortran element j=(i-1)*lab+1.
    fprintf(stderr, "[D1] sum-seg nstate=%d lab=%d vectci=%zu wf=%zu wmat=%zu\n",
        nstate, lab, vectci.size(), wf.size(), wmat.size()); fflush(stderr);
    double sum = 0.0;
    for (int i = 1; i <= nstate; ++i) {
        int j = (i - 1) * lab + 1;
        // vectci is 1-based-padding (meci writes vectci[1..lab]); F90 reads
        // vectci(j) with j=1 => physical element 1. supdot's g[] is
        // physically 0-based (g[1] reads physical element 1), so pass
        // &vectci[j-1] to make g[1] == vectci(j) == physical j.
        supdot(wf.data(), wmat.data(), &vectci[j - 1], lab);
        if (number == 3) {
            fprintf(stderr, "[D1] supdot_after wf1..4=%+.6f %+.6f %+.6f %+.6f\n",
                    wf[1], wf[2], wf[3], wf[4]); fflush(stderr);
            // manual check: g=vectci[1..4], h=wmat packed lower-tri physical 0
            fprintf(stderr, "[D1] manual g1..4=%+.8f %+.8f %+.8f %+.8f\n",
                    vectci[1], vectci[2], vectci[3], vectci[4]); fflush(stderr);
            fprintf(stderr, "[D1] manual h0..9=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
                    wmat[0], wmat[1], wmat[2], wmat[3], wmat[4],
                    wmat[5], wmat[6], wmat[7], wmat[8], wmat[9]); fflush(stderr);
            fprintf(stderr, "[D1] manual s3terms t1=%.6f t2=%.6f t3=%.6f t4=%.6f h8=%.6f\n",
                    vectci[1]*wmat[3], vectci[2]*wmat[4], vectci[3]*wmat[5],
                    vectci[4]*wmat[8], wmat[8]); fflush(stderr);
            double ms2 = vectci[1]*wmat[1] + vectci[2]*wmat[2] + vectci[3]*wmat[4] + vectci[4]*wmat[7];
            double ms3 = vectci[1]*wmat[3] + vectci[2]*wmat[4] + vectci[3]*wmat[5] + vectci[4]*wmat[8];
            fprintf(stderr, "[D1] manual s2=%.8f s3=%.8f\n", ms2, ms3); fflush(stderr);
        }
        for (int k = 1; k <= lab; ++k)
            sum += vectci[j + k - 1] * wf[k];
    }
    if (number == 3) {
        fprintf(stderr, "[D1] vectci1..4=%+.6f %+.6f %+.6f %+.6f wmat00=%+.6f wmat(1,1)=%+.6f (4,1)=%+.6f (4,4)=%+.6f\n",
                vectci[1], vectci[2], vectci[3], vectci[4],
                wmat[0], wmat[0], wmat[6], wmat[9]); fflush(stderr);
    }
    {
        fprintf(stderr, "[D1] sumin number=%d vectci0..3=%+.6f %+.6f %+.6f %+.6f wf1..4=%+.6f %+.6f %+.6f %+.6f wmat1..16=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
                number, vectci[0], vectci[1], vectci[2], vectci[3],
                wf[1], wf[2], wf[3], wf[4],
                wmat[1], wmat[2], wmat[3], wmat[4], wmat[5], wmat[6], wmat[7], wmat[8],
                wmat[9], wmat[10], wmat[11], wmat[12], wmat[13], wmat[14], wmat[15], wmat[16]); fflush(stderr);
    }
    grad = (grad + sum / nstate) * const_step;
    {
        fprintf(stderr, "[D1] number=%d grad_before=%.8f sum=%.8f nstate=%d grad_total=%.8f\n",
                number, grad, sum, nstate, grad); fflush(stderr);
    }

    // Copy the flat work buffer back into the out-parameter (1-based rows).
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= norbs; ++j)
            work[i][j] = wf[(j - 1) * norbs + i];
}
