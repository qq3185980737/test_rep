// deri2.cpp — C++ translation of MOPAC 2016 "deri2.F90" (complete).
// Relaxation part of non-variational energy derivatives. Three steps:
//   STEP 1: initial orthonormal basis B from F (deri21).
//   STEP 2: iterative relaxation in the diagonal metric SCALAR; project the
//           electronic Hessian (D-A) onto the growing basis, invert, and
//           extend the basis by the largest residual vector until
//           convergence or storage limits.
//   STEP 3: per geometric variable, unpack the CI-active MO derivatives
//           (deri23), relax the 2-electron integrals (dijkl2), build the
//           CI matrix derivative (mecid/mecih) and accumulate the
//           relaxation correction into dxyzr.
//
// Layout: f/fd/fci are [col][row] 1-based containers; internal flat buffers
// use the 1-based-padding convention (physical offset == Fortran element
// number, element 0 unused). b/ab/fb/bcoef are column-major supervectors in
// the same convention, so every routine is called with the F90 element
// offset (b(1,j) -> &b[(j-1)*minear+1], etc.).
#include "deri2.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "common_arrays_C.h"
#include "deri21.h"
#include "deri22.h"
#include "deri23.h"
#include "dijkl2.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "mecid.h"
#include "mecih.h"
#include "molkst_C.h"
#include "mxm.h"
#include "mxmt.h"
#include "mtxm.h"
#include "osinv.h"
#include "supdot.h"

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;
using funcon_C::fpc_9;

// xy(a,b,c,d) column-major index (F90 xy(nmos,nmos,nmos,nmos)).
static int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

void deri2(int minear, std::vector<std::vector<double>>& f,
           std::vector<std::vector<double>>& fd,
           std::vector<std::vector<double>>& fci, int ninear, int nvar_nvo,
           std::vector<double>& dxyzr, double throld,
           std::vector<double>& diag, std::vector<double>& scalar,
           std::vector<double>& work) {
    static int icalcn = 0;
    static int maxite = 0, ifirst = 0, limci = 0, ib = 0;
    static double cnst = 0.0;
    static bool debug = false;
    static std::vector<double> b_buf, ab_buf, fb_buf;
    // F90 ab is a module-level buffer shared across deri2 calls within the
    // same numcal; deri23/dijkl2 overwrite only the active block, so the
    // frozen part keeps the previous step's residues. cmo2d plays that role.
    static std::vector<std::vector<double>> cmo2d;

    if (icalcn != numcal) {
        limci = nmos * norbs * 20 / std::max(1, ninear);
        ib = std::max(lab, 20);
        debug = keywrd.find(" DERI2") != std::string::npos;
        cnst = fpc_9;
        icalcn = numcal;
        ifirst = std::min(nvar_nvo, 1 + maxite / 4);
        int k = std::max((lab * (lab + 1)) / 2, 20 * mpack);
        int i = std::max(k, (maxci * (maxci + 1)) / 2);
        b_buf.assign((size_t)minear * (i / minear) + 1, 0.0);
        i = std::max(ninear, i / minear);
        fb_buf.assign((size_t)std::max(i, nvar_nvo * ib) + 1, 0.0);
        i = std::max(20 * minear, (lab * (lab + 1)) / 2);
        i = std::max(i, mpack);
        ab_buf.assign((size_t)minear * (i * 2 / minear) + 1, 0.0);
        int csz = std::max(minear, norbs) + nmos + 2;
        cmo2d.assign((size_t)csz, std::vector<double>((size_t)csz, 0.0));
        maxite = std::min(std::min(60, (int)std::sqrt((double)nmeci * nmeci * nmeci)),
                          10000 * 2 / nvar_nvo);
        maxite = std::min(maxite, i / minear);
    }

    // 2D work matrix used by deri22 (same buffer for all columns).
    std::vector<std::vector<double>> work2d(
        (size_t)norbs + 1, std::vector<double>((size_t)norbs + 1, 0.0));
    std::vector<double> foc2((size_t)norbs * norbs + 1, 0.0);

    // Flatten the [col][row] containers to column-major 1-based-padding.
    std::vector<double> fflat((size_t)minear * nvar_nvo + 1, 0.0);
    std::vector<double> fdflat((size_t)ninear * nvar_nvo + 1, 0.0);
    const int ncolmax = nmeci * (nmeci + 1) * 10 / std::max(1, ninear);
    std::vector<double> fciflat((size_t)ninear * ncolmax + 1, 0.0);
    for (int iv = 1; iv <= nvar_nvo; ++iv) {
        for (int i = 1; i <= minear; ++i) fflat[(iv - 1) * minear + i] = f[iv][i];
        for (int i = 1; i <= ninear; ++i) fdflat[(iv - 1) * ninear + i] = fd[iv][i];
    }

    // STEP 1: initial orthonormal basis.
    int ilast = 0;
    deri21(&fflat[1], nvar_nvo, minear, ifirst, work.data(),
           &work[nvar_nvo * nvar_nvo], &b_buf[1], ilast);
    bool lbab = false;
    ifirst = 1;  // nbsize == 0 at this point
    ilast = ilast;
    std::vector<char> lconv(nvar_nvo + 1, 0);
    std::vector<double> bab(70 * 70 + 1, 0.0), babinv(4900, 0.0);
    std::vector<double> bcoef(1001, 0.0);  // 1-based padding (F90 bcoef(1000))
    bool fail = false;
    int nres = 0, nadd = 0;
    double test2 = 0.0, deter = 0.0;

    // STEP 2: relaxation loop.
    for (;;) {
        if (ilast > limci || ilast > ib) {
            // Analytical derivatives not possible; NOANCI path requested in F90.
            fail = true;
            goto done;
        }
        for (int j = ifirst; j <= ilast; ++j) {
            deri22(c, b_buf, j, work2d, foc2, ab_buf, j, minear, fciflat, j,
                   w, diag, scalar, ninear);
            mxm(&ab_buf[(j - 1) * minear + 1], 1, &b_buf[1], minear,
                &bab[(j - 1) * 70 + 1], ilast);
            for (int kk = 1; kk <= ifirst - 1; ++kk)
                bab[(kk - 1) * 70 + j] = bab[(j - 1) * 70 + kk];
        }
    invert:
        int l = 0;
        for (int j = 1; j <= ilast; ++j) {
            for (int kk = 1; kk <= ilast; ++kk)
                babinv[l + kk - 1] = bab[(kk - 1) * 70 + j];
            l += ilast;
        }
        osinv(babinv.data(), ilast, deter);
        if (std::fabs(deter) < 1e-20) {
            if (ilast != 1) {
                // BAB singular; stop relaxation here (F90 prints to iw).
                lbab = true;
                --ilast;
                goto invert;
            }
            lbab = true;
            --ilast;
            goto invert;
        }
        if (!lbab)
            mtxm(&fflat[1], nvar_nvo, &b_buf[(ifirst - 1) * minear + 1],
                 minear, &fb_buf[ifirst], ilast - ifirst + 1);
        if (ilast != 0)
            mxmt(babinv.data(), ilast, &fb_buf[1], ilast, &bcoef[1],
                 nvar_nvo);
        if (lbab) goto converged;

        // Select next basis vector as the largest residual.
        nres = 0;
        test2 = 0.0;
        for (int ivar = 1; ivar <= nvar_nvo; ++ivar) {
            if (lconv[ivar]) continue;
            mxm(&ab_buf[1], minear, &bcoef[ilast * (ivar - 1) + 1], ilast,
                work.data(), 1);
            double test = 0.0;
            for (int i = 1; i <= minear; ++i) {
                work[i - 1] = fflat[(ivar - 1) * minear + i] - work[i - 1];
                test = std::max(std::fabs(work[i - 1]), test);
            }
            test2 = std::max(test2, test);
            if (test <= throld) {
                lconv[ivar] = 1;
                if (nvar_nvo == 1) goto converged;
                continue;
            } else if (ilast + nres == maxite) {
                fail = (nres == 0);
                break;
            } else if (ilast + nres >= maxite - 1) {
                if (test <= std::max(0.01, throld * 2)) {
                    lconv[ivar] = 1;
                    continue;
                }
            } else {
                ++nres;
                for (int i = 1; i <= minear; ++i)
                    ab_buf[(ilast + nres - 1) * minear + i] = work[i - 1];
            }
        }
        if (nres == 0) goto converged;
        ifirst = ilast + 1;
        deri21(&ab_buf[(ifirst - 1) * minear + 1], nres, minear, nres,
               work.data(), &work[nres * nres], &b_buf[(ifirst - 1) * minear + 1],
               nadd);
        ilast += nadd;
    }

converged:
    if (fail) {
        // Analytical derivatives too inaccurate; job stopped (F90 prints).
        goto done;
    }
    // Unscaled solution supervectors, stored in f.
    if (ilast != 0) mxm(&b_buf[1], minear, &bcoef[1], ilast, &fflat[1], nvar_nvo);
    for (int j = 1; j <= nvar_nvo; ++j)
        for (int i = 1; i <= minear; ++i)
            fflat[(j - 1) * minear + i] *= scalar[i];
    // Fock matrix diagonal blocks over CI-active MOs, stored in fb.
    if (ilast != 0)
        mxm(&fciflat[1], ninear, &bcoef[1], ilast, &fb_buf[1], nvar_nvo);

    // STEP 3: final loop on geometric variables.
    {
        std::vector<double> cmoflat((size_t)(lab * (lab + 1)) / 2 + 1, 0.0);
        for (int ivar = 1; ivar <= nvar_nvo; ++ivar) {
            // CI-active MO derivatives into MO basis (cmo2d), eigenvalues
            // derivatives into bcoef.
            deri23(fflat, fdflat, eigs, fb_buf, cmo2d, bcoef, minear, ninear,
                   ivar);
            // F90 cmo is (norbs,norbs) written through the shared ab buffer
            // with 6-column stride: cmo(i,j) -> ab linear (j-1)*norbs+i.
            // Copy the fresh block over ab_buf; frozen positions (linear
            // 37..54, i.e. ab cols 5..7) keep deri22 residues across calls.
            for (int r = 0; r < norbs; ++r)
                for (int c = 0; c < norbs; ++c)
                    ab_buf[c * norbs + r + 1] = cmo2d[r][c];
            int iindex = (norbs * nelec + 1) % minear;
            int jindex = (norbs * nelec + 1) / minear + 1;
            if (iindex == 0) {
                iindex = minear;
                --jindex;
            }
            // 2-electron integral relaxation over CI-active MOs, in xy.
            {
                std::vector<std::vector<double>> dc(
                    (size_t)norbs + 1, std::vector<double>((size_t)norbs + 1, 0.0));
                // F90 dijkl2 receives dc(norbs,nmos): column j = MO derivative
                // vector, i.e. dc(p,j). The actual argument is ab(iindex,jindex)
                // with ab(minear,*), and the dc dummy has leading dim norbs,
                // so in column-major layout (1-based, Fortran):
                //   dc(p,j) -> ab linear (jindex-1)*minear + iindex
                //                          + (j-1)*norbs + (p-1)
                //          -> cmo linear offset X = that linear index - 1,
                // where cmo is the deri23 output written through the same ab
                // buffer (leading dim norbs). cmo2d here is [row][col] 0-based
                // with linear = col*norbs+row. C++ dijkl2 reads dc[j-1] (a row)
                // as column j, so:
                for (int j = 1; j <= nmos; ++j)
                    for (int p = 1; p <= norbs; ++p) {
                        int L = (jindex - 1) * minear + iindex + (j - 1) * norbs + (p - 1);
                        dc[j - 1][p - 1] = ab_buf[L];
                    }
                dijkl2(dc);
                // F90 dijkl2 overwrites ab(iindex,jindex) onward; those
                // residues are read back by the NEXT ivar's dc assembly
                // (frozen part), so write the result back into cmo2d at the
                // same linear offsets.
                for (int jj = 1; jj <= nmos; ++jj)
                    for (int pp = 1; pp <= norbs; ++pp) {
                        int L2 = (jindex - 1) * minear + iindex + (jj - 1) * norbs + (pp - 1);
                        ab_buf[L2] = dc[jj - 1][pp - 1];
                    }
            }
            // Build CI matrix derivative (stored in cmoflat packed).
            double gse = 0.0;
            mecid(&bcoef[nelec + 1], gse, &work[lab], work.data(), xy.data());
            mecih(work.data(), cmoflat.data(), nmos, &lab, xy.data());
            // F90 mecih writes cimat into ab(1,1) onward; those residues are
            // visible to the next deri2 entry (D2EN) and to supdot below.
            for (int k = 1; k <= (lab * (lab + 1)) / 2; ++k) ab_buf[k] = cmoflat[k];
            double sum = 0.0;
            if (gse > 1e10) sum = -sum;  // dummy use of gse (F90 keeps it)
            // supdot writes its result into s[1..lab] (s[0] is padding), so
            // pass a 1-based-padding buffer and read wrk[kk] (Fortran work(kk)).
            std::vector<double> wrk((size_t)lab + 1, 0.0);
            for (int i = 1; i <= nstate; ++i) {
                int j = (i - 1) * lab + 1;
                supdot(wrk.data(), cmoflat.data(), &vectci[j - 1], lab);
                for (int kk = 1; kk <= lab; ++kk)
                    sum += vectci[j + kk - 1] * wrk[kk];
            }
            dxyzr[ivar - 1] += sum * cnst / nstate;
        }
        // Copy flat outputs back into the [col][row] containers.
        for (int iv = 1; iv <= nvar_nvo; ++iv) {
            for (int i = 1; i <= minear; ++i) f[iv][i] = fflat[(iv - 1) * minear + i];
            for (int i = 1; i <= ninear; ++i) fd[iv][i] = fdflat[(iv - 1) * ninear + i];
        }
        for (int jj = 1; jj <= ncolmax; ++jj)
            for (int i = 1; i <= ninear; ++i)
                fci[jj][i] = fciflat[(jj - 1) * ninear + i];
    }

done:
    return;
}
