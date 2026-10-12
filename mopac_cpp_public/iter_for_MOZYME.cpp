// iter_for_MOZYME.cpp — C++ translation of MOPAC 2016 "iter_for_MOZYME.F90".
//
// The MOZYME SCF driver: runs the big SCF loop over localized MOs (LMOs).
//   - initialization (keyword parsing: THRESH/RELTHR/PL/REORTH/TIMES/DENOUT/
//     ITRY/RE-LOC/1SCF/BIGSCF/OLDEN/LEWIS/DEBUG/DENS/FOCK/ITER)
//   - tidy() expansion loop (occupied + virtual sets, re-allocate 60% larger)
//   - density_for_MOZYME / setupk / buildf / helecz / check / eimp / diagg /
//     addhb (every 3 iterations) / cnvgz (three-point extrapolation) /
//     isitsc convergence test / level-shift logic / reorth / panic (SHUT).

#include "iter_for_MOZYME.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "funcon_C.h"
#include "common_arrays_C.h"
#include "iter_C.h"
#include "cosmo_C.h"
#include "linear_cosmo.h"
#include "tidy.h"
#include "pinout.h"
#include "density_for_MOZYME.h"
#include "makvec.h"
#include "local_for_MOZYME.h"
#include "scfcri.h"
#include "setupk.h"
#include "buildf.h"
#include "helecz.h"
#include "check.h"
#include "eimp.h"
#include "diagg.h"
#include "addhb.h"
#include "cnvgz.h"
#include "chrge_for_MOZYME.h"
#include "isitsc.h"
#include "reorth.h"
#include "geout.h"
#include "mopend.h"
#include "memory_error.h"
#include "timer.h"
#include "vecprt_for_MOZYME.h"
#include "reada.h"
#include "ijbo.h"
#include "to_screen.h"

using namespace MOZYME_C;
using namespace molkst_C;
using namespace chanel_C;
using namespace funcon_C;
using namespace common_arrays_C;
using namespace iter_C;
using namespace cosmo_C;
using namespace linear_cosmo;

namespace {

// Fortran INDEX(): 1-based position of substring, 0 if not found.
int Index(const std::string& s, const std::string& sub) {
    std::string::size_type p = s.find(sub);
    return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// Fortran NINT().
int Nint(double x) { return static_cast<int>(std::lround(x)); }

}  // namespace

void iter_for_MOZYME(double& ee) {
    static char xchar = ' ';
    static bool prtpls = false, orthog = false, times = false, scf1 = false,
                bigscf = false, debug = false, prtden = false, prtfok = false,
                okscf = false, opend = false, panic = true, store_useps = false;
    static int i = 0, idnout = 0, nocc1 = 0, nvir1 = 0, istabl = 0,
               idiagg = 0, nhb = 0, itrmax = 0, niter = 0, lno = 0,
               iemax = 0, iemin = 0, lnv = 0, mn = 0, indi = 0, j = 0, k = 0,
               l = 0, nij = 0, re_local = 0, icalcn = 0, imol = 0, nmol = 0,
               lstart = 0;
    static double selcon = 0.0, eold = 0.0, sum = 0.0;

    if (nmol != numcal) {
        //
        //  INITIALIZE
        //
        //     THRESH          ERROR IN H.O.F., RELATIVE TO THRESH=1.D-13
        //
        //     1.D-8             0.21??     SCF unstable
        //     1.D-9             0.0294     SCF stable
        //     1.D-10            0.0112     SCF stable
        //     1.D-11            0.0040     SCF stable
        //     1.D-12            0.0011     SCF stable
        //
        //   IF THRESH IS CHANGED, THEN ALSO CHANGE DEFAULT IN WRTKEY
        //
        thresh = 1.e-13;
        scfref = 0.0;
        i = Index(keywrd, " RELTHR");
        if (i != 0) thresh = reada(keywrd, i) * thresh;
        i = Index(keywrd, " THRESH");
        if (i != 0) thresh = reada(keywrd, i);
        thresh = std::max(thresh, 1.e-25);
        prtpls = (Index(keywrd, " PL ") != 0);
        orthog = (Index(keywrd, " REORTH") != 0);
        times = (Index(keywrd, " TIMES") != 0);
        i = Index(keywrd, " DENOUT=");
        if (i != 0) {
            idnout = Nint(reada(keywrd, i + 8));
        } else {
            idnout = 10000000;
        }
        itrmax = 2000;
        if (Index(keywrd, " ITRY") != 0) {
            itrmax = Nint(reada(keywrd, Index(keywrd, " ITRY")));
        }
        re_local = 100000000;
        i = Index(keywrd, " RE-LOC");
        if (i != 0) {
            j = Index(keywrd.substr(i + 7 - 1), " ") + i + 7;
            if (Index(keywrd.substr(i - 1, j - i + 1), "=") != 0) {
                // Allow for RE-LOC=, RE-LOCAL=, etc.
                i = Index(keywrd.substr(i - 1, j - i + 1), "=") + i;
                re_local = Nint(reada(keywrd, i));
            }
        }
        nocc1 = 0;
        nvir1 = 0;
        scf1 = (Index(keywrd, " 1SCF") != 0);
        bigscf = (Index(keywrd, " BIGSCF") != 0);
        if (Index(keywrd, " OLDEN") == 0) bigscf = true;
        debug = (Index(keywrd, " DEBUG") != 0);
        prtden = (Index(keywrd, " DENS") != 0 && debug);
        prtfok = (Index(keywrd, " FOCK") != 0 && debug);
        debug = (Index(keywrd, " ITER") != 0);
        idiagg = 0;
        nocc1 = noccupied;
        store_useps = useps;
        if (Index(keywrd, " LEWIS") != 0) {
            makvec();
            return;
        }
        if (Index(keywrd, " OLDEN") != 0) {
            pinout(0);
            if (moperr) return;
            density_for_MOZYME(p, 0, noccupied, partp);
            partp = p;
        } else {
            useps = false;
            if (Index(keywrd, "OLD_SCF") == 0) makvec();
        }
        if (moperr) return;
        //
        //  A NEW MOLECULE, THEREFORE SEARCH FOR ALL WEAK INTERACTIONS.
        //
        nhb = 0;
    } else {
        //
        //  A MODIFIED GEOMETRY, THEREFORE SEARCH ONLY FOR NEW
        //  OR VERY WEAK INTERACTIONS.
        //
        nhb = 3;
    }

    if ((nscf + 1) % re_local == 0) {
        std::printf("\n%10s\n\n", "  LMOs being Re-Localized");
        local_for_MOZYME("OCCUPIED");
        local_for_MOZYME("VIRTUAL");
    }
    useps = store_useps;
    if (lpka) useps = true;
    if (icalcn != step_num) {
        istabl = 0;
        eold = 0.0;
        ovmax = 0.0;
        scfcri(selcon);
        idiagg = 0;
        tiny = 0.0;
    }
    //
    // Fill the array IDIAG pointing to the diagonal elements of P
    //
    l = 0;
    for (i = 1; i <= numat; ++i) {
        j = ijbo(i, i);
        for (k = 1; k <= iorbs[i]; ++k) {
            l = l + 1;
            j = j + k;
            idiag[l] = j;
        }
    }
    //
    //  Zero out POLD and P1 so that the new calculation is not affected
    //  by a previous one.
    //
    std::fill(pold.begin(), pold.end(), 0.0);
    std::fill(p1.begin(), p1.end(), 0.0);
    nscf = nscf + 1;
    //
    //  Force IDIAGG to be even - this is to ensure that DIAGG1
    //  re-builds the interaction list, in case the temporary space
    //  has been used between the calls to ITER.
    //
    if (idiagg % 2 == 1) idiagg = idiagg + 1;
    iemin = 0;
    iemax = 0;
    niter = 0;
    shift = 0.0;
    use_three_point_extrap = true;
    if (times) timer(" At start of ITER");

    //***********************************************************
    //
    //   Everything is now set up to allow the SCF to be run
    //
    //***********************************************************
    for (;;) {  // Big loop to run the SCF
        //----------------
        for (;;) {
            tidy(noccupied, ncf.data(), icocc.data(), icocc_dim, cocc.data(),
                 cocc_dim, nncf.data(), ncocc.data(), lno, mn, 1);
            if (moperr) {
                //
                //  During a run of "tidy", the amount of expansion space for
                //  the LMO's to use had become small.  The LMOs were stored to
                //  disc.  The old arrays will now be deleted and re-created
                //  60% larger than before.  Then the LMOs are read off disc.
                //
                icocc.clear();
                cocc.clear();
                icvir.clear();
                cvir.clear();
                //
                //  Re-allocate more memory
                //
                icocc_dim = Nint(icocc_dim * 1.6);
                cocc_dim = Nint(cocc_dim * 1.6);
                icvir_dim = Nint(icvir_dim * 1.6);
                cvir_dim = Nint(cvir_dim * 1.6);
                try {
                    icocc.resize(icocc_dim + 1);
                    cocc.resize(cocc_dim + 1);
                    icvir.resize(icvir_dim + 1);
                    cvir.resize(cvir_dim + 1);
                } catch (...) {
                    memory_error(" iter_for MOZYME");
                    return;
                }
                //
                //  Read in old density
                //
                pinout(0);
                moperr = false;
            } else {
                break;
            }
        }
        for (;;) {
            tidy(nvirtual, nce.data(), icvir.data(), icvir_dim, cvir.data(),
                 cvir_dim, nnce.data(), ncvir.data(), lnv, mn, 2);
            if (moperr) {
                //  Delete old memory
                icocc.clear();
                cocc.clear();
                icvir.clear();
                cvir.clear();
                //
                //  Re-allocate more memory
                //
                icocc_dim = Nint(icocc_dim * 1.6);
                cocc_dim = Nint(cocc_dim * 1.6);
                icvir_dim = Nint(icvir_dim * 1.6);
                cvir_dim = Nint(cvir_dim * 1.6);
                try {
                    icocc.resize(icocc_dim + 1);
                    cocc.resize(cocc_dim + 1);
                    icvir.resize(icvir_dim + 1);
                    cvir.resize(cvir_dim + 1);
                } catch (...) {
                    memory_error(" iter_for MOZYME");
                    return;
                }
                //
                //  Read in old density
                //
                pinout(0);
                moperr = false;
            } else {
                break;
            }
        }
        //----------------
        nocc1 = nelred / 2;
        nvir1 = norred - nocc1;
        //
        //   REMOVE ELECTRON DENSITY DUE TO LMO'S INVOLVED IN SCF FROM
        //   THE DENSITY MATRIX
        //
        if (imol != numcal || (icalcn != step_num && numat > numred + 1)) {
            //---------------------------------------------------------
            //
            //   THIS PART IS ONLY RUN WHEN ICALCN IS INCREMENTED
            //
            density_for_MOZYME(p, 0, noccupied, partp);  // Build the whole density matrix
            if (times) timer(" After DENSIT");
            if (prtden) {
                std::printf(" DENSITY MATRIX TO GO INTO PARTP\n");
                vecprt_for_MOZYME(p.data(), norbs);
            }
            setupk(nocc1);  // Work out the atom list to be used in the SCF
            if (times) timer(" After SETUPK");
            if (imol == numcal && numat > numred + 1) {
                // Remove density due to atoms to be used in the SCF
                density_for_MOZYME(partp, -1, nocc1, p);
                if (prtden) {
                    std::printf(" DENSITY MATRIX IN PARTP\n");
                    vecprt_for_MOZYME(partp.data(), norbs);
                }
            }
            buildf(f, partf, 0);
            if (prtfok) {
                std::printf(" FOCK MATRIX AT START OF ITER\n");
                vecprt_for_MOZYME(f.data(), norbs);
            }
            if (icalcn != step_num) {
                if (times) timer(" After BUILDF");
                ee = helecz();
                if (times) timer(" After HELEC");
                escf = (ee + enuclr) * fpc_9 + atheat;
                if (useps) escf = escf + solv_energy * fpc_9;
                if (prtpls)
                    std::printf("\n%10s%16.6f%3s\n",
                                " PLS ESCF USING THE OLD LMOs:", escf,
                                " KCAL/MOL");
            }
            if (imol == numcal && numat > numred + 1) buildf(partf, f, -1);
            icalcn = step_num;
            imol = numcal;
        }
        //
        //  Correct any small errors in normalization
        //
        check(nocc1, nncf, ncf, icocc, iorbs, ncocc, cocc);
        check(nvir1, nnce, nce, icvir, iorbs, ncvir, cvir);
        if (moperr) return;
        if ((niter + 1) % idnout == 0) {
            std::printf(" .den FILE TO BE WRITTEN OUT\n");
            pinout(1);
            std::printf(" .den FILE WRITTEN OUT\n");
        }
        eimp();
        if (nmol == numcal && numat > numred + 1) {
            indi = 1;
        } else {
            indi = 0;
        }

        if (bigscf || numcal != 1) {
            diagg(f, nocc1, nvir1, idiagg, partp, indi);
            idiagg = idiagg + 1;
        } else {
            density_for_MOZYME(p, mode, nocc1, partp);
            bigscf = true;
        }
        niter = niter + 1;
        if (niter % 3 == 0 && nhb < 4) {
            nhb = nhb + 1;
            //
            //   Check for missed hydrogen bonds and other unusual bonds
            //
            addhb(nocc1, nvir1, idiagg, nij, nhb);
            //
            //   If hydrogen bonds have been made, set IDIAGG even for DIAGG
            //   to make new interactions.
            //
            if (nij != 0 && idiagg % 2 == 1) {
                idiagg = idiagg + 1;
            }
        }
        if (times) {
            if (i == 1) {
                timer(" After DENS+1");
            } else {
                timer(" After DENSIT");
            }
        }
        if (use_three_point_extrap)
            cnvgz(p, pold, p1, p2, p3, niter, idiag);
        if (times) timer(" After CNVG");
        if (prtden) {
            std::printf(" DENSITY MATRIX ON ITERATION%4d\n", niter);
            vecprt_for_MOZYME(p.data(), norbs);
        }
        if (nmol == numcal && numat > numred + 1) {
            buildf(f, partf, 1);
        } else {
            buildf(f, partf, 0);
        }
        if (itrmax < 3) return;
        if (times) timer(" After BUILDF");
        if (prtfok) {
            std::printf(" FOCK    MATRIX ON ITERATION%4d\n", niter);
            vecprt_for_MOZYME(f.data(), norbs);
        }
        ee = helecz();
        escf = (ee + enuclr) * fpc_9 + atheat;
        if (times) timer(" After HELEC");
        if (useps) escf = escf + solv_energy * fpc_9;
        energy_diff = escf - eold;
        eold = escf;
        if (std::abs(ovmax) < 5.0 * selcon) {
            c_proc = 1.0;
        } else {
            c_proc = 5.0 * selcon / std::abs(ovmax);
        }
        if (prtpls || (debug && niter > itrmax - 20)) {
            escf = std::max(-999999.0, std::min(999999.0, escf));
            if (std::abs(energy_diff) > 9999.0) energy_diff = 0.0;
            char linebuf[160];
            std::snprintf(linebuf, sizeof(linebuf),
                          "ITER.%7d PLS=%10.3e%10s ENERGY %13.5f DELTAE%13.7f",
                          niter, ovmax, "", escf, energy_diff);
            line = linebuf;
            std::printf("%s\n", line.c_str());
            to_screen(line);
            if (debug) {
                std::printf("TINY:%9.6f SUMT:%7.1f OVMAX:%9.6f SUMB:%8.2f "
                            "DIFF:%11.3f IJ:%7d\n",
                            tiny, sumt, ovmax, sumb, energy_diff, ijc);
            }
        }
        if (debug) {
            chrge_for_MOZYME(p.data(), ws.data(), numat, iorbs.data());
            sum = 0.0;
            for (i = 1; i <= numat; ++i) sum = sum + std::abs(ws[i]);
            std::printf(" Atomic Electron Population on Iteration:%4d\n",
                        niter);
            for (i = 1; i <= numat; ++i) {
                std::printf("%8.4f", ws[i]);
                if (i % 10 == 0) std::printf("\n");
            }
            std::printf("\n Variance:%12.6f\n", sum);
        }
        isitsc(escf, selcon, emin, iemin, iemax, okscf, niter, itrmax);
        if (!bigscf && numcal == 1) {
            break;
        } else if (okscf && niter > 1 && (emin != 0.0 || niter > 3)) {
            break;
        }
        if (use_three_point_extrap) {
            if (niter % 3 == 2 && std::abs(energy_diff) < 0.1) {
                use_three_point_extrap = false;
                lstart = niter;
                shift = 0.0;
            }
        } else {
            if (energy_diff > 0.0 && shift < 11.0 && niter > lstart + 2) {
                shift = shift + 2.0;
                lstart = niter;
            }
        }
    }
    //************************************************************
    //
    //   The SCF equations are now solved
    //
    //************************************************************
    if (istabl > 100000 && escf - emin > 200.0 && panic) {
        std::printf(" Something disastrous has happened.  "
                    "The LMOs are probably corrupt.\n");
        std::printf(" The job should be restarted using RESTART but"
                    " OLDENS should NOT be used.\n");
        std::printf("\n%10s%16.6f\n%10s%16.6f\n",
                    " Value of previous Heat of Formation:", emin,
                    " Value of current Heat of Formation: ", escf);
        //
        //          Start the SHUT command.
        //
        {
            std::FILE* fp = std::fopen(end_fn.c_str(), "a");
            if (fp) {
                std::fputc(xchar, fp);
                std::fclose(fp);
            } else {
                //  The SHUT command is faulty.  Stop everything.
                geout(iw);
                mopend("Severe fault in RESTART");
            }
        }
        panic = false;
    }
    if (std::abs(escf - emin) < 1.0) {
        istabl = istabl + 1;
    } else {
        istabl = 0;
    }
    if (escf < emin || emin == 0.0) {
        emin = escf;
    }
    if (!scf1) {
        if (numat > numred + 1) {
            density_for_MOZYME(p, 1, nocc1, partp);
        } else {
            density_for_MOZYME(p, 0, nocc1, partp);
        }
    }
    icalcn = step_num;
    imol = numcal;
    if (orthog && (nscf + 1) % 10 == 1) {
        reorth(ws.data());  //   Re-orthogonalize the LMO's
        density_for_MOZYME(p, 0, noccupied, partp);
        buildf(f, partf, 0);
        ee = helecz();
        escf = (ee + enuclr) * fpc_9 + atheat;
        if (useps) escf = escf + solv_energy * fpc_9;
        std::printf("\n%10s%16.6f%3s\n",
                    " ESCF AFTER RE-ORTHOGONALIZING LMOs", escf,
                    " KCAL/MOL");
    }
    nmol = numcal;
    return;
}
