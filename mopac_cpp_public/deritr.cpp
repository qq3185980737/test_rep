// deritr.cpp — C++ translation of MOPAC 2016 "deritr.F90" (complete).
// Finite-difference derivatives of the energy with respect to internal
// coordinates using full SCF calculations (only used when no other
// derivative calculation will do). errfn(i) = (E(x) - E(x-delta))*const/step
// (backward difference), or centered (E(x+d)-E(x-d))/2d with PRECISE.
#include "deritr.h"

#include <cmath>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "gmetry.h"
#include "hcore.h"
#include "iter.h"
#include "molkst_C.h"
#include "post_scf_corrections.h"
#include "reada.h"
#include "symtry.h"

// MOZYME SCF path (mozyme=false in tests; stubs elsewhere).
void hcore_for_MOZYME();
void iter_for_MOZYME(double& ee);

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;

namespace {
bool debug = false, precise = false;
int icalcn = 0;
double xderiv = 0.0, delta = 0.0, cnst = 0.0;
}  // namespace

void deritr() {
    if (icalcn != numcal) {
        debug = keywrd.find("DERITR") != std::string::npos;
        precise = keywrd.find("PRECISE") != std::string::npos;
        bool PM6_H = method_pm7 ||
                     ((keywrd.find("PM6-D") != std::string::npos) +
                          (keywrd.find("PM6-H") != std::string::npos) !=
                      0);
        icalcn = numcal;
        cnst = fpc_9;
        delta = 0.0;
        int i = (int)keywrd.find(" DELTA");
        if (i > 0) delta = reada(keywrd, i);
        if (precise) {
            if (delta < 1e-10) delta = 0.001;
            xderiv = 0.5 / delta;
        } else {
            if (delta < 1e-10) delta = 0.0002;
            xderiv = 1.0 / delta;
        }
        (void)PM6_H;
    }
    std::vector<double> xparam(3 * natoms + 1, 0.0);
    for (int i = 1; i <= nvar; ++i)
        xparam[i] = geo[loc[2][i]][loc[1][i]];
    double escf_store = escf;
    double enuclr_store = enuclr;
    double elect_store = elect;
    double aa = 0.0;
    if (!precise) {
        // Establish the energy at the current point.
        if (ndep != 0) symtry();
        gmetry(geo, coord);
        if (norbs * nelecs > 0) {
            if (mozyme) {
                hcore_for_MOZYME();
                if (moperr) return;
                aa = 0.0;
                iter_for_MOZYME(aa);
            } else {
                hcore();
                iter(aa, true, true);
            }
        } else {
            aa = 0.0;
        }
        if (method_pm7 ||
            ((keywrd.find("PM6-D") != std::string::npos) +
                 (keywrd.find("PM6-H") != std::string::npos) !=
             0)) {
            double sum = 0.0;
            post_scf_corrections(sum, false);
            aa = aa + sum / cnst;
        }
    }
    // Restore the density matrix.
    for (size_t ii = 0; ii < pa.size(); ++ii) p[ii] = pa[ii] * 2.0;
    aa = aa + enuclr;
    for (int i = 1; i <= nvar; ++i) {
        int k = loc[1][i];
        int l = loc[2][i];
        double xstore = xparam[i];
        for (int j = 1; j <= nvar; ++j)
            geo[loc[2][j]][loc[1][j]] = xparam[j];
        double ee = 0.0;
        if (precise) {
            geo[l][k] = xstore + delta;
            if (ndep != 0) symtry();
            gmetry(geo, coord);
            if (norbs * nelecs > 0) {
                if (mozyme) {
                    hcore_for_MOZYME();
                    if (moperr) return;
                    aa = 0.0;
                    iter_for_MOZYME(aa);
                } else {
                    hcore();
                    iter(aa, true, true);
                }
                if (method_pm7 ||
                    ((keywrd.find("PM6-D") != std::string::npos) +
                         (keywrd.find("PM6-H") != std::string::npos) !=
                     0)) {
                    double sum = 0.0;
                    post_scf_corrections(sum, false);
                    aa = aa + sum / cnst;
                }
            } else {
                aa = 0.0;
            }
            aa = aa + enuclr;
        }
        geo[l][k] = xstore - delta;
        if (ndep != 0) symtry();
        gmetry(geo, coord);
        if (norbs * nelecs > 0) {
            if (mozyme) {
                hcore_for_MOZYME();
                if (moperr) return;
                ee = 0.0;
                iter_for_MOZYME(ee);
            } else {
                hcore();
                iter(ee, true, true);
            }
            if (method_pm7 ||
                ((keywrd.find("PM6-D") != std::string::npos) +
                     (keywrd.find("PM6-H") != std::string::npos) !=
                 0)) {
                double sum = 0.0;
                post_scf_corrections(sum, false);
                ee = ee + sum / cnst;
            }
        } else {
            ee = 0.0;
        }
        ee = ee + enuclr;
        errfn[i] = (aa - ee) * cnst * xderiv;
    }
    escf = escf_store;
    enuclr = enuclr_store;
    elect = elect_store;
    (void)chanel_C::iw;
}
