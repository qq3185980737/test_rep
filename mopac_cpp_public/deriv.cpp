// deriv.cpp — C++ translation of MOPAC 2016 "deriv.F90" (complete).
// Computes the derivatives of the energy with respect to the internal
// coordinates. Main arrays: loc (internal-coordinate address table), geo
// (internal coordinates), gradnt (derivatives on exit).
// Path: Cartesian derivatives from dernvo (analytical CI / half-electron)
// or dcart (variationally optimized) -> optional field / pressure terms ->
// jcarin builds the Jacobian of the coordinate transformation -> mxm maps
// Cartesian derivatives to internal-coordinate gradients.
#include "deriv.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "dcart.h"
#include "derivs_C.h"
#include "dernvo.h"
#include "deritr.h"
#include "dfield.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "gmetry.h"
#include "jcarin.h"
#include "molkst_C.h"
#include "mopend.h"
#include "mxm.h"
#include "post_scf_corrections.h"
#include "symmetry_C.h"
#include "symtry.h"
#include "upcase.h"
#include "volume.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace derivs_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace symmetry_C;

namespace {
bool aifrst = true, saddle = false, debug = false, field = false,
     large = false, precis = false, DH_correction = false;
bool scf1 = false, halfe = false, slow = false, intn = true,
     geochk = false, ci = false, aic = false, noanci = false;
int icalcn = 0, idelta = 0, nw2 = 0;
double grlim = 0.0;
double change[3] = {0.0, 0.0, 0.0};
}  // namespace

static double dot3(const double* a, const double* b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

void deriv(const std::vector<std::vector<double>>& geo,
           std::vector<double>& gradnt) {
    double sum = 0.0;
    if (icalcn != numcal) {
        aifrst = keywrd.find(" RESTART") == std::string::npos;
        saddle = keywrd.find(" SADDLE") != std::string::npos;
        debug = keywrd.find(" DERIV") != std::string::npos;
        field = keywrd.find(" FIELD") != std::string::npos;
        large = keywrd.find("LARGE") != std::string::npos;
        precis = keywrd.find(" PREC") != std::string::npos;
        DH_correction =
            ((keywrd.find(" PM6-D") != std::string::npos) +
             (keywrd.find(" PM6-H") != std::string::npos) != 0) ||
            method_pm7;
        intn = keywrd.find("  XYZ") == std::string::npos;
        if (saddle)
            nw2 = std::max(mpack, 9 * natoms * natoms);
        else
            nw2 = std::max(mpack, 6 * natoms * l123);
        errfn.assign(nvar + 1, 0.0);
        aidref.assign(nvar + 1, 0.0);
        work2.assign(nw2 + 1, 0.0);

        geochk = ((keywrd.find(" TS") != std::string::npos) +
                  (keywrd.find(" NLLSQ") != std::string::npos) +
                  (keywrd.find(" SIGMA") != std::string::npos)) != 0;
        geochk = geochk && intn && nvar >= numat * 3 - 6 && id == 0 &&
                 keywrd.find("GEO-OK") == std::string::npos;
        geochk = geochk && keywrd.find(" XYZ") == std::string::npos;
        ci = keywrd.find(" C.I.") != std::string::npos;
        scf1 = keywrd.find(" 1SCF") != std::string::npos;
        aic = keywrd.find("AIDER") != std::string::npos;
        if (aic && aifrst) {
            // Read ab-initio derivatives from the input job file.
            std::ifstream fin(job_fn, std::ios::in);
            std::string line;
            bool found = false;
            while (std::getline(fin, line)) {
                upcase(line, 80);
                if (line.find("AIDER") != std::string::npos) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                mopend("KEYWORD \"AIDER\" SPECIFIED, BUT NOT PRESENT AFTER "
                       "Z-MATRIX.  JOB STOPPED");
                return;
            }
            isok = false;
            const int j = natoms > 2 ? 3 * natoms - 6 : 1;
            int pos_guard = 0;
            while (std::getline(fin, line)) {
                std::istringstream iss(line);
                int got = 0;
                for (int i = 1; i <= j; ++i) {
                    double v;
                    if (iss >> v) {
                        aidref[i] = v;
                        ++got;
                    }
                }
                if (got == j) break;
                if (++pos_guard > 10000) {
                    mopend("FAULT IN READ OF AB INITIO DERIVATIVES");
                    return;
                }
            }
            for (int i = 1; i <= nvar; ++i) {
                int jc;
                if (loc[1][i] > 3)
                    jc = 3 * loc[1][i] + loc[2][i] - 9;
                else if (loc[1][i] == 3)
                    jc = loc[2][i] + 1;
                else
                    jc = 1;
                aidref[i] = aidref[jc];
            }
            if (ndep != 0) {
                for (int i = 1; i <= nvar; ++i) {
                    double s = aidref[i];
                    int cnt = 0;
                    for (int k = 1; k <= ndep; ++k)
                        if (loc[1][i] == locpar[k] &&
                            (loc[2][i] == idepfn[k] ||
                             (loc[2][i] == 3 && idepfn[k] == 14)))
                            ++cnt;
                    aidref[i] = aidref[i] + cnt * s;
                }
            }
        }
        grlim = precis ? 0.0001 : 0.01;
        halfe = (nopen > nclose && std::fabs(fract - 2.0) > 1e-20 &&
                 std::fabs(fract) > 1e-20) ||
                ci;
        idelta = -7;
        change[0] = std::pow(10.0, idelta);
        change[1] = std::pow(10.0, idelta);
        change[2] = std::pow(10.0, idelta);
    }
    if (nvar == 0) return;

    double gnorm = 0.0;
    std::vector<double> gold(nvar + 1, 0.0), xparam(nvar + 1, 0.0);
    for (int i = 1; i <= nvar; ++i) {
        gold[i] = gradnt[i];
        xparam[i] = geo[loc[2][i]][loc[1][i]];
        gnorm += gradnt[i] * gradnt[i];
    }
    gnorm = std::sqrt(gnorm);
    slow = false;
    noanci = false;
    if (halfe) {
        noanci = keywrd.find("NOANCI") != std::string::npos || nopen == norbs;
        slow = noanci && (gnorm < grlim || scf1);
    } else {
        slow = keywrd.find("NOANCI") != std::string::npos;
    }
    if (ndep != 0) symtry();
    std::vector<std::vector<double>> coord(4,
                                           std::vector<double>(natoms + 1, 0.0));
    gmetry(const_cast<std::vector<std::vector<double>>&>(geo), coord);
    // COORD now holds the Cartesian coordinates.

    if (halfe && !noanci && numat > 1) {
        dernvo();
        if (moperr) return;
    } else {
        // Variationally optimized derivatives.
        std::vector<std::vector<double>> dxyz2d(
            4, std::vector<double>(numat + 1, 0.0));
        dcart(coord, dxyz2d);
        for (int k = 1; k <= 3; ++k)
            for (int i = 1; i <= numat; ++i) dxyz[3 * (i - 1) + k] = dxyz2d[k][i];
    }
    if (DH_correction) post_scf_corrections(sum, true);
    if (field) dfield();

    if (std::fabs(pressure) > 1e-4) {
        if (id == 1) {
            int i = 3 * l123 / 2;
            double press =
                pressure / std::sqrt(dot3(&tvec[1][1], &tvec[1][1], 3));
            for (int j = 1; j <= 3; ++j) {
                dxyz[j + i] += tvec[j][1] * press;
                dxyz[j + i + 3] -= tvec[j][1] * press;
            }
        } else if (id == 3) {
            double tvec3[9];
            for (int c = 0; c < 3; ++c)
                for (int r = 1; r <= 3; ++r)
                    tvec3[c * 3 + (r - 1)] = tvec[r][c + 1];
            double summ = volume(tvec3, 3);
            double press1 = summ / dot3(&tvec[1][1], &tvec[1][1], 3) * pressure;
            double press2 = summ / dot3(&tvec[1][2], &tvec[1][2], 3) * pressure;
            double press3 = summ / dot3(&tvec[1][3], &tvec[1][3], 3) * pressure;
            int i = 3 * (l1u * (2 * l2u + 1) * (2 * l3u + 1) +
                         l2u * (2 * l3u + 1) + l3u - 1);
            for (int j = 1; j <= 3; ++j) {
                dxyz[j + i] += tvec[j][1] * press1;
                dxyz[j + i] += tvec[j][2] * press2;
                dxyz[j + i] += tvec[j][3] * press3;
            }
            i = 3 * (l1u * (2 * l2u + 1) * (2 * l3u + 1) +
                     l2u * (2 * l3u + 1) + l3u);
            for (int j = 1; j <= 3; ++j) dxyz[j + i] -= tvec[j][3] * press3;
            i = 3 * (l1u * (2 * l2u + 1) * (2 * l3u + 1) +
                     (l2u + 1) * (2 * l3u + 1) + l3u - 1);
            for (int j = 1; j <= 3; ++j) dxyz[j + i] -= tvec[j][2] * press2;
            i = 3 * ((l1u + 1) * (2 * l2u + 1) * (2 * l3u + 1) +
                     l2u * (2 * l3u + 1) + l3u - 1);
            for (int j = 1; j <= 3; ++j) dxyz[j + i] -= tvec[j][1] * press1;
        }
    }

    double step = change[0];
    int nstep = nw2 / (3 * numat * l123);
    for (int i = 1; i <= nvar; i += nstep) {
        int j = std::min(i + nstep - 1, nvar);
        int ncol = 0;
        jcarin(xparam.data(), step, precis, work2.data(), ncol, i, j);
        mxm(work2.data(), j - i + 1, &dxyz[1], ncol, &gradnt[i], 1);
        std::fprintf(stderr, "[DGD] i=%d j=%d nvar=%d ncol=%d nstep=%d step=%.9f\n", i, j, nvar, ncol, nstep, step);
        for (int r = 1; r <= (j - i + 1); ++r) {
            std::fprintf(stderr, "[DGD-R%d]", r);
            for (int cc = 1; cc <= ncol; ++cc)
                std::fprintf(stderr, " %+.6f", work2[(cc - 1) * (j - i + 1) + (r - 1)]);
            std::fprintf(stderr, "\n");
        }
    }
    step = precis ? 0.5 / step : 1.0 / step;
    for (int i = 1; i <= nvar; ++i) gradnt[i] *= step;
    std::fprintf(stderr, "[DGD-GD] step2=%.9f\n", step);
    for (int i = 1; i <= nvar; ++i)
        std::fprintf(stderr, "[DGD-GD%d] %+.9f\n", i, gradnt[i]);
    std::fprintf(stderr, "[DGD-DXYZ]");
    for (int i = 1; i <= 3 * numat; ++i)
        std::fprintf(stderr, " %+.9f", dxyz[i]);
    std::fprintf(stderr, "\n");
        if (geochk) {
            sum = dot3(&gradnt[1], &gradnt[1], nvar);
            if (sum < 2.0 &&
                dot3(&dxyz[1], &dxyz[1], 3 * numat) > std::max(4.0, sum * 4.0)) {
                for (int i = 1; i <= nvar; ++i) {
                    int j = (int)(xparam[i] / 3.141);
                    if (!(loc[2][i] == 2 && loc[1][i] > 3 &&
                          std::fabs(xparam[i] - j * pi) < 0.005))
                        continue;
                    mopend(" INTERNAL COORDINATE DERIVATIVES DO NOT REFLECT "
                           "CARTESIAN COORDINATE DERIVATIVES");
                    return;
                }
            }
        }
    if (slow) {
        deritr();
        icalcn = numcal;
        for (int i = 1; i <= nvar; ++i) errfn[i] = errfn[i] - gradnt[i];
    }
    cosine = dot3(&gradnt[1], &gold[1], nvar) /
             std::sqrt(dot3(&gradnt[1], &gradnt[1], nvar) *
                           dot3(&gold[1], &gold[1], nvar) +
                       1e-20);
    if (slow)
        for (int i = 1; i <= nvar; ++i) gradnt[i] = gradnt[i] + errfn[i];
    if (aic) {
        if (aifrst) {
            aifrst = false;
            for (int i = 1; i <= nvar; ++i)
                aicorr[i] = (-aidref[i]) - gradnt[i];
        }
        for (int i = 1; i <= nvar; ++i) gradnt[i] = gradnt[i] + aicorr[i];
    }
    (void)chanel_C::iw;
}
