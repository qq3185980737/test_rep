// thermo.cpp — C++ translation of MOPAC 2016 "thermo.F90".
// Calculated thermodynamic properties (partition functions, HOF, H, Cp, S).
#include "thermo.h"
#include "molkst_C.h"
#include "funcon_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace molkst_C {
extern std::string keywrd, title, koment;
extern double mol_weight;
extern double temp_1, temp_2;
extern int ilim;
}
namespace common_arrays_C {
extern std::vector<double> T_range, HOF_tot, H_tot, Cp_tot, S_tot;
}
namespace chanel_C {
extern int iw;
}

extern double reada(const std::string&, int);
extern void web_message(int, const char*);

void thermo(double a, double b, double c, int linear, double sym,
            double* vibs, int nvibs, double escf) {
    double R = funcon_C::fpc_5, h = funcon_C::fpc_6, ak = funcon_C::fpc_7,
           ac = funcon_C::fpc_8;
    int it1 = 200, it2 = 400, istep = 10;
    std::string tmpkey = molkst_C::keywrd;
    std::string trimmed = tmpkey;
    // Trim trailing whitespace for the index search.
    while (!trimmed.empty() && (trimmed.back() == ' ' || trimmed.back() == '\t')) trimmed.pop_back();
    tmpkey = trimmed;
    int i = (int)tmpkey.find("THERMO(") + (int)tmpkey.find("THERMO=(");
    if (i != 0) {
        if (i == (int)std::string::npos) i = (int)tmpkey.find("THERMO(");
        if (i == (int)std::string::npos) i = (int)tmpkey.find("THERMO=(");
        // Erase all text except THERMO data.
        for (int k = 0; k < i; ++k) tmpkey[k] = ' ';
        int jp = (int)tmpkey.find(')');
        if (jp != (int)std::string::npos)
            for (int k = jp; k < (int)tmpkey.size(); ++k) tmpkey[k] = ' ';
        it1 = (int)std::lround(reada(tmpkey, i));
        if (it1 < 100) {
            std::fprintf(stdout, "\n\n          TEMPERATURE RANGE STARTS TOO LOW, LOWER BOUND IS RESET TO 30K\n");
            it1 = 100;
        }
        i = (int)tmpkey.find(',');
        if (i != (int)std::string::npos) {
            tmpkey[i] = ' ';
            it2 = (int)std::lround(reada(tmpkey, i));
            if (it2 < it1) {
                it2 = it1 + 200;
                istep = 10;
                goto thermo_loop;
            }
            i = (int)tmpkey.find(',');
            if (i != (int)std::string::npos) {
                tmpkey[i] = ' ';
                istep = (int)std::lround(reada(tmpkey, i));
                istep = std::max(1, istep);
            } else {
                istep = (it2 - it1) / 20;
                if (istep == 0) istep = 1;
                if (istep >= 2 && istep < 5) istep = 2;
                if (istep >= 5 && istep < 10) istep = 5;
                if (istep >= 10 && istep < 20) istep = 10;
                if (istep >= 20 && istep < 50) istep = 20;
                if (istep >= 50 && istep < 100) istep = 50;
                istep = std::min(100, istep);
            }
        } else {
            it2 = it1 + 200;
        }
    }
thermo_loop:
    std::fprintf(stdout, "\n\n%s\n", molkst_C::title.c_str());
    std::fprintf(stdout, "%s\n", molkst_C::koment.c_str());
    if (linear)
        std::fprintf(stdout, "\n          MOLECULE IS LINEAR\n");
    else
        std::fprintf(stdout, "\n          MOLECULE IS NOT LINEAR\n");
    if (nvibs > 99)
        std::fprintf(stdout, "\n          THERE ARE %5d GENUINE VIBRATIONS IN THIS SYSTEM\n", nvibs);
    else
        std::fprintf(stdout, "\n          THERE ARE %3d GENUINE VIBRATIONS IN THIS SYSTEM\n", nvibs);
    std::fprintf(stdout, "          THIS THERMODYNAMICS CALCULATION IS LIMITED TO\n");
    std::fprintf(stdout, "          MOLECULES WHICH HAVE NO INTERNAL ROTATIONS\n\n");
    std::fprintf(stdout, "\n                    CALCULATED THERMODYNAMIC PROPERTIES\n");
    std::fprintf(stdout, "                                          *\n");
    std::fprintf(stdout, "   TEMP. (K)   PARTITION FUNCTION   H.O.F.    ENTHALPY  HEAT CAPACITY ENTROPY\n");
    std::fprintf(stdout, "                                   KCAL/MOL   CAL/MOLE    CAL/K/MOL  CAL/K/MOL\n");
    for (int k = 1; k <= nvibs; ++k) vibs[k] = std::fabs(vibs[k]);
    molkst_C::ilim = 1;
    for (int itemp = it1; itemp <= it2; itemp += istep) {
        molkst_C::ilim += 1;
        common_arrays_C::T_range[molkst_C::ilim] = itemp;
    }
    common_arrays_C::T_range[1] = 298.0;
    double qtr2 = 2.0 * funcon_C::pi * ac * molkst_C::mol_weight / funcon_C::fpc_10;
    if (molkst_C::ilim == 2 &&
        std::fabs(common_arrays_C::T_range[1] - common_arrays_C::T_range[2]) < 1.e-4)
        molkst_C::ilim = 1;
    double cptot = 0.0, stot = 0.0;
    for (int ir = 1; ir <= molkst_C::ilim; ++ir) {
        int itemp = (int)common_arrays_C::T_range[ir];
        double T = itemp;
        double c1 = h * ac / ak / T;
        double qv = 1.0, hv = 0.0, cpv = 0.0, sv1 = 0.0;
        for (int k = 1; k <= nvibs; ++k) {
            double wi = vibs[k];
            if (wi < 1.e-4) continue;
            double ewj = std::exp(-wi * c1);
            double ewjr = 1.0 - ewj;
            qv = qv / ewjr;
            double e0 = wi * ewj / ewjr;
            hv = hv + e0;
            cpv = cpv + e0 * e0 / ewj;
            sv1 = sv1 + std::log(ewjr);
        }
        double e0 = R * c1;
        double sv = hv * e0 - R * sv1;
        hv = hv * e0 * T;
        cpv = cpv * e0 * c1;
        double qr = 0.0, cpr = 0.0, hr = 0.0, sr = 0.0;
        if (nvibs == 0) {
            qr = 0.0; cpr = 0.0; hr = cpr * T; sr = 0.0;
        } else if (!linear) {
            e0 = funcon_C::pi / (a * b * c);
            qr = std::sqrt(e0 / c1) / c1 / sym;
            cpr = 1.5 * R;
            hr = cpr * T;
            sr = R * (1.5 * std::log(1.0 / c1) - std::log(sym) + std::log(e0) / 2.0 + 1.5);
        } else {
            qr = 1.0 / (c1 * a * sym);
            cpr = R;
            hr = cpr * T;
            sr = R * std::log(qr) + R;
        }
        double qint = qv * qr;
        double hint = hv + hr;
        double cpint = cpv + cpr;
        double sint = sv + sr;
        double qtr = std::sqrt(qtr2 / c1 / h);
        qtr = qtr * qtr * qtr;
        double cptr = 2.5 * R;
        double htr = cptr * T;
        double str = 4.96804 * (std::log(T) + 0.6 * std::log(molkst_C::mol_weight)) - 2.31482;
        cptot = cptr + cpint;
        stot = str + sint;
        double htot = htr + hint;
        double h298 = htot;
        if (ir == 1) h298 = htot;
        std::fprintf(stdout, "\n%7.2f  VIB.%15.4e          %17.4f %11.4f %11.4f\n", T, qv, hv, cpv, sv);
        std::fprintf(stdout, "        ROT.%15.4e          %17.4f %11.4f %11.4f\n", qr, hr, cpr, sr);
        std::fprintf(stdout, "        INT.%15.4e          %17.4f %11.4f %11.4f\n", qint, hint, cpint, sint);
        std::fprintf(stdout, "        TRA.%15.4e          %17.4f %11.4f %11.4f\n", qtr, htr, cptr, str);
        std::fprintf(stdout, "        TOT.                  %17.3f %13.4f %11.4f %11.4f\n",
                     escf + (htot - h298) / 1000.0, htot, cptot, stot);
        common_arrays_C::HOF_tot[ir] = escf + (htot - h298) / 1000.0;
        common_arrays_C::H_tot[ir] = htot;
        common_arrays_C::Cp_tot[ir] = cptot;
        common_arrays_C::S_tot[ir] = stot;
    }
    if (molkst_C::ilim == 1 && std::fabs(common_arrays_C::T_range[1] - 298) < 1.e1) {
        molkst_C::temp_1 = cptot;
        molkst_C::temp_2 = stot;
    }
    std::fprintf(stdout, "\n  *: NOTE: HEATS OF FORMATION ARE RELATIVE TO THE\n");
    std::fprintf(stdout, "            ELEMENTS IN THEIR STANDARD STATE AT 298K\n");
    std::fprintf(stdout, "            (=  Standard Enthalpy of Formation).\n");
    std::fprintf(stdout, "\n            Hvib:  Zero-point energy is not included.\n");
    std::fprintf(stdout, "                   frequencies of less than zero cm-1 are not included.\n");
    std::fprintf(stdout, "            Hrot = (3/2)RT\n");
    std::fprintf(stdout, "            Htra = (3/2)RT + pV = (5/2)RT\n");
    std::fprintf(stdout, "\n            Heat capacity is Cp, not Cv\n");
    web_message(0, "thermochemistry.html");
}
