// calpar.cpp — C++ translation of MOPAC 2016 "calpar.F90".
// sp_two_electron / inid / create_parameters_for_PMx_C are external stubs.

#include "calpar.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace parameters_C;

namespace {
void sp_two_electron() {}
void create_parameters_for_PMx_C(const std::string&, char) {}
}
// inid(): external linkage - real implementation in mndod.cpp (mndod.F90 inid).
void inid();

// nspqn(i): principal quantum number. 2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0.
static const int nspqn[108] = {
    0,
    /*1-2*/1,1,
    /*3-10*/2,2,2,2,2,2,2,2,
    /*11-18*/3,3,3,3,3,3,3,3,
    /*19-36*/4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    /*37-54*/5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
    /*55-86*/6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    /*87-107*/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

void calpar() {
    double p = 2.0;
    double p4 = p * p * p * p;
    sp_two_electron();

    for (int i = 1; i <= 107; ++i) am[i] = 0.0;

    for (int i = 2; i <= 97; ++i) {
        double gssc = std::max(ios[i] - 1, 0);
        int k = iop[i];
        double gspc = (double)(ios[i] * k);
        int l = std::min(k, 6 - k);
        double gp2c = (k * (k - 1)) / 2.0 + 0.5 * ((l * (l - 1)) / 2.0);
        double gppc = -0.5 * ((l * (l - 1)) / 2.0);
        double hspc = -k * ios[i] * 0.5;

        if (zp[i] < 1e-4 && zs[i] < 1e-4) continue;

        if (zp[i] < 0.3) zp[i] = 0.3;
        double hpp = 0.5 * (gpp[i] - gp2[i]);
        if (hpp < 0.1) hpp = 0.1;
        if (hsp[i] < 1e-7) hsp[i] = 1e-7;

        eisol[i] = uss[i] * ios[i] + upp[i] * iop[i] + udd[i] * iod[i] +
                   gss[i] * gssc + gpp[i] * gppc + gsp[i] * gspc +
                   gp2[i] * gp2c + hsp[i] * hspc;

        double qn = nspqn[i];
        dd[i] = (2 * qn + 1) * pow(4 * zs[i] * zp[i], qn + 0.5) /
                    pow(zs[i] + zp[i], 2 * qn + 2) / sqrt(3.0);
        qq[i] = sqrt((4 * qn * qn + 6 * qn + 2) / 20.0) / zp[i];

        const int jmax = 5;
        double gdd1 = pow(hsp[i] / (funcon_C::ev * dd[i] * dd[i]), 1.0 / 3.0);
        double d1 = gdd1;
        double d2 = gdd1 + 0.04;
        for (int j = 1; j <= jmax; ++j) {
            double df = d2 - d1;
            double hsp1 = 0.5 * d1 - 0.5 / sqrt(4 * dd[i] * dd[i] + 1 / (d1 * d1));
            double hsp2 = 0.5 * d2 - 0.5 / sqrt(4 * dd[i] * dd[i] + 1 / (d2 * d2));
            if (std::abs(hsp2 - hsp1) < 1e-25) break;
            double d3 = d1 + df * (hsp[i] / funcon_C::ev - hsp1) / (hsp2 - hsp1);
            d1 = d2;
            d2 = d3;
        }
        double gqq = pow(p4 * hpp / (funcon_C::ev * 48.0 * pow(qq[i], 4)), 0.2);
        double q1 = gqq;
        double q2 = gqq + 0.04;
        for (int j = 1; j <= jmax; ++j) {
            double qf = q2 - q1;
            double hpp1 = 0.25 * q1 - 0.5 / sqrt(4 * qq[i] * qq[i] + 1 / (q1 * q1)) +
                          0.25 / sqrt(8 * qq[i] * qq[i] + 1 / (q1 * q1));
            double hpp2 = 0.25 * q2 - 0.5 / sqrt(4 * qq[i] * qq[i] + 1 / (q2 * q2)) +
                          0.25 / sqrt(8 * qq[i] * qq[i] + 1 / (q2 * q2));
            if (std::abs(hpp2 - hpp1) < 1e-25) break;
            double q3 = q1 + qf * (hpp / funcon_C::ev - hpp1) / (hpp2 - hpp1);
            q1 = q2;
            q2 = q3;
        }
        am[i] = gss[i] / funcon_C::ev;
        ad[i] = d2;
        aq[i] = q2;
    }

    for (int i = 1; i <= 107; ++i) {
        if (am[i] < 1e-20) {
            if (gss[i] > 1e-20) am[i] = gss[i] / funcon_C::ev;
            else am[i] = 1.0;
        }
    }
    eisol[1] = uss[1];
    am[1] = gss[1] / funcon_C::ev;
    ad[1] = am[1];
    aq[1] = am[1];

    for (int i = 1; i <= 100; ++i) {
        if (f0sd_store[i] < 1e-20) f0sd[i] = 0.0;
        if (g2sd_store[i] < 1e-20) g2sd[i] = 0.0;
    }
    inid();
    am[102] = 1e-10;

    if (molkst_C::keywrd.find(" DEP ") == std::string::npos) return;
    if (!molkst_C::in_house_only) return;
    char num = molkst_C::method_pm7 ? '7' : '6';
    std::string file_name = "parameters_for_PM7_TS_C.F90";
    create_parameters_for_PMx_C(file_name, num);
}
