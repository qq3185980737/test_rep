// prtpka.cpp — C++ translation of MOPAC 2016 "prtpka.F90".
// Calculates and prints the pKa values for an organic molecule with an
// ionizable hydrogen attached to an oxygen atom.
#include "prtpka.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "parameters_for_AM1_C.h"
#include "chanel_C.h"
#include "to_screen.h"

using common_arrays_C::grad;
using common_arrays_C::xparam;
using common_arrays_C::p;
using common_arrays_C::nat;
using common_arrays_C::coord;
using molkst_C::numat;
using molkst_C::escf;
using molkst_C::mozyme;
using molkst_C::emin;
using molkst_C::mpack;
using molkst_C::moperr;
using parameters_C::tore;
using parameters_C::uss;
using parameters_C::upp;
using parameters_C::udd;
using parameters_C::zs;
using parameters_C::zp;
using parameters_C::zd;
using parameters_C::betas;
using parameters_C::betap;
using parameters_C::betad;
using parameters_C::gss;
using parameters_C::gsp;
using parameters_C::gpp;
using parameters_C::gp2;
using parameters_C::hsp;
using parameters_C::zsn;
using parameters_C::zpn;
using parameters_C::zdn;
using cosmo_C::iseps;
using cosmo_C::useps;
using cosmo_C::lpka;
using chanel_C::iw;

// Heavy external dependencies — provided by the host executable / tests.
extern void cosini(bool);
extern void coscav();
extern void mkbmat();
extern void moldat(int);
extern void calpar();
extern void switch_method();
extern void compfg(const std::vector<double>&, bool, double&, bool,
                   std::vector<double>&, bool);
extern void chrge(const std::vector<double>&, std::vector<double>&);
namespace linear_cosmo {
extern void ini_linear_cosmo();
extern void coscavz(std::vector<double>&, int);
}  // namespace linear_cosmo

// ---------------------------------------------------------------------
// Parameters_for_PKA — pKa parameter set + overrides the NDDO parameter
// tables for H, C, N, O, F, Cl, Br, I.
// ---------------------------------------------------------------------
void Parameters_for_PKA(double& c1, double& c2, double& c3) {
    c1 = -288.0530769;
    c2 = 28.6888717;
    c3 = 89.1172382;
    // Element 1 — Hydrogen
    uss[1] = -9.3555403;
    betas[1] = -2.6789813;
    zs[1] = 1.2458568;
    gss[1] = 14.5964443;
    // Element 6 — Carbon
    uss[6] = -49.9390029;
    upp[6] = -43.8233347;
    betas[6] = -12.6446141;
    betap[6] = -9.4614222;
    zs[6] = 1.6412075;
    zp[6] = 1.5810984;
    gss[6] = 16.7616538;
    gsp[6] = 12.3925395;
    gpp[6] = 11.2679258;
    gp2[6] = 11.0158335;
    hsp[6] = 0.4200022;
    // Element 7 — Nitrogen
    uss[7] = -55.4307466;
    upp[7] = -49.8247298;
    betas[7] = -21.7612447;
    betap[7] = -16.9335869;
    zs[7] = 1.4860676;
    zp[7] = 2.0420069;
    gss[7] = 7.1149144;
    gsp[7] = 7.2242656;
    gpp[7] = 14.8103679;
    gp2[7] = 11.4500819;
    hsp[7] = 4.5873409;
    // Element 8 — Oxygen
    uss[8] = -89.9472647;
    upp[8] = -70.7249961;
    betas[8] = -66.7497177;
    betap[8] = -21.7795050;
    zs[8] = 4.3991035;
    zp[8] = 2.1614926;
    gss[8] = 16.3875625;
    gsp[8] = 16.0293626;
    gpp[8] = 16.6098992;
    gp2[8] = 11.0339705;
    hsp[8] = 4.7920603;
    // Element 9 — Fluorine
    uss[9] = -141.5218991;
    upp[9] = -98.4596662;
    betas[9] = -69.7628706;
    betap[9] = -30.1591819;
    zs[9] = 5.9393528;
    zp[9] = 4.7675390;
    gss[9] = 12.5528331;
    gsp[9] = 19.8940208;
    gpp[9] = 8.8875462;
    gp2[9] = 12.5543486;
    hsp[9] = 3.5419742;
    // Element 17 — Chlorine
    uss[17] = -61.5736619;
    upp[17] = -54.5254717;
    udd[17] = -38.2581550;
    betas[17] = -0.5351767;
    betap[17] = -11.7700638;
    betad[17] = -4.0377510;
    zs[17] = 4.2327248;
    zp[17] = 1.2135427;
    zd[17] = 1.3240330;
    zsn[17] = 0.9562970;
    zpn[17] = 2.4640670;
    zdn[17] = 6.4103250;
    gss[17] = 9.4258401;
    gsp[17] = 5.3649222;
    gpp[17] = 10.2118976;
    gp2[17] = 9.1156399;
    hsp[17] = 4.9183864;
    // Element 35 — Bromine
    uss[35] = -46.3749516;
    upp[35] = -50.1929356;
    udd[35] = 7.0867380;
    betas[35] = -31.5606655;
    betap[35] = -9.1283278;
    betad[35] = -9.8391240;
    zs[35] = 5.0204513;
    zp[35] = 2.2579169;
    zd[35] = 1.5210310;
    zsn[35] = 3.0947770;
    zpn[35] = 3.0657640;
    zdn[35] = 2.8200030;
    gss[35] = 8.1108378;
    gsp[35] = 5.2179143;
    gpp[35] = 10.3707055;
    gp2[35] = 8.5089799;
    hsp[35] = 4.8307027;
    // Element 53 — Iodine
    uss[53] = -59.6492236;
    upp[53] = -56.2042321;
    udd[53] = -28.8226030;
    betas[53] = -30.0143361;
    betap[53] = -5.5364107;
    betad[53] = -7.6761070;
    zs[53] = 4.9453882;
    zp[53] = 2.4202135;
    zd[53] = 1.8751750;
    zsn[53] = 9.1352440;
    zpn[53] = 6.8881910;
    zdn[53] = 3.7915230;
    gss[53] = 7.7877000;
    gsp[53] = 9.4756338;
    gpp[53] = 10.5149673;
    gp2[53] = 8.1832762;
    hsp[53] = 4.7971922;
}

void prtpka(int* ipKa_sorted, double* pKa_sorted, int* ipKa_unsorted,
            double* pKa_unsorted, int& no) {
    double sum, c1, c2, c3;
    std::vector<double> dist(numat + 1, 0.0);
    std::vector<double> q(numat + 1, 0.0);
    int i, j, k, nh, loop, i_min = 0;
    double sum_min, sum_H_min;
    bool store_iseps, store_useps;
    for (i = 1; i <= numat; ++i) {
        ipKa_sorted[i] = 0;
        ipKa_unsorted[i] = 0;
        pKa_sorted[i] = 0.0;
        pKa_unsorted[i] = 0.0;
    }
    std::printf("METH pre: am1=%d pm7=%d pm6=%d mndo=%d rm1=%d mndod=%d pm3=%d pm7ts=%d\n",
                molkst_C::method_am1, molkst_C::method_pm7, molkst_C::method_pm6, molkst_C::method_mndo, molkst_C::method_rm1, molkst_C::method_mndod, molkst_C::method_pm3, molkst_C::method_pm7_ts);
    std::printf("METH uss8=%f upp8=%f gss8=%f uss1=%f ussam1_8=%f\n",
                parameters_C::uss[8], parameters_C::upp[8], parameters_C::gss[8],
                parameters_C::uss[1], parameters_for_AM1_C::ussam1[8]);
    std::fflush(stdout);
    // First, switch in parameters specific to the pKa calculation.
    Parameters_for_PKA(c1, c2, c3);
    // Re-calculate the charges
    store_iseps = iseps;
    store_useps = useps;
    iseps = true;
    useps = true;
    emin = 0.0;
    cosini(true);
    if (mozyme) linear_cosmo::ini_linear_cosmo();
    if (mozyme) {
        std::vector<double> cf(3 * numat + 1, 0.0);
        for (i = 1; i <= numat; ++i)
            for (k = 1; k <= 3; ++k) cf[3 * (i - 1) + k] = coord[k-1][i];
        linear_cosmo::coscavz(cf, numat);
        lpka = true;
    } else {
        coscav();
        mkbmat();
    }
    i = mpack;
    moldat(1);
    mpack = i;
    moperr = false;  // an error in moldat is not important here
    calpar();
    compfg(xparam, true, escf, true, grad, false);
    nh = 0;
    if (moperr) return;
    lpka = false;
    chrge(p, q);
    for (i = 1; i <= numat; ++i) q[i] = tore[nat[i]] - q[i];
    // Find the ionizable hydrogens
    for (i = 1; i <= numat; ++i) {
        if (nat[i] != 1) continue;
        sum_min = 2.0;
        for (j = 1; j <= numat; ++j) {
            if (nat[j] != 8) continue;
            sum = (coord[0][i] - coord[0][j]) * (coord[0][i] - coord[0][j]) +
                  (coord[1][i] - coord[1][j]) * (coord[1][i] - coord[1][j]) +
                  (coord[2][i] - coord[2][j]) * (coord[2][i] - coord[2][j]);
            if (sum < sum_min) {  // H within 1.4 A of an oxygen is attached
                sum_min = sum;
                i_min = j;
            }
        }
        if (sum_min < 1.9999) {
            // Exclude water
            sum_H_min = 2.0;
            for (k = 1; k <= numat; ++k) {
                if (nat[k] != 1 || k == i) continue;
                sum = (coord[0][k] - coord[0][i_min]) *
                          (coord[0][k] - coord[0][i_min]) +
                      (coord[1][k] - coord[1][i_min]) *
                          (coord[1][k] - coord[1][i_min]) +
                      (coord[2][k] - coord[2][i_min]) *
                          (coord[2][k] - coord[2][i_min]);
                if (sum < sum_H_min) sum_H_min = sum;
            }
            if (sum_H_min < 1.9999 && numat > 3) continue;
            ++nh;
            ipKa_unsorted[nh] = i;
            pKa_unsorted[nh] = q[i];
            dist[nh] = std::sqrt(sum_min);
        }
    }
    // Now sort into descending order
    for (i = 1; i <= numat; ++i) {
        ipKa_sorted[i] = 0;
        pKa_sorted[i] = 0.0;
    }
    if (nh == 0) {
        if (iw > 0) {
            to_screen(
                "\n"
                "          A request was made to print the pKa values for "
                "this system,\n"
                "          but there are no hydrogen atoms attached to an "
                "oxygen atom,\n"
                "          so the pKa calculation cannot be completed.");
        }
        no = 0;
        return;
    }
    // Convert charges and bond-lengths into pKa values
    for (i = 1; i <= nh; ++i)
        pKa_unsorted[i] = pKa_unsorted[i] * c1 + dist[i] * c2 + c3;
    // Sort pKa into order, select lowest four
    loop = 0;
    for (i = 1; i <= nh; ++i) {
        sum = 100.0;
        k = 0;
        for (j = 1; j <= nh; ++j)
            if (pKa_unsorted[j] < sum) {
                k = j;
                sum = pKa_unsorted[j];
            }
        if (k == 0) break;
        j = 1;
        if (loop > 0) {
            for (j = 1; j <= loop; ++j)
                if (ipKa_sorted[j] == ipKa_unsorted[k] &&
                    std::fabs(pKa_sorted[j] - pKa_unsorted[k]) < 0.2)
                    break;
        }
        if (j <= loop) continue;
        ++loop;
        ipKa_sorted[loop] = ipKa_unsorted[k];
        pKa_sorted[loop] = pKa_unsorted[k];
        pKa_unsorted[k] = 200.0 + pKa_unsorted[k];
    }
    for (i = 1; i <= loop; ++i) pKa_unsorted[i] -= 200.0;
    no = loop;
    std::printf("METH pre-switch uss8=%f\n", parameters_C::uss[8]); std::fflush(stdout);
    switch_method();
    std::printf("METH post-switch uss8=%f gss8=%f uss1=%f\n",
                parameters_C::uss[8], parameters_C::gss[8], parameters_C::uss[1]);
    std::fflush(stdout);
    iseps = store_iseps;
    useps = store_useps;
    i = mpack;
    moldat(1);
    mpack = i;
    moperr = false;
    calpar();
    emin = 0.0;
    compfg(xparam, true, escf, true, grad, false);
}
