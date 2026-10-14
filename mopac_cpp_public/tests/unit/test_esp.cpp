// test_esp.cpp — M10 esp.F90 batch: setup_esp/pdgrid/surfac/elesn/espfit/
// espplane.  All PASS + ASan clean is the acceptance gate.
#include "esp.h"
#include "esp_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "funcon_C.h"
#include "parameters_C.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using esp_C::nesp;
using esp_C::idip;
using esp_C::iz;
using esp_C::ipx;
using esp_C::isc;
using esp_C::is_esp;
using esp_C::icd;
using esp_C::ipe;
using esp_C::npr;
using esp_C::ncc;
using esp_C::scale;
using esp_C::cf;
using esp_C::rms;
using esp_C::rrms;
using esp_C::dx;
using esp_C::dy;
using esp_C::dz;
using esp_C::den;
using esp_C::es;
using esp_C::b_esp;
using esp_C::esp_array;
using esp_C::cespm;
using esp_C::co;
using esp_C::potpt;
using esp_C::qsc;
using esp_C::rad;
using esp_C::cc;
using esp_C::ex;
using esp_C::cen;
using esp_C::iam;
using esp_C::indc;

using common_arrays_C::geo;
using common_arrays_C::coord;
using common_arrays_C::labels;
using common_arrays_C::nat;
using common_arrays_C::na;
using common_arrays_C::nb;
using common_arrays_C::nc;
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using common_arrays_C::q;
using common_arrays_C::c;
using molkst_C::keywrd;
using molkst_C::moperr;
using molkst_C::numat;
using molkst_C::natoms;
using molkst_C::norbs;
using molkst_C::nopen;
using molkst_C::fract;
using molkst_C::nclose;
using molkst_C::method_mndo;
using funcon_C::a0;
using parameters_C::tore;
using parameters_C::zs;
using parameters_C::zp;

static int fails = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (cond) {                                                          \
            std::printf("PASS: %s\n", msg);                                  \
        } else {                                                             \
            std::printf("FAIL: %s (line %d)\n", msg, __LINE__);              \
            ++fails;                                                         \
        }                                                                    \
    } while (0)

// 1) setup_esp(1) allocates and initialises the esp_C arrays.
void t_setup_esp() {
    norbs = 6;
    numat = 3;
    setup_esp(1);
    CHECK(es.size() == 5001, "setup_esp: es sized 5001");
    CHECK(std::fabs(es[1] - 1.0e9) < 1.0, "setup_esp: es init 1e9");
    CHECK(co.size() == 4 && co[1].size() == 4, "setup_esp: co(3,numat)");
    CHECK(std::fabs(co[1][1] - 1.0e9) < 1.0, "setup_esp: co init 1e9");
    CHECK(potpt.size() == 4 && potpt[1].size() == 5001, "setup_esp: potpt(3,5000)");
    CHECK(std::fabs(potpt[1][1]) < 1e-12, "setup_esp: potpt init 0");
    CHECK(iam.size() == 37 && iam[1].size() == 3, "setup_esp: iam(6*norbs,2)");
    CHECK(indc.size() == 7, "setup_esp: indc(norbs)");
}

// 2) espfit: 1 atom + 3 points, no dipole constraint.
void t_espfit() {
    numat = 1;
    nesp = 3;
    iz = 0;
    idip = 0;
    dx = dy = dz = 0.0;
    // co(1 atom at origin)
    co[1][1] = 0.0;
    co[2][1] = 0.0;
    co[3][1] = 0.0;
    // three ESP probe points
    potpt[1][1] = 2.0; potpt[2][1] = 0.0; potpt[3][1] = 0.0;
    potpt[1][2] = 0.0; potpt[2][2] = 2.0; potpt[3][2] = 0.0;
    potpt[1][3] = 0.0; potpt[2][3] = 0.0; potpt[3][3] = 2.0;
    esp_array[1] = 0.5;
    esp_array[2] = -0.25;
    esp_array[3] = 0.125;
    espfit();
    CHECK(std::isfinite(cf) && cf > 0.0, "espfit: cf = fpc1*a0*fpc8*1e-10");
    CHECK(std::isfinite(q[1]), "espfit: q(1) finite");
    CHECK(std::fabs(q[1]) < 1.0e4, "espfit: q(1) bounded");
    CHECK(std::isfinite(rms) && std::isfinite(rrms), "espfit: rms/rrms finite");
    CHECK(rms > 0.0, "espfit: rms > 0");
    CHECK(std::fabs(rrms) < 10.0, "espfit: rrms bounded");
    // charge constraint row/col: q sums to iz == 0
    double sumq = q[1];
    CHECK(std::fabs(sumq - 0.0) < 1.0e-6, "espfit: charge constraint sum q = iz");
}

// 3) pdgrid: Williams surface on water (3 atoms).
void t_pdgrid() {
    natoms = 3;
    numat = 3;
    norbs = 6;
    labels.assign(4, 0);
    labels[1] = 8;
    labels[2] = 1;
    labels[3] = 1;
    nat.assign(4, 0);
    nat[1] = 8;
    nat[2] = 1;
    nat[3] = 1;
    na.assign(8, 0);
    nb.assign(8, 0);
    nc.assign(8, 0);
    na[2] = 1;
    na[3] = 1;
    geo.assign(4, std::vector<double>(4, 0.0));
    coord.assign(4, std::vector<double>(4, 0.0));
    // Z-matrix: O at origin; H1 0.96 A along +x (na[2]=1);
    // H2 0.96 A at 104.5 deg bond angle (na[3]=1, angle in radians).
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    geo[1][2] = 0.96; geo[2][2] = 0.0; geo[3][2] = 0.0;
    double ang = 104.5 * 3.14159265358979323846 / 180.0;
    geo[1][3] = 0.96;
    geo[2][3] = ang;
    geo[3][3] = 0.0;
    setup_esp(1);
    pdgrid();
    CHECK(nesp > 0, "pdgrid: nesp > 0");
    CHECK(nesp < 5000, "pdgrid: nesp within potpt bound");
    std::printf("DBG pdgrid: nesp=%d co O=(%g,%g,%g) H1=(%g,%g,%g) H2=(%g,%g,%g)\n",
                nesp, co[1][1], co[2][1], co[3][1], co[1][2], co[2][2],
                co[3][2], co[1][3], co[2][3], co[3][3]);
    if (nesp > 0) {
        // pdgrid uses MOPAC's own vdW table (2x Bondi-like radii):
        // vderw[8]=2.6 (O), vderw[1]=2.4 (H); shell=1.2.
        const double r8 = 2.6;
        const double r1 = 2.4;
        bool ok = true;
        int nfail = 0;
        for (int p = 1; p <= nesp; ++p) {
            bool in_shell = false;
            for (int a = 1; a <= 3; ++a) {
                double d = std::sqrt((co[1][a] - potpt[1][p]) * (co[1][a] - potpt[1][p]) +
                                     (co[2][a] - potpt[2][p]) * (co[2][a] - potpt[2][p]) +
                                     (co[3][a] - potpt[3][p]) * (co[3][a] - potpt[3][p]));
                double rv = (a == 1) ? r8 : r1;
                if (d >= rv - 1e-9 && d <= rv + 1.2 + 1e-9) in_shell = true;
            }
            if (!in_shell) {
                if (nfail < 3) {
                    std::printf("DBG bad pt p=%d (%g,%g,%g) d=", p, potpt[1][p],
                                potpt[2][p], potpt[3][p]);
                    for (int a = 1; a <= 3; ++a) {
                        double d = std::sqrt((co[1][a] - potpt[1][p]) * (co[1][a] - potpt[1][p]) +
                                             (co[2][a] - potpt[2][p]) * (co[2][a] - potpt[2][p]) +
                                             (co[3][a] - potpt[3][p]) * (co[3][a] - potpt[3][p]));
                        std::printf("%g ", d);
                    }
                    std::printf("\n");
                }
                ++nfail;
                ok = false;
            }
        }
        std::printf("DBG pdgrid bad count=%d/%d\n", nfail, nesp);
        CHECK(ok, "pdgrid: every point within a vdW shell");
    }
}

// 4) surfac: single carbon atom, contact points.
void t_surfac() {
    natoms = 1;
    numat = 1;
    norbs = 4;
    labels.assign(2, 0);
    labels[1] = 6;
    nat.assign(2, 0);
    nat[1] = 6;
    geo.assign(4, std::vector<double>(2, 0.0));
    coord.assign(4, std::vector<double>(2, 0.0));
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    setup_esp(1);
    scale = 1.4;
    den = 1.0;
    nesp = 0;
    surfac();
    std::printf("DBG surfac: nesp=%d rad[1]=%g\n", nesp, rad[1]);
    CHECK(nesp > 40, "surfac: nesp > 40 for single C (49 expected)");
    if (nesp > 0) {
        double r = 1.50 * 1.4;  // vander[C]=1.50
        bool ok = true;
        for (int p = 1; p <= nesp; ++p) {
            double d = std::sqrt(potpt[1][p] * potpt[1][p] +
                                 potpt[2][p] * potpt[2][p] +
                                 potpt[3][p] * potpt[3][p]);
            if (std::fabs(d - r) > 1e-6) ok = false;
        }
        CHECK(ok, "surfac: all points at vdW radius");
    }
}

// 5) elesn: single hydrogen atom, minimal full-chain run.
void t_elesn() {
    natoms = 1;
    numat = 1;
    norbs = 6;
    nclose = 1;
    nopen = 0;
    fract = 0.0;
    nfirst.assign(8, 0);
    nlast.assign(8, 0);
    nfirst[1] = 1;
    nlast[1] = 1;
    zs[1] = 1.0;  // H: Slater 1s exponent
    tore[1] = 1.0;  // H nuclear charge
    setup_esp(1);
    labels.assign(2, 0);
    labels[1] = 1;
    nat.assign(2, 0);
    nat[1] = 1;
    co[1][1] = 0.0; co[2][1] = 0.0; co[3][1] = 0.0;
    nesp = 2;
    potpt[1][1] = 1.0; potpt[2][1] = 0.0; potpt[3][1] = 0.0;
    potpt[1][2] = -1.0; potpt[2][2] = 0.0; potpt[3][2] = 0.0;
    // MO coefficient for one doubly occupied H 1s (density P=2*c*cT)
    c.assign(2, std::vector<double>(2, 0.0));
    c[1][1] = 1.0;
    keywrd = "PM7";
    elesn();
    std::printf("DBG elesn: es[1]=%g es[2]=%g esp_arr[1]=%g b_esp[1]=%g\n",
                es[1], es[2], esp_array[1], b_esp[1]);
    CHECK(std::isfinite(es[1]) && std::isfinite(es[2]), "elesn: es finite");
    CHECK(std::isfinite(esp_array[1]) && std::isfinite(esp_array[2]),
          "elesn: esp_array finite");
    CHECK(std::fabs(esp_array[1]) < 1.0e3, "elesn: esp_array bounded");
    CHECK(std::isfinite(b_esp[1]), "elesn: b_esp finite");
    CHECK(std::fabs(b_esp[1]) > 1.0e-9, "elesn: b_esp nonzero");
    // sanity: |es| should not be astronomically large for H
    CHECK(std::fabs(es[1]) < 1.0e5, "elesn: es bounded");
}

// 6) espplane: fills a plane grid when nesp==0.
void t_espplane() {
    nesp = 0;
    setup_esp(1);
    double xmin[4] = {0.0, 0.0, 0.0, 0.0};
    double step[4] = {0.0, 0.5, 0.25, 0.5};
    espplane(3, xmin, step, 4, 5);
    CHECK(nesp == 20, "espplane: 4*5 = 20 points");
    if (nesp == 20) {
        double z = 0.0 + 0.5 * (3 - 1);
        bool ok = true;
        for (int p = 1; p <= 20; ++p)
            if (std::fabs(potpt[3][p] - z) > 1e-12) ok = false;
        CHECK(ok, "espplane: all points on plane iplane");
    }
}

int main() {
    std::printf("== M10 esp batch ==\n");
    na.assign(16, 0);
    nb.assign(16, 0);
    nc.assign(16, 0);
    std::printf("[t1]\n"); std::fflush(stdout);
    t_setup_esp();
    std::printf("[t2]\n"); std::fflush(stdout);
    t_espfit();
    std::printf("[t3]\n"); std::fflush(stdout);
    t_pdgrid();
    std::printf("[t4]\n"); std::fflush(stdout);
    t_surfac();
    std::printf("[t5]\n"); std::fflush(stdout);
    t_elesn();
    std::printf("[t6]\n"); std::fflush(stdout);
    t_espplane();
    if (fails == 0) {
        std::printf("ALL PASS: M10 esp group\n");
        return 0;
    }
    std::printf("FAILURES: %d\n", fails);
    return 1;
}
