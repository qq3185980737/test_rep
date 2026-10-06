// test_ef_pm7.cpp - PM7 geometry optimisation of H2 via ef (P-RFO/QA driver).
// Chain: switch_method(PM7)->fordd->moldat->calpar->hcore->ef(xparam,escf)
// Official PM7 optimisation baseline: h2_pm7_opt.out (HOF=-32.01048 kcal/mol,
// TOTAL=-28.04697, EE=-42.79009, CC=14.74313).
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "molkst_C.h"
#include "iter_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "parameters_C.h"
#include "compfg.h"
#include "ef.h"
#include "moldat.h"

using molkst_C::numat;
using molkst_C::norbs;
using molkst_C::mpack;
using molkst_C::nvar;
using molkst_C::keywrd;
using molkst_C::mozyme;
using molkst_C::moperr;
using molkst_C::natoms;
using molkst_C::lm61;
using molkst_C::numcal;
using molkst_C::method_pm7;
using common_arrays_C::nat;
using common_arrays_C::coord;
using common_arrays_C::xparam;
using common_arrays_C::grad;
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using common_arrays_C::labels;
using common_arrays_C::geo;
using common_arrays_C::c;
using common_arrays_C::pdiag;
using common_arrays_C::uspd;
using common_arrays_C::h;
using common_arrays_C::w;
using common_arrays_C::f;
using common_arrays_C::fb;
using common_arrays_C::pa;
using common_arrays_C::pb;
using common_arrays_C::p;
using common_arrays_C::eigs;
using common_arrays_C::eigb;
using common_arrays_C::atmass;
using common_arrays_C::na;
using common_arrays_C::nb;
using common_arrays_C::nc;
using common_arrays_C::loc;
using parameters_C::tore;
using cosmo_C::nspa;
using cosmo_C::nppa;
using cosmo_C::ioldcv;

extern void switch_method();
extern void calpar();
extern void fordd();
extern void hcore();
using common_arrays_C::dxyz;

int main() {
    std::printf("efpm7 s0 enter\n"); std::fflush(stdout);
    numat = 2;
    norbs = 2;
    mpack = 3;
    nvar = 1;
    dxyz.assign(3 * numat + 1, 0.0);   // setup_mopac_arrays equivalent (deriv path)
    molkst_C::tleft = 1e6;             // F90 tleft=0 relies on tstep~0; ASan build is slow
    nat.assign(3, 0);
    nat[1] = 1; nat[2] = 1;
    nfirst.assign(3, 0); nlast.assign(3, 0);
    nfirst[1] = 1; nfirst[2] = 2;
    nlast[1] = 1; nlast[2] = 2;
    coord.assign(4, std::vector<double>(numat + 1, 0.0));
    coord[0][1] = 0.0;  coord[1][1] = 0.0;  coord[2][1] = 0.0;   // H1
    coord[0][2] = 0.755; coord[1][2] = 0.0;  coord[2][2] = 0.0;   // H2
    keywrd = "PM7";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;
    natoms = 2;
    labels.assign(3, 0);
    labels[1] = 1; labels[2] = 1;
    geo.assign(4, std::vector<double>(4, 0.0));
    geo[1][2] = 0.755; geo[2][2] = 0.0;
    na.assign(3, 0); nb.assign(3, 0); nc.assign(3, 0);
    na[2] = 1; nb[2] = 0; nc[2] = 0;
    c.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
    pdiag.assign(norbs + 1, 0.0);
    uspd.assign(norbs + 1, 0.0);
    h.assign(mpack + 1, 0.0);
    f.assign(mpack + 1, 0.0);
    fb.assign(mpack + 1, 0.0);
    pa.assign(mpack + 1, 0.0);
    pb.assign(mpack + 1, 0.0);
    p.assign(mpack + 1, 0.0);
    eigs.assign(norbs + 1, 0.0);
    eigb.assign(norbs + 1, 0.0);
    w.assign(132 + 2025 + 1, 0.0);
    using iter_C::pold; using iter_C::pold2; using iter_C::pold3;
    using iter_C::pbold; using iter_C::pbold2; using iter_C::pbold3;
    pold.assign(6 * mpack + 1, 0.0); pold2.assign(6 * mpack + 1, 0.0); pold3.assign(6 * mpack + 1, 0.0);
    pbold.assign(6 * mpack + 1, 0.0); pbold2.assign(6 * mpack + 1, 0.0); pbold3.assign(6 * mpack + 1, 0.0);
    atmass.assign(numat + 1, 0.0);
    xparam.assign(nvar + 1, 0.0);
    grad.assign(nvar + 1, 0.0);
    loc.assign(3, std::vector<int>(3 * natoms + 1, 0));
    loc[1][1] = 2; loc[2][1] = 1;
    tore[1] = 1;
    xparam[1] = 0.755;

    method_pm7 = true;
    switch_method();
    fordd();
    moldat(1);
    calpar();

    double escf = 0.0;
    molkst_C::atheat = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) molkst_C::atheat += parameters_C::eheat[nat[i]];
    double eat_sum = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) eat_sum += parameters_C::eisol[nat[i]];
    molkst_C::atheat = molkst_C::atheat - eat_sum * 23.061;
    hcore();

    std::printf("efpm7 s1 before ef: xparam[1]=%.5f escf_init=%.5f\n", xparam[1], escf);
    std::fflush(stdout);
    ef(xparam, escf);
    std::printf("efpm7 s2 after ef: escf=%.5f kcal/mol  (official PM7 opt: -32.01048)\n", escf);
    std::printf("efpm7 s3 xparam[1]=%.5f A  (bond length)\n", xparam[1]);
    std::printf("efpm7 s4 EE=%.5f CC=%.5f TOTAL=%.5f eV  (official: EE=-42.79009 CC=14.74313 TOTAL=-28.04697)\n",
                molkst_C::elect, molkst_C::enuclr, molkst_C::elect + molkst_C::enuclr);
    std::printf("efpm7 s5 grad="); for (int i = 1; i <= nvar; ++i) std::printf("%.5f ", grad[i]);
    std::printf("\n");
    return 0;
}
