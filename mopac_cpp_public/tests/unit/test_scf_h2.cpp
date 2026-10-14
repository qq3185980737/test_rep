// test_scf_h2.cpp - clean AM1 single-point SCF chain for H2 (2 atoms, 2 electrons)
// Same chain as test_scf_am1 (main->switch_method->calpar->hcore->compfg->iter).
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
using molkst_C::method_am1;
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

int main() {
    std::printf("h2 s0 enter\n"); std::fflush(stdout);
    numat = 2;
    norbs = 2;
    mpack = 3;
    nvar = 1;
    nat.assign(3, 0);
    nat[1] = 1; nat[2] = 1;
    nfirst.assign(3, 0); nlast.assign(3, 0);
    nfirst[1] = 1; nfirst[2] = 2;
    nlast[1] = 1; nlast[2] = 2;
    coord.assign(4, std::vector<double>(numat + 1, 0.0));
    coord[0][1] = 0.0;  coord[1][1] = 0.0;  coord[2][1] = 0.0;   // H1
    coord[0][2] = 0.74; coord[1][2] = 0.0;  coord[2][2] = 0.0;   // H2
    keywrd = "AM1 HCORE";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;
    natoms = 2;
    labels.assign(3, 0);
    labels[1] = 1; labels[2] = 1;
    geo.assign(4, std::vector<double>(4, 0.0));
    geo[1][2] = 0.74; geo[2][2] = 0.0;
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
    xparam[1] = 0.74;

    method_am1 = true;
    switch_method();
    fordd();
    std::printf("h2 s1 uss1=%.5f gss1=%.5f betas1=%.5f\n",
                parameters_C::uss[1], parameters_C::gss[1], parameters_C::betas[1]);
    std::fflush(stdout);

    moldat(1);
    std::printf("h2 s2 uspd="); for (int i = 1; i <= 2; ++i) std::printf("%.5f ", uspd[i]);
    std::printf("\n"); std::fflush(stdout);
    calpar();
    std::printf("h2 s2.6 after calpar\n"); std::fflush(stdout);

    extern void hcore();
    std::printf("h2 sPO po1H=%.6f po7H=%.6f\n",
        common_arrays_C::po[1][1], common_arrays_C::po[7][1]);
    std::fflush(stdout);
    hcore();
    std::printf("h2 sH h[1..3]="); for (int i = 1; i <= 3; ++i) std::printf("%.6f ", h[i]);
    std::printf("\n"); std::fflush(stdout);

    // atheat = sum(eheat) - sum(eisol)*fpc_9   (run_mopac.F90 / run_mopac.cpp:556-562)
    molkst_C::atheat = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) molkst_C::atheat += parameters_C::eheat[nat[i]];
    double eat_sum = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) eat_sum += parameters_C::eisol[nat[i]];
    molkst_C::atheat = molkst_C::atheat - eat_sum * 23.061;
    std::printf("h2 s9b atheat=%.5f kcal/mol\n", molkst_C::atheat);

    double escf = 0.0;
    compfg(xparam, true, escf, true, grad, false);
    std::printf("h2 s3 escf(HOF)=%.5f kcal/mol\n", escf);
    std::printf("h2 s4 eigs="); for (int i = 1; i <= 2; ++i) std::printf("%.4f ", eigs[i]);
    std::printf("\n");
    std::printf("h2 s5 p="); for (int i = 1; i <= 3; ++i) std::printf("%.4f ", p[i]);
    std::printf("\n");
    std::printf("h2 s6 f="); for (int i = 1; i <= 3; ++i) std::printf("%.4f ", f[i]);
    std::printf("\n");
    std::printf("h2 s6 f="); for (int i = 1; i <= 3; ++i) std::printf("%.4f ", f[i]);
    std::printf("\n");
    std::printf("h2 s7 enuclr(CC)=%.5f eV elect(EE)=%.5f eV  (official 0.74: CC=13.37667 EE=-40.84818)\n", molkst_C::enuclr, molkst_C::elect);
    std::printf("h2 s8 TOTAL=%.5f eV  (official 0.74: -27.47151)\n", molkst_C::elect + molkst_C::enuclr);
    std::printf("h2 s9 escf-HOF=%.5f kcal/mol  (official 0.74: -3.68820)\n",
                (molkst_C::elect + molkst_C::enuclr) * 23.061 + molkst_C::atheat);
    return 0;
}
