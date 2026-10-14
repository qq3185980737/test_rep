// test_scf_am1.cpp - clean AM1 single-point SCF chain (no pKa overrides)
// Goal: reproduce official Mopac_64.exe water AM1 1SCF result:
//   FINAL HEAT OF FORMATION = -59.23271 kcal/mol
//   ELECTRONIC ENERGY      = -493.33914 eV
//   CORE-CORE REPULSION    =  144.77675 eV
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
extern void fordd();   // run_mopac.F90:102 - fills indpp/inddd/inddp index tables

int main() {
    std::printf("s0 enter\n"); std::fflush(stdout);
    numat = 3;
    norbs = 6;
    mpack = 21;
    nvar = 3;
    nat.assign(4, 0);
    nat[1] = 8; nat[2] = 1; nat[3] = 1;
    nfirst.assign(4, 0); nlast.assign(4, 0);
    nfirst[1] = 1; nfirst[2] = 5; nfirst[3] = 6;
    nlast[1] = 4; nlast[2] = 5; nlast[3] = 6;
    coord.assign(4, std::vector<double>(numat + 1, 0.0));
    coord[0][1] = 0.0;  coord[1][1] = 0.0;  coord[2][1] = 0.0;   // O
    coord[0][2] = 0.96; coord[1][2] = 0.0;  coord[2][2] = 0.0;   // H1
    coord[0][3] = -0.24; coord[1][3] = 0.93; coord[2][3] = 0.0; // H2
    keywrd = "AM1 HCORE";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;
    natoms = 3;
    labels.assign(4, 0);
    labels[1] = 8; labels[2] = 1; labels[3] = 1;
    geo.assign(4, std::vector<double>(4, 0.0));
    geo[1][2] = 0.96; geo[2][2] = 0.0; geo[3][2] = 0.0;
    geo[1][3] = 0.96; geo[2][3] = 104.5 * 3.141592653589793 / 180.0; geo[3][3] = 0.0;
    na.assign(4, 0); nb.assign(4, 0); nc.assign(4, 0);
    na[2] = 1; nb[2] = 0; nc[2] = 0;
    na[3] = 1; nb[3] = 2; nc[3] = 0;
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
    // iter_C scratch (setup_mopac_arrays.cpp:108 normally does this)
    using iter_C::pold; using iter_C::pold2; using iter_C::pold3;
    using iter_C::pbold; using iter_C::pbold2; using iter_C::pbold3;
    pold.assign(6 * mpack + 1, 0.0); pold2.assign(6 * mpack + 1, 0.0); pold3.assign(6 * mpack + 1, 0.0);
    pbold.assign(6 * mpack + 1, 0.0); pbold2.assign(6 * mpack + 1, 0.0); pbold3.assign(6 * mpack + 1, 0.0);
    atmass.assign(numat + 1, 0.0);
    xparam.assign(nvar + 1, 0.0);
    grad.assign(nvar + 1, 0.0);
    loc.assign(3, std::vector<int>(3 * natoms + 1, 0));
    loc[1][1] = 2; loc[2][1] = 1;
    loc[1][2] = 3; loc[2][2] = 1;
    loc[1][3] = 3; loc[2][3] = 2;
    tore[1] = 1; tore[8] = 6;
    xparam[1] = 0.96; xparam[2] = 0.96; xparam[3] = 104.5 * 3.141592653589793 / 180.0;

    // --- AM1 method selection (as run_mopac does: parse keywrd, then switch) ---
    method_am1 = true;               // keywrd contains "AM1"
    switch_method();                 // copies AM1 parameter tables into parameters_C
    fordd();                         // index tables for two-center d integrals
    std::printf("s1 uss8=%.5f upp8=%.5f gss8=%.5f uss1=%.5f\n",
                parameters_C::uss[8], parameters_C::upp[8],
                parameters_C::gss[8], parameters_C::uss[1]);
    std::fflush(stdout);

    moldat(1);
    std::printf("s2 uspd="); for (int i = 1; i <= 6; ++i) std::printf("%.5f ", uspd[i]);
    std::printf("\n"); std::fflush(stdout);
    std::printf("s2.5 before calpar\n"); std::fflush(stdout);
    calpar();
    std::printf("s2.6 after calpar\n"); std::fflush(stdout);

    extern void hcore();
    std::printf("sPO po1H=%.6f po9H=%.6f po1O=%.6f po9O=%.6f po7H=%.6f\n",
        common_arrays_C::po[1][1], common_arrays_C::po[9][1],
        common_arrays_C::po[1][8], common_arrays_C::po[9][8], common_arrays_C::po[7][1]);
    std::fflush(stdout);
    hcore();
    std::printf("sH h[1..21]="); for (int i = 1; i <= 21; ++i) std::printf("%.4f ", h[i]);
    std::printf("\n"); std::fflush(stdout);

    // atheat = sum(eheat) - sum(eisol)*fpc_9   (run_mopac.F90 / run_mopac.cpp:556-562)
    molkst_C::atheat = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) molkst_C::atheat += parameters_C::eheat[nat[i]];
    double eat_sum = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) eat_sum += parameters_C::eisol[nat[i]];
    molkst_C::atheat = molkst_C::atheat - eat_sum * 23.061;
    std::printf("s3b atheat=%.5f kcal/mol\n", molkst_C::atheat);

    double escf = 0.0;
    std::printf("s2.7 before compfg\n"); std::fflush(stdout);
    compfg(xparam, true, escf, true, grad, false);
    std::printf("s2.8 after compfg\n"); std::fflush(stdout);
    std::printf("sH h[1..21]="); for (int i = 1; i <= 21; ++i) std::printf("%.4f ", h[i]);
    std::printf("\n"); std::fflush(stdout);
    std::printf("sW w[1..23]="); for (int i = 1; i <= 23; ++i) std::printf("%.4f ", w[i]);
    std::printf("\n"); std::fflush(stdout);
    std::printf("s3 escf(HOF)=%.5f kcal/mol\n", escf);
    std::printf("s4 eigs="); for (int i = 1; i <= 6; ++i) std::printf("%.4f ", eigs[i]);
    std::printf("\n");
    std::printf("s5 pdiag="); for (int i = 1; i <= 6; ++i) std::printf("%.4f ", pdiag[i]);
    std::printf("\n");
    return 0;
}
