// test_o2_rhf.cpp - PM7 geometry optimisation of O2 (RHF singlet) via ef.
// Chain: switch_method(PM7)->fordd->moldat->calpar->hcore->ef(xparam,escf)
// Official PM7 RHF O2 baseline (o2_rhf.out):
//   HOF=25.67492 kcal/mol, TOTAL=-585.91459, EE=-949.10711, CC=363.19252 eV,
//   opt O-O bond = 1.135136 A (started 1.2 A).
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
    std::printf("efo2 s0 enter\n"); std::fflush(stdout);
    numat = 2;
    norbs = 8;        // O: s+px+py+pz = 4 AO per atom
    mpack = 36;       // 8*9/2
    nvar = 1;
    dxyz.assign(3 * numat + 1, 0.0);
    molkst_C::tleft = 1e6;
    nat.assign(3, 0);
    nat[1] = 8; nat[2] = 8;
    nfirst.assign(3, 0); nlast.assign(3, 0);
    nfirst[1] = 1; nfirst[2] = 5;
    nlast[1] = 4; nlast[2] = 8;
    coord.assign(4, std::vector<double>(numat + 1, 0.0));
    coord[0][1] = 0.0; coord[1][1] = 0.0; coord[2][1] = 0.0;   // O1
    coord[0][2] = 0.0; coord[1][2] = 0.0; coord[2][2] = 1.2;   // O2
    keywrd = "PM7";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;
    natoms = 2;
    labels.assign(3, 0);
    labels[1] = 8; labels[2] = 8;
    geo.assign(4, std::vector<double>(4, 0.0));
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    geo[1][2] = 0.0; geo[2][2] = 0.0; geo[3][2] = 1.2;
    na.assign(3, 0); nb.assign(3, 0); nc.assign(3, 0);
    na[2] = 0; nb[2] = 0; nc[2] = 0;
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
    loc[1][1] = 3; loc[2][1] = 2;   // xparam[1]=geo[3][2] (O2 z = bond length, along z like official)
    tore[1] = 1; tore[8] = 6;       // Zcore(O)=6 (MOPAC tore(8))
    xparam[1] = 1.2;

    method_pm7 = true;
    switch_method();
    fordd();
    moldat(1);
    calpar();

    std::printf("efo2 p1 calpar: uss[8]=%.6f upp[8]=%.6f betas[8]=%.6f betap[8]=%.6f\n", parameters_C::uss[8], parameters_C::upp[8], parameters_C::betas[8], parameters_C::betap[8]);
    std::printf("efo2 p2 zs[8]=%.6f zp[8]=%.6f gss[8]=%.6f gsp[8]=%.6f gpp[8]=%.6f gp2[8]=%.6f hsp[8]=%.6f\n", parameters_C::zs[8], parameters_C::zp[8], parameters_C::gss[8], parameters_C::gsp[8], parameters_C::gpp[8], parameters_C::gp2[8], parameters_C::hsp[8]);
    std::printf("efo2 p3 tore[8]=%.4f ios[8]=%d iop[8]=%d iod[8]=%d npq(8)=%d,%d,%d\n", parameters_C::tore[8], parameters_C::ios[8], parameters_C::iop[8], parameters_C::iod[8], parameters_C::npq[8][1], parameters_C::npq[8][2], parameters_C::npq[8][3]);

    double escf = 0.0;
    molkst_C::atheat = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) molkst_C::atheat += parameters_C::eheat[nat[i]];
    double eat_sum = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) eat_sum += parameters_C::eisol[nat[i]];
    molkst_C::atheat = molkst_C::atheat - eat_sum * 23.061;
    hcore();

    std::printf("efo2 s1 before ef: xparam[1]=%.5f escf_init=%.5f\n", xparam[1], escf);
    std::printf("efo2 h h[1..10]=%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10]);
    std::printf("efo2 pd pdiag[1..8]=%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", pdiag[1],pdiag[2],pdiag[3],pdiag[4],pdiag[5],pdiag[6],pdiag[7],pdiag[8]);
    std::printf("efo2 hcore h[1..12]="); for (int i = 1; i <= 12; ++i) std::printf("%.4f ", h[i]); std::printf("\n");
    std::fflush(stdout);
    std::vector<double> grad_sp(nvar + 1, 0.0);
    double escf_sp = 0.0;
    compfg(xparam, true, escf_sp, true, grad_sp, false);
    std::printf("efo2 sp singlepoint: escf=%.5f kcal/mol EE=%.5f CC=%.5f TOTAL=%.5f eV  (official single-point EE=-934.89742 CC=349.12843 TOTAL=-585.76899)\n",
                escf_sp, molkst_C::elect, molkst_C::enuclr, molkst_C::elect + molkst_C::enuclr);
    std::printf("efo2 fock f[1..36]="); for (int i = 1; i <= 36; ++i) std::printf("%.6f ", f[i]); std::printf("\n");
    std::printf("efo2 dens p[1..36]="); for (int i = 1; i <= 36; ++i) std::printf("%.6f ", p[i]); std::printf("\n");
    std::printf("efo2 h h[1..36]="); for (int i = 1; i <= 36; ++i) std::printf("%.6f ", h[i]); std::printf("\n");
    std::fflush(stdout);
    ef(xparam, escf);
    std::printf("efo2 s2 after ef: escf=%.5f kcal/mol  (official RHF O2: 25.67492)\n", escf);
    std::printf("efo2 s3 xparam[1]=%.5f A  (bond length; official opt 1.13514)\n", xparam[1]);
    std::printf("efo2 s4 EE=%.5f CC=%.5f TOTAL=%.5f eV  (official: EE=-949.10711 CC=363.19252 TOTAL=-585.91459)\n",
                molkst_C::elect, molkst_C::enuclr, molkst_C::elect + molkst_C::enuclr);
    std::printf("efo2 s5 grad="); for (int i = 1; i <= nvar; ++i) std::printf("%.5f ", grad[i]);
    std::printf("\n");
    return 0;
}
