// test_periodic.cpp - periodic-table smoke test: X2 homonuclear dimer, PM7 RHF 1SCF single point.
// Usage: test_periodic <Z> [R_angstrom]   (default R=1.5)
// Chain: switch_method(PM7)->fordd->moldat->calpar->hcore->compfg(fulscf single point)
// Official baselines produced by Mopac_64.exe on identical .mop inputs (1SCF, R=1.5 A).
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "molkst_C.h"
#include "iter_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
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
using cosmo_C::nspa;
using cosmo_C::nppa;
using cosmo_C::ioldcv;

extern void switch_method();
extern void calpar();
extern void fordd();
extern void hcore();
extern void fbx();  // factorials + Pascal's triangle; must run before aijm/ddpo (F90: run_mopac setup)
using common_arrays_C::dxyz;

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: test_periodic <Z> [R]\n"); return 1; }
    int Z = atoi(argv[1]);
    double R = (argc > 2) ? atof(argv[2]) : 1.5;
    std::printf("pperiodic Z=%d R=%.4f\n", Z, R); std::fflush(stdout);

    natoms = 2;
    numat = 2;
    norbs = 0; mpack = 0;
    nvar = 1;
    dxyz.assign(3 * numat + 1, 0.0);
    molkst_C::tleft = 1e6;
    nat.assign(3, 0);
    nat[1] = Z; nat[2] = Z;
    nfirst.assign(3, 0); nlast.assign(3, 0);
    keywrd = "PM7";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;
    labels.assign(3, 0);
    labels[1] = Z; labels[2] = Z;
    geo.assign(4, std::vector<double>(4, 0.0));
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    geo[1][2] = 0.0; geo[2][2] = 0.0; geo[3][2] = R;
    na.assign(3, 0); nb.assign(3, 0); nc.assign(3, 0);
    na[2] = 0; nb[2] = 0; nc[2] = 0;
    coord.assign(4, std::vector<double>(numat + 1, 0.0));
    coord[0][1] = 0.0; coord[1][1] = 0.0; coord[2][1] = 0.0;
    coord[0][2] = 0.0; coord[1][2] = 0.0; coord[2][2] = R;
    xparam.assign(nvar + 1, 0.0);
    grad.assign(nvar + 1, 0.0);
    loc.assign(3, std::vector<int>(3 * natoms + 1, 0));
    loc[1][1] = 3; loc[2][1] = 2;   // xparam[1] = geo[3][2] (bond length along z)
    xparam[1] = R;

    method_pm7 = true;
    switch_method();   // loads PM7 parameter tables (uss/gss/zs/zp/zd/betas...)

    // Basis size per atom: mirror moldat natorb logic (zd -> 9, zp -> 4, zs -> 1)
    int natb = 1;
    if (parameters_C::zd[Z] > 1e-8) natb = 9;
    else if (parameters_C::zp[Z] > 1e-20) natb = 4;
    else if (parameters_C::zs[Z] > 1e-20) natb = 1;
    else { std::printf("pperiodic ERROR: no basis parameters for Z=%d\n", Z); return 2; }
    nfirst[1] = 1; nlast[1] = natb;
    nfirst[2] = natb + 1; nlast[2] = 2 * natb;
    norbs = 2 * natb;
    mpack = (norbs * (norbs + 1)) / 2;

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
    w.assign(6 * mpack * mpack + 100, 0.0);
    using iter_C::pold; using iter_C::pold2; using iter_C::pold3;
    using iter_C::pbold; using iter_C::pbold2; using iter_C::pbold3;
    pold.assign(6 * mpack + 1, 0.0); pold2.assign(6 * mpack + 1, 0.0); pold3.assign(6 * mpack + 1, 0.0);
    pbold.assign(6 * mpack + 1, 0.0); pbold2.assign(6 * mpack + 1, 0.0); pbold3.assign(6 * mpack + 1, 0.0);
    atmass.assign(numat + 1, 0.0);

    fbx();
    fordd();
    moldat(1);
    calpar();
    std::printf("pperiodic norbs=%d mpack=%d natorb[Z]=%d nclose=%d nelecs=%d\n",
                norbs, mpack, natb, molkst_C::nclose, molkst_C::nelecs);

    molkst_C::atheat = 0.0;
    for (int i = 1; i <= numat; ++i) molkst_C::atheat += parameters_C::eheat[nat[i]];
    double eat_sum = 0.0;
    for (int i = 1; i <= numat; ++i) eat_sum += parameters_C::eisol[nat[i]];
    molkst_C::atheat = molkst_C::atheat - eat_sum * funcon_C::fpc_9;
    hcore();

    std::vector<double> grad_sp(nvar + 1, 0.0);
    double escf_sp = 0.0;
    compfg(xparam, true, escf_sp, true, grad_sp, false);
    std::printf("pperiodic RESULT escf=%.5f kcal/mol EE=%.5f CC=%.5f TOTAL=%.5f eV\n",
                escf_sp, molkst_C::elect, molkst_C::enuclr, molkst_C::elect + molkst_C::enuclr);
    return 0;
}
