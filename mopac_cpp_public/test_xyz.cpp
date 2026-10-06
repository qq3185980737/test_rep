// test_xyz.cpp - general N-atom PM7 RHF 1SCF single point.
// Usage: test_xyz N Z1 x1 y1 z1 Z2 x2 y2 z2 ...   (coordinates in Angstrom)
// Chain: switch_method(PM7)->fordd->moldat->calpar->hcore->compfg(fulscf single point)
#include <cstdio>
#include <cmath>
#include <cstdlib>
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
extern void fbx();
using common_arrays_C::dxyz;

int main(int argc, char** argv) {
    if (argc < 6) { std::printf("usage: test_xyz N Z1 x1 y1 z1 Z2 x2 y2 z2 ...\n"); return 1; }
    int N = atoi(argv[1]);
    if (argc < 2 + 4 * N) { std::printf("test_xyz: expected %d atom blocks, got %d args\n", N, argc - 2); return 2; }
    std::vector<int> Zs(N + 1, 0);
    std::vector<std::vector<double>> xyz(N + 1, std::vector<double>(3, 0.0));
    for (int i = 1; i <= N; ++i) {
        Zs[i] = atoi(argv[2 + 4 * (i - 1)]);
        xyz[i][0] = atof(argv[3 + 4 * (i - 1)]);
        xyz[i][1] = atof(argv[4 + 4 * (i - 1)]);
        xyz[i][2] = atof(argv[5 + 4 * (i - 1)]);
    }
    std::printf("pxyz N=%d:", N); std::fflush(stdout);
    for (int i = 1; i <= N; ++i) std::printf(" %d(%.3f,%.3f,%.3f)", Zs[i], xyz[i][0], xyz[i][1], xyz[i][2]);
    std::printf("\n"); std::fflush(stdout);

    natoms = N;
    numat = N;
    norbs = 0; mpack = 0;
    nvar = 1;
    dxyz.assign(3 * numat + 1, 0.0);
    molkst_C::tleft = 1e6;
    nat.assign(N + 1, 0);
    nfirst.assign(N + 1, 0); nlast.assign(N + 1, 0);
    keywrd = "PM7";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;
fprintf(stderr,"M1 labels\n");
    labels.assign(N + 1, 0);
fprintf(stderr,"M2 labels\n");
fprintf(stderr,"M3 geo\n");
    geo.assign(4, std::vector<double>(N + 1, 0.0));   // geo[component][atom] (matches gmetry)
fprintf(stderr,"M4 geo\n");
fprintf(stderr,"M5 coord\n");
    coord.assign(4, std::vector<double>(N + 1, 0.0));
fprintf(stderr,"M6 coord\n");
    na.assign(N + 1, 0); nb.assign(N + 1, 0); nc.assign(N + 1, 0);



    fprintf(stderr,"M11 pre-switch\n");
    method_pm7 = true;
    switch_method();
    fprintf(stderr,"M12 post-switch\n");
    fbx();
    fprintf(stderr,"M13 post-fbx\n");
    fordd();
    fprintf(stderr,"M14 post-fordd\n");
    fprintf(stderr,"M15 pre-ibasis nfirst.size=%zu coord.size=%zu geo.size=%zu\n", nfirst.size(), coord.size(), geo.size());
    int ibasis = 0;
    std::vector<int> natb(N + 1, 0);
    for (int i = 1; i <= N; ++i) {
        int Z = Zs[i];
        fprintf(stderr,"M16 i=%d Z=%d\n", i, Z);
        if (i >= (int)nat.size() || i >= (int)labels.size()) { fprintf(stderr,"OOB nat/labels i=%d size=%zu\n", i, nat.size()); return 9; }
        nat[i] = Z; labels[i] = Z;
        for (int a = 0; a < 3; ++a) { coord[a][i] = xyz[i][a]; geo[a + 1][i] = xyz[i][a]; }
        int nb_ = 1;
        if (parameters_C::zd[Z] > 1e-8) nb_ = 9;
        else if (parameters_C::zp[Z] > 1e-20) nb_ = 4;
        else if (parameters_C::zs[Z] > 1e-20) nb_ = 1;
        else { std::printf("pxyz ERROR: no basis for Z=%d\n", Z); return 3; }
        natb[i] = nb_;
        nfirst[i] = ibasis + 1; nlast[i] = ibasis + nb_;
        ibasis += nb_;
    }
    fprintf(stderr,"M17 pre-norbs ibasis=%d\n", ibasis);
    norbs = ibasis;
    mpack = (norbs * (norbs + 1)) / 2;
    fprintf(stderr,"M18 post-mpack norbs=%d mpack=%d\n", norbs, mpack);
fprintf(stderr,"M7 c\n");
    c.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
fprintf(stderr,"M8 c\n");
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
    xparam.assign(nvar + 1, 0.0);
    grad.assign(nvar + 1, 0.0);
fprintf(stderr,"M9 loc\n");
    loc.assign(3, std::vector<int>(3 * natoms + 1, 0));
fprintf(stderr,"M10 loc\n");
    xparam[1] = 1.0;   // dummy; moldat fills actual coordinates
    moldat(1);
    calpar();
    std::printf("pxyz norbs=%d mpack=%d nclose=%d nelecs=%d\n", norbs, mpack, molkst_C::nclose, molkst_C::nelecs);

    molkst_C::atheat = 0.0;
    for (int i = 1; i <= numat; ++i) molkst_C::atheat += parameters_C::eheat[nat[i]];
    double eat_sum = 0.0;
    for (int i = 1; i <= numat; ++i) eat_sum += parameters_C::eisol[nat[i]];
    molkst_C::atheat = molkst_C::atheat - eat_sum * funcon_C::fpc_9;
    hcore();

    std::vector<double> grad_sp(nvar + 1, 0.0);
    double escf_sp = 0.0;
    compfg(xparam, true, escf_sp, true, grad_sp, false);
    std::printf("pxyz RESULT escf=%.5f kcal/mol EE=%.5f CC=%.5f TOTAL=%.5f eV\n",
                escf_sp, molkst_C::elect, molkst_C::enuclr, molkst_C::elect + molkst_C::enuclr);
    return 0;
}
