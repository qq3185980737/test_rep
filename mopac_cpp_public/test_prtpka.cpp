// test_prtpka.cpp - smoke test: minimal water molecule through prtpka
#include <cstdio>
#include <string>
#include <vector>
#include "prtpka.h"
#include "molkst_C.h"
#include "parameters_for_AM1_C.h"
#include "parameters_for_PM7_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "iter_C.h"
#include "parameters_C.h"

using molkst_C::numat;
using molkst_C::norbs;
using molkst_C::mpack;
using molkst_C::nvar;
using molkst_C::keywrd;
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
using common_arrays_C::loc;
using parameters_C::tore;
using common_arrays_C::nb;
using common_arrays_C::nc;
using common_arrays_C::p;
using cosmo_C::nspa;
using cosmo_C::nppa;
using molkst_C::lm61;
using molkst_C::n2elec;
using molkst_C::numcal;
using iter_C::pold;
using iter_C::pold2;
using iter_C::pold3;
using iter_C::pbold;
using iter_C::pbold2;
using iter_C::pbold3;
using cosmo_C::ioldcv;
using molkst_C::mozyme;
using molkst_C::moperr;
using molkst_C::natoms;

int main() {
    std::printf("t0 main enter\n"); std::fflush(stdout);
    // Minimal water molecule (O + 2 H), 1-based indexing.
    numat = 3;
    norbs = 6;
    mpack = 21;
    nvar = 3;
    nat.assign(4, 0);
    nat[1] = 8; nat[2] = 1; nat[3] = 1;
    nfirst.assign(4, 0); nlast.assign(4, 0);
    nfirst[1] = 1; nfirst[2] = 5; nfirst[3] = 6;
    nlast[1] = 4; nlast[2] = 5; nlast[3] = 6;
    // coord layout: 0-based row = x/y/z, 1-based atom index (cosmo.cpp convention)
    coord.assign(4, std::vector<double>(numat + 1, 0.0));
    coord[0][1] = 0.0;  coord[1][1] = 0.0;  coord[2][1] = 0.0;   // O
    coord[0][2] = 0.96; coord[1][2] = 0.0;  coord[2][2] = 0.0;   // H1
    coord[0][3] = -0.24; coord[1][3] = 0.93; coord[2][3] = 0.0; // H2
    keywrd = "AM1";
    mozyme = false; moperr = false;
    nspa = 42; nppa = 42; lm61 = 21; ioldcv = 0;
    numcal = 1;   // hcore/compfg gate their first-call init on numcal
    // moldat() needs more molecule-level state.
    natoms = 3;
    labels.assign(4, 0);
    labels[1] = 8; labels[2] = 1; labels[3] = 1;
    // geo layout: geo[param][atom]; params 1=length,2=angle(rad),3=dihedral
    geo.assign(4, std::vector<double>(4, 0.0));
    geo[1][2] = 0.96; geo[2][2] = 0.0; geo[3][2] = 0.0;
    geo[1][3] = 0.96; geo[2][3] = 104.5 * 3.141592653589793 / 180.0; geo[3][3] = 0.0;
    // Z-matrix connectivity: atom 2 attached to 1; atom 3 attached to 1, angle ref 2
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
    pold.assign(mpack + 1, 0.0);
    pold2.assign(mpack + 1, 0.0);
    pold3.assign(mpack + 1, 0.0);
    pbold.assign(mpack + 1, 0.0);
    pbold2.assign(mpack + 1, 0.0);
    pbold3.assign(mpack + 1, 0.0);
    // w is allocated in setup_mopac_arrays.F90 as (n2elec + 2025); n2elec is
    // computed inside moldat() so allocate after it runs inside prtpka().
    // Water: n9=0, n4=1, n1=2 -> n2elec=132 -> w needs >= 132+2025+1.
    w.assign(132 + 2025 + 1, 0.0);
    atmass.assign(numat + 1, 0.0);
    xparam.assign(nvar + 1, 0.0);
    grad.assign(nvar + 1, 0.0);
    // loc table: loc[1][i] = atom, loc[2][i] = geo parameter row
    loc.assign(3, std::vector<int>(3 * natoms + 1, 0));  // F90 allocates loc(2,3*natoms) in setup_mopac_arrays; geout reads loc(1,n) up to nvar+1
    loc[1][1] = 2; loc[2][1] = 1;   // atom2 bond length
    loc[1][2] = 3; loc[2][2] = 1;   // atom3 bond length
    loc[1][3] = 3; loc[2][3] = 2;   // atom3 angle
    tore[1] = 1; tore[8] = 6;       // valence electrons for SCF electron count
    p.assign(mpack + 1, 0.0);

    int s[5] = {0}, u[5] = {0};
    double sp[5] = {0}, up[5] = {0};
    int no = 0;
    std::printf("t1 calling prtpka\n"); std::fflush(stdout);
    std::printf("AM1REF ussam1[1]=%.6f ussam1[8]=%.6f uppam1[8]=%.6f gssam1[8]=%.6f\n",
                parameters_for_AM1_C::ussam1[1], parameters_for_AM1_C::ussam1[8],
                parameters_for_AM1_C::uppam1[8], parameters_for_AM1_C::gssam1[8]);
    std::printf("PM7REF uss7[1]=%.6f uss7[8]=%.6f upp7[8]=%.6f gss7[8]=%.6f\n",
                parameters_for_PM7_C::uss7[1], parameters_for_PM7_C::uss7[8],
                parameters_for_PM7_C::upp7[8], parameters_for_PM7_C::gss7[8]);
    std::printf("METH now: am1=%d pm7=%d pm6=%d mndo=%d\n",
                molkst_C::method_am1, molkst_C::method_pm7, molkst_C::method_pm6, molkst_C::method_mndo);
    prtpka(s, sp, u, up, no);
    std::printf("prtpka no=%d PASS\n", no);
    std::printf("DIAG uspd: "); for (int i = 1; i <= 6; ++i) std::printf("%.5f ", uspd[i]);
    std::printf("\nDIAG eigs: "); for (int i = 1; i <= 6; ++i) std::printf("%.5f ", eigs[i]);
    std::printf("\nDIAG h: "); for (int i = 1; i <= 21; ++i) std::printf("%.4f ", h[i]);
    std::printf("\nDIAG p: "); for (int i = 1; i <= 21; ++i) std::printf("%.4f ", p[i]);
    std::printf("\nDIAG w1-20: "); for (int i = 1; i <= 20; ++i) std::printf("%.4f ", w[i]);
    std::printf("\nDIAG w100-119: "); for (int i = 100; i <= 119; ++i) std::printf("%.4f ", w[i]);
    std::printf("\n");
    return 0;
}
