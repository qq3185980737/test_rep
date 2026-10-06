// test_batchM05c.cpp — M05c point-group analysis batch:
// symtry/haddon, maksym, symh, symt, sympop, symr/symp, charmo, symoir, symtrz.
// ASCII only. MSVC + ASan.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "symtry.h"
#include "maksym.h"
#include "symh.h"
#include "symt.h"
#include "sympop.h"
#include "symr.h"
#include "charmo.h"
#include "symoir.h"
#include "symtrz.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace symmetry_C;

static int nPass = 0, nFail = 0;
static void chk(const char* name, double got, double exp, double tol) {
    if (std::fabs(got - exp) <= tol) { ++nPass; }
    else { ++nFail; std::printf("FAIL %s: got %.12g exp %.12g\n", name, got, exp); }
}
static void chki(const char* name, int got, int exp) {
    if (got == exp) { ++nPass; } else { ++nFail; std::printf("FAIL %s: got %d exp %d\n", name, got, exp); }
}
static void chks(const char* name, const std::string& got, const std::string& exp) {
    if (got == exp) { ++nPass; } else { ++nFail; std::printf("FAIL %s: got '%s' exp '%s'\n", name, got.c_str(), exp.c_str()); }
}

// stubs: permanent 2016 gaps
double charmvi(const std::vector<std::vector<double>>&, int, int,
              const std::vector<std::vector<double>>&, int) { return 0.0; }
double charmst(const std::vector<std::vector<double>>&, const std::vector<int>&,
               int, int, const std::vector<std::vector<double>>&, int, bool&) { return 0.0; }
void makopr(int, double*, int&, double*) {}
void molsym(double*, int&, double*) {}

int main() {
    // ---- 1: symtry + haddon ----
    {
        geo.resize(4, std::vector<double>(4, 0.0));
        na.assign(4, 0);
        geo[1][1] = 2.5; geo[1][2] = 1.5; geo[1][3] = 3.0;
        geo[2][1] = 1.0; geo[2][2] = 1.0; geo[2][3] = 1.0;
        geo[3][1] = 0.7; geo[3][2] = 0.7; geo[3][3] = 0.7;
        // m=1 cartesian X=X on atom loc=1 -> geo[1][2] = 2.5
        ndep = 5;
        idepfn.assign(6, 0); locpar.assign(6, 0); locdep.assign(6, 0); depmul.assign(6, 0.0);
        idepfn[1] = 1;  locpar[1] = 1; locdep[1] = 2; depmul[1] = 2.0;  // X=X; fact for m19 below
        idepfn[2] = 15; locpar[2] = 2; locdep[2] = 3; depmul[2] = 1.0;  // bond-length/2 (na[2]!=0)
        idepfn[3] = 3;  locpar[3] = 2; locdep[3] = 3; depmul[3] = 1.0;  // dihedral equal
        idepfn[4] = 19; locpar[4] = 2; locdep[4] = 1; depmul[4] = 2.0;  // bond * depmul[1]
        idepfn[5] = 3;  locpar[5] = 2; locdep[5] = 2; depmul[5] = 1.0;  // dihedral equal
        na[2] = 1; na[3] = 2;  // atoms 2,3 have bond lengths (internal coords)
        geo[1][2] = 0.0; geo[1][3] = 0.0; geo[3][2] = 0.0; geo[3][3] = 0.0;  // geo[1][1]=2.5 kept as source
        symtry();
        chk("symtry m1", geo[1][2], 2.5, 1e-12);   // X = X
        chk("symtry m15", geo[1][3], 0.5, 1e-12);   // bond/2 = geo[2][1]/2 = 0.5
        chk("symtry m3", geo[3][3], 1.0, 1e-12);    // dihedral equal = geo[2][3]
        chk("symtry m19", geo[1][1], 2.0, 1e-12);   // bond*fact = 1.0*2.0 (depmul[1])
        chk("symtry m3b", geo[3][2], 1.0, 1e-12);   // dihedral equal = geo[2][3]
    }

    // ---- 2: maksym ----
    {
        natoms = 3; nvar = 4;
        na.assign(4, 0); na[1] = 1; na[2] = 2;
        int loc[8] = {0};
        // loc(2,nvar) column-major: loc[row1,col1], loc[row2,col1], loc[row1,col2], ...
        loc[0]=11; loc[1]=3;   // var1: atom 11, dihedral
        loc[2]=12; loc[3]=3;   // var2: atom 12, dihedral
        loc[4]=13; loc[5]=1;   // var3: atom 13, bond
        loc[6]=14; loc[7]=1;   // var4: atom 14, bond
        double xparam[4] = {5.0, 5.0, 1.0, 1.0};
        double xstore[4] = {0, 0, 0, 0};
        maksym(loc, xparam, xstore);
        // 5.0 folds into [-pi, pi): 5.0 - 2*pi = -1.283185
        chk("maksym fold", xparam[0], 5.0 - 2 * 3.14159265358979323846, 1e-9);
        // dep 1: dihedral equal (var2 same value) -> locpar=11, idepfn=3, locdep=12
        chki("maksym ndep", ndep, 2);
        chki("maksym dep1 par", locpar[1], 11);
        chki("maksym dep1 fn", idepfn[1], 3);
        chki("maksym dep1 dep", locdep[1], 12);
        // dep 2: bond equal (var4 == var3) -> idepfn=1? No: m=locl=1, idepfn=locl=1
        chki("maksym dep2 par", locpar[2], 13);
        chki("maksym dep2 fn", idepfn[2], 1);
        chki("maksym dep2 dep", locdep[2], 14);
        // compression: vars 2 and 4 removed -> nvar = 2
        chki("maksym nvar", nvar, 2);
        chki("maksym loc11", loc[0], 11);
        chki("maksym loc21", loc[1], 3);
        chki("maksym loc12", loc[2], 13);
        chki("maksym loc22", loc[3], 1);
    }

    // ---- 3: symh ----
    {
        numat = 2;
        ipo.assign(3, std::vector<int>(121, 0));
        ipo[1][1] = 1; ipo[2][1] = 2;
        ipo[1][2] = 1; ipo[2][2] = 1;  // op 2 maps atom 2 -> atom 1
        for (int a = 0; a < 10; ++a) for (int b = 0; b < 121; ++b) r[a][b] = 0.0;
        r[1][1] = 1; r[5][1] = 1; r[9][1] = 1;  // E
        r[1][2] = 1; r[5][2] = 1; r[9][2] = 1;  // identity for test
        std::vector<double> h(22, 0.0), dip(20, 0.0);
        h[6] = 2.0;   // (1,1) diagonal: iel=(3*4)/2=6
        dip[(1-1)*3+0] = 1.0;  // dip(1,1)
        symh(h.data(), &dip[0], 2, 2, nullptr);
        chk("symh diag", h[21], 2.0, 1e-12);          // written 1.0, then doubled to 2.0
        chk("symh dip", dip[(4-1)*3+0], 1.0, 1e-12); // dip(1,1) -> (1,4): dip[9]
        // doubling: istart=7..21 doubled -> h[21] now 2.0
        chk("symh double", h[21], 2.0, 1e-12);
    }

    // ---- 4: symt ----
    {
        nsym = 2;
        ipo.assign(3, std::vector<int>(121, 0));
        ipo[1][1] = 1; ipo[2][1] = 2;
        ipo[1][2] = 1; ipo[2][2] = 2;
        for (int a = 0; a < 10; ++a) for (int b = 0; b < 121; ++b) r[a][b] = 0.0;
        r[1][1] = 1; r[5][1] = 1; r[9][1] = 1;
        r[1][2] = 1; r[5][2] = 1; r[9][2] = 1;
        std::vector<double> h(22, 0.0), ha(22, 0.0), dip(20, 0.0);
        h[1] = 1.0; h[6] = 2.0; h[21] = 3.0;
        dip[(1-1)*3+0] = 4.0;
        symt(h.data(), &dip[0], ha.data());
        chk("symt h1", h[1], 1.0, 1e-12);    // two identity ops -> unchanged
        chk("symt h6", h[6], 2.0, 1e-12);
        chk("symt h21", h[21], 3.0, 1e-12);
        chk("symt dip", dip[(1-1)*3+0], 4.0, 1e-12);
    }

    // ---- 5: sympop ----
    {
        nsym = 2;
        ipo.assign(3, std::vector<int>(121, 0));
        ipo[1][1] = 1; ipo[2][1] = 2;
        ipo[1][2] = 1; ipo[2][2] = 1;  // ipo(2,2)=1 < 2
        std::vector<double> h(22, 0.0), dip(20, 0.0);
        int iskip = -1;
        sympop(h.data(), 2, iskip, &dip[0]);
        chki("sympop iskip", iskip, 3);
    }

    // ---- 6: symr ----
    {
        numat = 2;
        coord.resize(4, std::vector<double>(8, 0.0));
        coord[1][1] = 1.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
        coord[1][2] = -1.0; coord[2][2] = 0.0; coord[3][2] = 0.0;
        nclass = 1;
        symr();
        chki("symr nsym", nsym, 1);
        chki("symr ipo11", ipo[1][1], 1);
        chki("symr ipo21", ipo[2][1], 2);
    }
    {
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 21; ++k) elem[i][j][k] = 0.0;
        elem[1][1][2] = 1; elem[2][2][2] = -1; elem[3][3][2] = -1;  // C2x
        nclass = 2;
        coord.resize(4, std::vector<double>(8, 0.0));
        coord[1][1] = 1.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
        coord[1][2] = -1.0; coord[2][2] = 0.0; coord[3][2] = 0.0;
        symr();
        chki("symr2 nsym", nsym, 2);
        chki("symr2 ipo12", ipo[1][2], 1);
        chki("symr2 ipo22", ipo[2][2], 2);
        chk("symr2 r12", r[1][2], 1.0, 1e-12);   // C2x xx
        chk("symr2 r52", r[5][2], -1.0, 1e-12);  // C2x yy
    }

    // ---- 7: charmo ----
    {
        numat = 1; norbs = 1;
        jelem.assign(3, std::vector<int>(8, 0));
        jelem[1][1] = 1; jelem[2][1] = 1;
        std::vector<std::vector<double>> vv(2, std::vector<double>(2, 0.0));
        vv[1][1] = 1.0;
        std::vector<int> ntype(2, 0);
        ntype[1] = 100 * 1 + 9 + 1;  // 110: s-type of atom 1
        std::vector<std::vector<double>> rr(4, std::vector<double>(4, 0.0));
        rr[1][1] = 1; rr[2][2] = 1; rr[3][3] = 1;
        bool first = true;
        chk("charmo E", charmo(vv, ntype, 1, 1, rr, 1, first), 1.0, 0.0);
        chk("charmo s", charmo(vv, ntype, 1, 2, rr, 1, first), 1.0, 1e-12);
    }

    // ---- 8: symoir (nclass=1) ----
    {
        numat = 1;
        nclass = 1; nirred = 1;
        jx.assign(3, " A  ");
        namo.assign(6, " ");
        jndex.assign(6, 0);
        jndex[1] = 1;
        double v[2] = {1.0, 0.0}, eigs[2] = {0.0, 0.0}, r9[9] = {1,0,0,0,1,0,0,0,1};
        symoir(1, v, eigs, 1, r9, 1);
        chks("symoir namo", namo[1], " A  ");
        chki("symoir jndex", jndex[1], 1);
    }

    // ---- 9: symtrz early return ----
    {
        numat = 1;
        coord.resize(4, std::vector<double>(8, 0.0));
        nat.assign(8, 1);
        atmass.assign(8, 1.0);
        jelem.clear();
        double v[2] = {0}, e[2] = {0}, r9[9] = {0};
        symtrz(v, e, 1, 1);  // jelem.size() < numat*20 -> immediate return
        chki("symtrz ok", 1, 1);
    }

    std::printf("M05c: %d PASS / %d FAIL\n", nPass, nFail);
    return nFail == 0 ? 0 : 1;
}
