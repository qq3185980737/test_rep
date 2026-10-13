// test_mozyme_interfaces.cpp — interface-level unit tests for the MOZYME
// SCF driver dependencies of iter_for_MOZYME.cpp.  Every expected value is
// hand-computed from the Fortran original (MOPAC2016).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "reada.h"
#include "cnvgz.h"
#include "isitsc.h"
#include "scfcri.h"
#include "eimp.h"
#include "check.h"
#include "density_for_MOZYME.h"
#include "tidy.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "overlaps_C.h"
#include "funcon_C.h"
#include "iter_C.h"
#include "cosmo_C.h"

using namespace molkst_C;
using namespace MOZYME_C;
using namespace common_arrays_C;
using namespace chanel_C;
using namespace overlaps_C;
using namespace funcon_C;
using namespace iter_C;
using namespace cosmo_C;

static int npass = 0, nfail = 0;
#define CHECK(name, cond)                              \
    do {                                               \
        bool _ok = (cond);                             \
        std::printf("[%s] %s\n", _ok ? "PASS" : "FAIL", name); \
        if (_ok) ++npass; else ++nfail;                \
    } while (0)
#define CHECK_EQ(name, a, b, tol)                      \
    do {                                               \
        double _a = (a), _b = (b);                     \
        bool _ok = std::fabs(_a - _b) < (tol);         \
        std::printf("[%s] %s got=%.9g want=%.9g\n", _ok ? "PASS" : "FAIL", \
                    name, _a, _b);                     \
        if (_ok) ++npass; else ++nfail;                \
    } while (0)

int main() {
    ws.resize(16, 0.0);   // check() writes ws[1..nvec]
    // ============ reada (Fortran reada.F90) ============
    CHECK_EQ("reada: 'PM7 1SCF' @5", reada("PM7 1SCF", 5), 1.0, 1e-12);
    CHECK_EQ("reada: THRESH=1.D-13 @8", reada("THRESH=1.D-13", 8), 1e-13, 1e-25);
    CHECK_EQ("reada: RELSCF=0.5 @2", reada(" RELSCF=0.5", 2), 0.5, 1e-12);
    CHECK_EQ("reada: no number -> 0", reada("abc", 1), 0.0, 1e-12);
    CHECK_EQ("reada: -0.25D+01", reada("-0.25D+01", 1), -2.5, 1e-12);
    CHECK_EQ("reada: 1.5e2", reada("1.5e2", 1), 150.0, 1e-9);
    CHECK_EQ("reada: 1.0E-03 @2", reada(" 1.0E-03", 2), 1e-3, 1e-15);
    CHECK_EQ("reada: 123 trailing x", reada("123xyz", 1), 123.0, 1e-12);

    // ============ cnvgz (Fortran cnvgz.F90) ============
    // NOTE: local names p/p1/p2/p3 shadow the globals; that is intended here
    // (cnvgz takes all its arrays as parameters).
    norbs = 2; mpack = 3;
    std::vector<int> idiag = {0, 1, 2};            // idiag[1]=1, idiag[2]=2
    use_three_point_extrap = true; pmax = 0.0;
    std::vector<double> pnew = {0, 1.0, 2.0, 3.0}, p = {0, 0.9, 2.1, 3.0},
                        p1 = {0, 0, 0, 0}, p2 = {0, 0, 0, 0}, p3 = {0, 0, 0, 0};
    cnvgz(pnew, p, p1, p2, p3, 3, idiag);
    // p3[i]=pnew[idiag[i]]: p3=[0,1,2,0]; p2=[0,0.9,2.1,0]; pmax=max(|1-0.9|,|2-2.1|,|3-3|)=0.1
    // faca=|1-.9|^2+|2-2.1|^2=0.02; facb=(1-1.8)^2+(2-4.2)^2=0.64+4.84=5.48
    // fac=sqrt(0.02/5.48); pnew[i]+=fac*(pnew[i]-p[i])
    double fac = std::sqrt(0.02 / 5.48);
    CHECK_EQ("cnvgz n=3: pnew[1]", pnew[1], 1.0 + fac * 0.1, 1e-12);
    CHECK_EQ("cnvgz n=3: pnew[2]", pnew[2], 2.0 - fac * 0.1, 1e-12);
    CHECK_EQ("cnvgz n=3: pnew[3]", pnew[3], 3.0, 1e-12);
    CHECK_EQ("cnvgz n=3: pmax", pmax, 0.1, 1e-12);
    CHECK_EQ("cnvgz n=3: p1==p2", p1[1], 0.9, 1e-12);
    CHECK_EQ("cnvgz n=3: p==pnew", p[1], pnew[1], 1e-12);

    // damp branch: niter=4, pmax=0.1>0.05, diag diffs 0.1,-0.1 both >0.05
    // Fortran: pnew(j)=p2(i)+Sign(damp,p3(i)-p2(i)); then pnew(j)=Min(2.0,Max(pnew(j),0.0))
    // j=1: 0.9+0.05=0.95 (kept); j=2: 2.1-0.05=2.05 -> clamped to 2.0
    std::vector<double> qnew = {0, 1.0, 2.0, 3.0}, q = {0, 0.9, 2.1, 3.0};
    p1 = {0, 0, 0, 0}; p2 = {0, 0, 0, 0}; p3 = {0, 0, 0, 0};
    use_three_point_extrap = true;
    cnvgz(qnew, q, p1, p2, p3, 4, idiag);
    CHECK_EQ("cnvgz damp: qnew[1]", qnew[1], 0.95, 1e-12);
    CHECK_EQ("cnvgz damp: qnew[2] clamped", qnew[2], 2.0, 1e-12);

    // ============ isitsc (Fortran isitsc.F90) ============
    ovmax = 1e-10; energy_diff = 1e-10; iscf = 0;
    int iemin = 0, iemax = 0; bool okscf = false;
    isitsc(-60.0, 1e-2, -50.0, iemin, iemax, okscf, 2, 2000);
    CHECK("isitsc: first call NOT converged", !okscf);
    isitsc(-61.0, 1e-2, -50.0, iemin, iemax, okscf, 3, 2000);
    CHECK("isitsc: second call converged via scf1", okscf);
    CHECK("isitsc: iscf==1 on rapid exit", iscf == 1);

    // ============ scfcri (Fortran scfcri.F90) ============
    numcal = 1; keywrd = "PM7 1SCF XYZ"; efield[1] = efield[2] = efield[3] = 0.0;
    double selcon = 0.0;
    scfcri(selcon);
    CHECK_EQ("scfcri: default 1e-2", selcon, 1e-2, 1e-15);
    scfcri(selcon);   // same numcal -> else branch, no change
    CHECK_EQ("scfcri: re-call unchanged", selcon, 1e-2, 1e-15);
    numcal = 2; keywrd = "PM7 RELSCF=0.5";
    scfcri(selcon);
    CHECK_EQ("scfcri: RELSCF=0.5 -> 5e-3", selcon, 5e-3, 1e-15);
    numcal = 3; keywrd = "PM7 SCFCRT=1.D-4";
    scfcri(selcon);
    CHECK_EQ("scfcri: SCFCRT=1e-4", selcon, 1e-4, 1e-18);

    // ============ eimp (Fortran eimp.F90) ============
    numat = 2; iorbs = {0, 1, 1};
    lijbo = true;
    nijbo.assign(3, std::vector<int>(3, 0));
    f = {0, 0.5, 0.6, 0.7, 0.8, 0.9};          // 1-based, size>=2
    common_arrays_C::p = {0, 0.0, 0.0};       // NOTE: local 'p' shadows global p; use qualified name
    eimp();
    // i=2,j=1: k=ijbo(2,1)=0; l=iorbs[2]*iorbs[1]=1; sum=f(1)^2=0.25; p(1)=sum
    // NOTE: local 'p' from the cnvgz test shadows the global; check the global.
    CHECK_EQ("eimp: p[1]=f[1]^2", common_arrays_C::p[1], 0.25, 1e-12);

    // ============ check (Fortran check.F90) ============
    std::vector<int> nnc = {0, 0}, nc = {0, 2}, icvec = {0, 1, 2},
                     ncvec = {0, 0};
    std::vector<double> cvec = {0, 0.5, 0.8};
    check(1, nnc, nc, icvec, iorbs, ncvec, cvec);
    // sum=0.25+0.64=0.89, scale=1/sqrt(0.89); cvec[i]*=scale
    double scale = 1.0 / std::sqrt(0.89);
    CHECK_EQ("check: cvec[1] renormalized", cvec[1], 0.5 * scale, 1e-12);
    CHECK_EQ("check: cvec[2] renormalized", cvec[2], 0.8 * scale, 1e-12);
    CHECK_EQ("check: ws[1]=sum", ws[1], 0.89, 1e-12);

    // ============ density_for_MOZYME (Fortran density_for_MOZYME.F90) ============
    // One occupied LMO over 2 atoms (one s orbital each): cocc=(0.6, 0.8).
    // Packed lower triangle addresses (addr-1): (1,1)->0, (2,1)->1, (2,2)->2.
    keywrd = "PM7";                       // no " DIAG " -> prnt=false
    mpack = 3;
    std::vector<double> pp = {0, 0.0, 0.0, 0.0};
    std::vector<double> partp0 = {0, 0.0, 0.0, 0.0};
    ncocc = {0, 0}; nncf = {0, 0}; ncf = {0, 2};
    icocc = {0, 1, 2};
    cocc = {0, 0.6, 0.8};
    lijbo = true;
    nijbo.assign(3, std::vector<int>(3, 0));
    nijbo[1][1] = 0; nijbo[2][1] = 1; nijbo[2][2] = 2;
    density_for_MOZYME(pp, 0, 1, partp0);
    // P = c c^T: P(1,1)=0.36 P(2,1)=0.48 P(2,2)=0.64; x2 for spin:
    // p[1]=0.72, p[2]=0.96, p[3]=1.28
    CHECK_EQ("density_for_MOZYME: p[1]", pp[1], 0.72, 1e-12);
    CHECK_EQ("density_for_MOZYME: p[2]", pp[2], 0.96, 1e-12);
    CHECK_EQ("density_for_MOZYME: p[3]", pp[3], 1.28, 1e-12);

    // ============ tidy (Fortran tidy.F90) ============
    // Fortran "integer :: nmol = 0" is implicitly SAVE: after a mode==2 call
    // nmol==numcal and subsequent calls skip the "new molecule" branch, so
    // imode[] is NOT reset and selmos() is reached when step_num changes.
    // A plain local nmol (old C++) re-entered that branch every call, which
    // reset imode each time and suppressed selmos().
    numat = 2; norbs = 2; numred = 0; thresh = 1e-6;
    numcal = 1; step_num = 2; keywrd = "PM7"; moperr = false;
    int nmos_loc = 2, n01 = 200, n02 = 200, ln = 0, mn = 0;
    std::vector<int> nc2 = {0, 2, 2}, ic2(201, 0), nnc_loc2 = {0, 0, 2},
                     ncmo2 = {0, 0, 2};
    std::vector<double> c2(201, 0.0);
    ic2[1] = 1; ic2[2] = 2; ic2[3] = 1; ic2[4] = 2;
    c2[1] = 0.5; c2[2] = 0.6; c2[3] = 0.7; c2[4] = 0.8;
    // call 1: mode=2, new molecule (nmol=0!=numcal); compress + reallocate
    tidy(nmos_loc, nc2.data(), ic2.data(), n01, c2.data(), n02,
         nnc_loc2.data(), ncmo2.data(), ln, mn, 2);
    // after compression ln=mn=4; free space (200-4)/2=96>numat -> ispace=jspace=0;
    // reorder moves LMO2 to the very end: nnc_loc2[2]=198, ncmo2[2]=198,
    // ic2[199]=1, ic2[200]=2, c2[199]=0.7, c2[200]=0.8
    CHECK("tidy c1: ln==4", ln == 4);
    CHECK("tidy c1: mn==4", mn == 4);
    CHECK("tidy c1: nnc_loc2[2]==198", nnc_loc2[2] == 198);
    CHECK("tidy c1: ncmo2[2]==198", ncmo2[2] == 198);
    CHECK("tidy c1: ic2[199]==1", ic2[199] == 1);
    CHECK("tidy c1: ic2[200]==2", ic2[200] == 2);
    CHECK_EQ("tidy c1: c2[199]", c2[199], 0.7, 1e-12);
    CHECK_EQ("tidy c1: c2[200]", c2[200], 0.8, 1e-12);
    // call 2: mode=1, step_num advanced to 3; with the fix nmol==numcal ->
    // skip "new molecule" -> imode[1] still 2 (set by call 1) ->
    // step_num(3)!=imode[1] -> isnew=2 -> selmos() runs (numred=0 so all
    // LMOs go to the "not selected" tail list); the big loop then iterates
    // again (isnew=2 != 0 -> exit) and re-compresses, so the final state is
    // the compressed form nnc_loc2=[0,0,2], ncmo2=[0,0,2] with the selmos
    // move of coefficients left at c2[197..200]={.5,.6,.7,.8}.
    // Without the fix nmol==0 resets imode[1]=3 -> step_num==imode exit ->
    // no selmos -> the plain reorder puts LMO2 at the very end:
    // nnc_loc2=[0,0,198], ncmo2=[0,2,198], c2[197..198]=0.
    step_num = 3;
    tidy(nmos_loc, nc2.data(), ic2.data(), n01, c2.data(), n02,
         nnc_loc2.data(), ncmo2.data(), ln, mn, 1);
    // Final state after call 2: the reorder block (after the big loop)
    // runs unconditionally, so nnc_loc2[2]/ncmo2[2] are 198 in BOTH the
    // fixed and broken builds.  The distinguishing trace of selmos() is
    // its coefficient move: c2[197..200]={.5,.6,.7,.8} survive the second
    // compression (which re-reads them into c2[1..4]) -- without selmos
    // those cells stay 0.
    CHECK("tidy c2: reorder anchor nnc_loc2[2]==198", nnc_loc2[2] == 198);
    CHECK("tidy c2: reorder anchor ncmo2[2]==198", ncmo2[2] == 198);
    CHECK_EQ("tidy c2: c2[197] moved by selmos", c2[197], 0.5, 1e-12);
    CHECK_EQ("tidy c2: c2[198] moved by selmos", c2[198], 0.6, 1e-12);

    std::printf("\n%d passed, %d failed\n", npass, nfail);
    return nfail == 0 ? 0 : 1;
}
