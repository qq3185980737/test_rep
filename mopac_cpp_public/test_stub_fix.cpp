// test_stub_fix.cpp — verify the stub-completion translations:
//   setcup       (moldat_helpers.cpp; moldat.F90 897-1035)
//   setup_nhco   (moldat_helpers.cpp; moldat.F90 1250-1323)
//   delete_MOZYME_arrays (set_up_MOZYME_arrays.cpp)
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molmec_C.h"
#include "molkst_C.h"
#include "moldat_helpers.h"
#include "reada.h"
#include "set_up_MOZYME_arrays.h"

// fillij.cpp calls the Fortran-symbol reada_; forward it to the C++ version.
extern "C" double reada_(const char* s, int* i, int len) {
    return reada(std::string(s, s + len), *i);
}

static int npass = 0, nfail = 0;
static void chk(bool ok, const char* what) {
    if (ok) { ++npass; std::printf("PASS %s\n", what); }
    else    { ++nfail; std::printf("FAIL %s\n", what); }
}
static bool near(double a, double b) { return std::fabs(a - b) < 1.0e-9; }

int main() {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace molmec_C;

    // ---- setcup: id = 0 (isolated molecule) — quick return ----
    id = 0; keywrd = " PM7";
    tvec.assign(4, std::vector<double>(4, 0.0));
    setcup();
    chk(near(cutofp, 1.0e10), "setcup id=0 cutofp=1e10");
    chk(l1u == 0 && l2u == 0 && l3u == 0, "setcup id=0 l1u/l2u/l3u=0");
    chk(l123 == 1, "setcup id=0 l123=1");

    // ---- setcup: id = 3 (3D crystal), no CUTOFP -> default 30 ----
    id = 3; keywrd = " PM7";
    // Orthogonal 10-Angstrom cell.
    tvec[1][1] = 10.0; tvec[2][1] = 0.0;  tvec[3][1] = 0.0;
    tvec[1][2] = 0.0;  tvec[2][2] = 10.0; tvec[3][2] = 0.0;
    tvec[1][3] = 0.0;  tvec[2][3] = 0.0;  tvec[3][3] = 10.0;
    setcup();
    chk(near(cutofp, 30.0), "setcup id=3 default cutofp=30");
    // Orthogonal cell: vol=1000, area23=100*sin(acos(0.5))=86.6025,
    // tv1=vol/area23=11.5470 -> l1u = int(30*4/3/11.5470)+1 = int(3.4641)+1 = 4
    // Orthogonal cell: r12=r13=r23=sqrt(200)=14.142 (body diagonals),
    // area12=area13=area23=100*sin(acos(0))=100 -> tv1=vol/area=10
    // -> l1u = int(30*4/3/10)+1 = int(4)+1 = 5; l123 = 11^3 = 1331
    chk(l1u == 5 && l2u == 5 && l3u == 5, "setcup id=3 l1u=l2u=l3u=5");
    chk(l123 == 1331, "setcup id=3 l123=1331");
    chk(l11 == 1 && l21 == 1 && l31 == 1, "setcup id=3 l11/l21/l31=1");

    // ---- setcup: id = 3 with CUTOFP ----
    keywrd = " PM7 CUTOFP=10.0";
    setcup();
    chk(near(cutofp, 10.0), "setcup CUTOFP=10.0 parsed");
    chk(l1u == 2, "setcup CUTOFP l1u=2 (int(10*4/3/11.547)+1)");

    // ---- setup_nhco: no O=C-N-H system ----
    keywrd = " PM7";
    method_pm7 = true; method_am1 = false; method_pm3 = false;
    method_mndo = false; method_rm1 = false; method_pm6 = false;
    nat.assign(3, 0);
    coord.assign(4, std::vector<double>(3, 0.0));  // coord[1..3][1..2]
    numat = 2; nat[1] = 6; nat[2] = 1;   // CH2-like, no C=O/N/H chain
    int ii = 0;
    setup_nhco(ii);
    chk(nnhco == 0, "setup_nhco no-system nnhco=0");
    chk(near(htype, 3.1595), "setup_nhco PM7 htype=3.1595");
    chk(ii == 0, "setup_nhco ii=0 without NOMM");

    // NOMM keyword -> ii=1
    keywrd = " PM7 NOMM";
    setup_nhco(ii);
    chk(ii == 1, "setup_nhco ii=1 with NOMM");

    // ---- delete_MOZYME_arrays: must not crash and clears state ----
    MOZYME_C::iorbs.assign(3, 0);
    MOZYME_C::nijbo.assign(3, std::vector<int>(3, 0));
    numat = 2;
    delete_MOZYME_arrays();
    chk(MOZYME_C::iorbs.empty(), "delete_MOZYME_arrays iorbs cleared");
    chk(MOZYME_C::nijbo.empty(), "delete_MOZYME_arrays nijbo cleared");

    std::printf("\n%d passed, %d failed\n", npass, nfail);
    return nfail == 0 ? 0 : 1;
}
