// test_pinout.cpp — round-trip write/read of pinout (mode=1 then mode=2).
#include "pinout.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "chanel_C.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace molkst_C {
extern int numat, nelecs, norbs;
extern std::string keywrd, line;
}
namespace chanel_C {
extern int iw;
extern std::string density_fn;
}
namespace common_arrays_C {
extern std::vector<int> nbonds;
extern std::vector<std::vector<int>> ibonds;
}
extern std::vector<int> nce, ncf, ncvir, nncf, nnce, ncocc, icocc, icvir;
extern std::vector<double> cocc, cvir;
extern int cocc_dim, cvir_dim, icocc_dim, icvir_dim;

int main() {
    using MOZYME_C::iorbs;
    bool ok = true;
    molkst_C::numat = 2;
    molkst_C::nelecs = 2;
    molkst_C::norbs = 2;
    molkst_C::keywrd = "";
    chanel_C::density_fn = "test_pinout.density";

    ncf.assign(2, 0); ncf[1] = 1;
    nce.assign(2, 0); nce[1] = 1;
    nncf.assign(2, 0); nncf[1] = 0;
    nnce.assign(2, 0); nnce[1] = 0;
    ncocc.assign(2, 0); ncocc[1] = 0;
    ncvir.assign(2, 0); ncvir[1] = 0;
    icocc.assign(2, 1); icocc[1] = 1;
    icvir.assign(2, 1); icvir[1] = 2;
    cocc.assign(4, 0.0); cocc[1] = 0.12345;
    cvir.assign(4, 0.0); cvir[1] = 0.67890;
    icocc_dim = 2; cocc_dim = 3; icvir_dim = 2; cvir_dim = 3;
    iorbs.assign(3, 1); iorbs[1] = 1; iorbs[2] = 1;
    common_arrays_C::nbonds.assign(3, 0);
    common_arrays_C::nbonds[1] = 1; common_arrays_C::nbonds[2] = 1;
    common_arrays_C::ibonds.assign(10, std::vector<int>(3, 0));
    common_arrays_C::ibonds[1][1] = 2; common_arrays_C::ibonds[1][2] = 1;

    pinout(1);

    // Reset incoming arrays, then read back.
    cocc.assign(4, -1.0); cvir.assign(4, -1.0);
    icocc.assign(2, 0); icvir.assign(2, 0);
    ncf.assign(2, 0); nce.assign(2, 0);
    iorbs.assign(3, 0);
    common_arrays_C::nbonds.assign(3, 0);
    pinout(2);

    if (ncf[1] != 1) { std::printf("FAIL ncf\n"); ok = false; }
    if (nce[1] != 1) { std::printf("FAIL nce\n"); ok = false; }
    if (icocc[1] != 1 || icvir[1] != 2) { std::printf("FAIL icocc/icvir\n"); ok = false; }
    if (std::fabs(cocc[1] - 0.12345) > 1e-9) { std::printf("FAIL cocc=%g\n", cocc[1]); ok = false; }
    if (std::fabs(cvir[1] - 0.67890) > 1e-9) { std::printf("FAIL cvir=%g\n", cvir[1]); ok = false; }
    if (iorbs[1] != 1 || iorbs[2] != 1) { std::printf("FAIL iorbs\n"); ok = false; }
    if (common_arrays_C::nbonds[1] != 1 || common_arrays_C::ibonds[1][2] != 1) {
        std::printf("FAIL bonds\n"); ok = false;
    }

    // Missing-file path (mode 0).
    chanel_C::density_fn = "test_pinout.nonexistent";
    pinout(0);

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
