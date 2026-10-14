// test_dftd3.cpp — numeric checks for dftd3_bits / gdisp / dftd3 entry.
#include "dftd3_bits.h"
#include "copyc6.h"
#include "dftd3.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace common_arrays_C;
using namespace molkst_C;

int main() {
    int ok = 1;
    // lin
    if (lin(1,1)!=1 || lin(2,1)!=2 || lin(2,2)!=3 || lin(3,3)!=6) { std::printf("FAIL lin\n"); ok=0; }
    // esym
    if (esym(1)!="h " || esym(6)!="c " || esym(8)!="o " || esym(94)!="pu") { std::printf("FAIL esym\n"); ok=0; }
    // hbpar
    if (hbpar(7)!=1 || hbpar(8)!=2 || hbpar(17)!=6 || hbpar(6)!=0) { std::printf("FAIL hbpar\n"); ok=0; }
    // ncoord: single atom -> 0
    std::vector<std::vector<double>> xyz(4, std::vector<double>(2, 0.0));
    std::vector<int> nat1(2, 0);
    std::vector<double> rcov(95, 0.0);
    rcov[1] = 1.0;
    std::vector<double> cn(2, 0.0);
    nat1[1] = 1;
    ncoord(1, rcov, nat1, xyz, cn);
    if (cn[1] != 0.0) { std::printf("FAIL ncoord\n"); ok=0; }
    // getc6 at exact reference point -> exact C6 for H-H
    {
        int maxc=6, max_elem=500;
        std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> c6ab(
          max_elem+1, std::vector<std::vector<std::vector<std::vector<double>>>>(max_elem+1,
          std::vector<std::vector<std::vector<double>>>(maxc+1,
          std::vector<std::vector<double>>(maxc+1, std::vector<double>(4, 0.0)))));
        std::vector<int> maxci(max_elem+1, 0);
        copyc6(maxc, max_elem, c6ab, maxci);
        double c6 = 0.0;
        getc6(maxc, max_elem, c6ab, maxci, 1, 1, 0.9118, 0.9118, c6);
        if (std::fabs(c6 - 3.0267) > 1e-6) { std::printf("FAIL getc6 H-H: %g\n", c6); ok=0; }
        // C-H pair: row (6,1) ref (0.0, 0.9118) C6=12.1402
        getc6(maxc, max_elem, c6ab, maxci, 6, 1, 0.0, 0.9118, c6);
        if (std::fabs(c6 - 12.1402) > 1e-4) { std::printf("FAIL getc6 C-H: %g\n", c6); ok=0; }
    }
    // eabh: linear A-H-B, energy must be negative and finite
    {
        std::vector<std::vector<double>> xyz2(4, std::vector<double>(5, 0.0));
        xyz2[1][1] = -2.0; xyz2[1][2] = 2.0; xyz2[1][3] = 0.0;
        double e = eabh(3, 1, 2, 3, xyz2, 4.0, 0.5);
        if (!(e < 0.0) || !std::isfinite(e)) { std::printf("FAIL eabh sign: %g\n", e); ok=0; }
    }
    // edisp: isolated atoms -> zero
    {
        std::vector<std::vector<double>> xyz3(4, std::vector<double>(3, 0.0));
        std::vector<int> nat3(3, 0);
        std::vector<double> rcov3(95, 0.0), r2r4(95, 1.0);
        std::vector<std::vector<double>> r0ab3(95, std::vector<double>(95, 1.0));
        std::vector<int> mxc3(95, 1);
        int maxc=1, max_elem=94;
        std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> c6ab3(
          max_elem+1, std::vector<std::vector<std::vector<std::vector<double>>>>(max_elem+1,
          std::vector<std::vector<std::vector<double>>>(maxc+1,
          std::vector<std::vector<double>>(maxc+1, std::vector<double>(4, 0.0)))));
        nat3[1]=1; nat3[2]=1;
        xyz3[1][2] = 100.0;  // far apart
        double e6=0, e8=0;
        edisp(max_elem, maxc, 2, xyz3, nat3, c6ab3, mxc3, r2r4, r0ab3, rcov3, 1.56, 1.0, 14.0, 16.0, e6, e8);
        if (e6 != 0.0 || e8 != 0.0) { std::printf("FAIL edisp far: %g %g\n", e6, e8); ok=0; }
    }
    // end-to-end dftd3 on H2 (0.74 A): must run, return finite <= 0
    {
        numcal = 1;
        keywrd = "";
        numat = 2;
        nat.assign(3, 0); nat[1]=1; nat[2]=1;
        coord.assign(4, std::vector<double>(3, 0.0));
        coord[1][2] = 0.0; coord[2][2] = 0.0; coord[3][2] = 0.74;
        std::vector<std::vector<double>> dxyz(4, std::vector<double>(3, 0.0));
        double e = dftd3(false, dxyz);
        if (!std::isfinite(e) || e > 0.0) { std::printf("FAIL dftd3 H2: %g\n", e); ok=0; }
        std::printf("dftd3(H2) = %12.6g kcal/mol\n", e);
    }
    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}