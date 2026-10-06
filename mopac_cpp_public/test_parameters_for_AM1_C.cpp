// test_parameters_for_AM1_C.cpp — spot-check AM1 parameter data (Fortran data statements).
#include "parameters_for_AM1_C.h"

#include <cstdio>

using namespace parameters_for_AM1_C;

int main() {
    bool ok = true;
    // H: uss = -11.3964270, betas = -6.1737870, zs = 1.1880780, alp = 2.8823240,
    //     gss = 12.8480000, polvol = 0.1719890
    auto near = [](double a, double b) { return (a - b) * (a - b) < 1e-8; };
    if (!near(ussam1[1], -11.3964270)) { std::printf("FAIL ussam1[1]=%.10f\n", ussam1[1]); ok = false; }
    if (!near(betasa[1], -6.1737870)) { std::printf("FAIL betasa[1]\n"); ok = false; }
    if (!near(zsam1[1], 1.1880780)) { std::printf("FAIL zsam1[1]\n"); ok = false; }
    if (!near(alpam1[1], 2.8823240)) { std::printf("FAIL alpam1[1]\n"); ok = false; }
    if (!near(gssam1[1], 12.8480000)) { std::printf("FAIL gssam1[1]\n"); ok = false; }
    if (!near(polvolam1[1], 0.1719890)) { std::printf("FAIL polvolam1[1]\n"); ok = false; }
    // He: uss = -35.2271054, upp = 9.9998070
    if (!near(ussam1[2], -35.2271054)) { std::printf("FAIL ussam1[2]\n"); ok = false; }
    if (!near(uppam1[2], 9.9998070)) { std::printf("FAIL uppam1[2]\n"); ok = false; }
    // guesa1(1,1)=0.1227960, guesa1(1,3)=-0.0183360; guesa3(1,1)=1.2, (1,3)=2.1
    if (!near(guesa1[1][1], 0.1227960)) { std::printf("FAIL guesa1[1][1]\n"); ok = false; }
    if (!near(guesa1[1][3], -0.0183360)) { std::printf("FAIL guesa1[1][3]\n"); ok = false; }
    if (!near(guesa3[1][1], 1.2000000)) { std::printf("FAIL guesa3[1][1]\n"); ok = false; }
    if (!near(guesa3[1][3], 2.1000000)) { std::printf("FAIL guesa3[1][3]\n"); ok = false; }
    // Li: uss = -4.9384384
    if (!near(ussam1[3], -4.9384384)) { std::printf("FAIL ussam1[3]\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
