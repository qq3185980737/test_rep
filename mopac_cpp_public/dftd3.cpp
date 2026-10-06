// dftd3.cpp — C++ translation of MOPAC 2016 "dftd3.F90".
// Grimme D3 dispersion: edisp/gdisp/hbsimple/setr0ab/copyc6 from dftd3_bits/copyc6.

#include "dftd3.h"
#include "dftd3_bits.h"
#include "copyc6.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace funcon_C;
using namespace molkst_C;

namespace {
const int max_elem = 94, maxc = 5;
std::vector<double> r2r4v = {0,8.0589,3.4698,29.0974,14.8517,11.8799,7.8715,5.5588,4.7566,3.8025,3.1036,
26.1552,17.2304,17.7210,12.7442,9.5361,8.1652,6.7463,5.6004,29.2012,22.3934,
19.0598,16.8590,15.4023,12.5589,13.4788,12.2309,11.2809,10.5569,10.1428,9.4907,
13.4606,10.8544,8.9386,8.1350,7.1251,6.1971,30.0162,24.4103,20.3537,17.4780,
13.5528,11.8451,11.0355,10.1997,9.5414,9.0061,8.6417,8.9975,14.0834,11.8333,
10.0179,9.3844,8.4110,7.5152,32.7622,27.5708,23.1671,21.6003,20.9615,20.4562,
20.1010,19.7475,19.4828,15.6013,19.2362,17.4717,17.8321,17.4237,17.1954,17.1631,
14.5716,15.8758,13.8989,12.4834,11.4421,10.2671,8.3549,7.8496,7.3278,7.4820,
13.5124,11.6554,10.0959,9.7340,8.8584,8.0125,29.8135,26.3157,19.1885,15.8542,
16.1305,15.6161,15.1226,16.1576};
std::vector<double> rcov = {0,0.32,0.46,1.20,0.94,0.77,0.75,0.71,0.63,0.64,0.67,
1.40,1.25,1.13,1.04,1.10,1.02,0.99,0.96,1.76,1.54,
1.33,1.22,1.21,1.10,1.07,1.04,1.00,0.99,1.01,1.09,
1.12,1.09,1.15,1.10,1.14,1.17,1.89,1.67,1.47,1.39,
1.32,1.24,1.15,1.13,1.13,1.08,1.15,1.23,1.28,1.26,
1.26,1.23,1.32,1.31,2.09,1.76,1.62,1.47,1.58,1.57,
1.56,1.55,1.51,1.52,1.51,1.50,1.49,1.49,1.48,1.53,
1.46,1.37,1.31,1.23,1.18,1.16,1.11,1.12,1.13,1.32,
1.30,1.30,1.36,1.31,1.38,1.42,2.01,1.81,1.67,1.58,
1.52,1.53,1.54,1.55};

// module-persistent D3 state (Fortran save)
int icalcn_d3 = -1;
c6ab_t c6ab(max_elem+1, std::vector<std::vector<std::vector<std::vector<double>>>>(max_elem+1,
    std::vector<std::vector<std::vector<double>>>(maxc+1,
    std::vector<std::vector<double>>(maxc+1, std::vector<double>(4, 0.0)))));
std::vector<int> mxc(max_elem+1, 0);
std::vector<std::vector<double>> r0ab(max_elem+1, std::vector<double>(max_elem+1, 0.0));
}

double dftd3(bool l_grad, std::vector<std::vector<double>>& dxyz) {
    static double s6, alp6, rs6, s18, alp8, rs8, hbscale, au_to_kcal;
    if (icalcn_d3 != numcal) {
        icalcn_d3 = numcal;
        for (int i = 1; i <= max_elem; ++i) rcov[i] = 4.0 / 3.0 * rcov[i] / a0;
        au_to_kcal = fpc_9 * fpc_2 / a0;
        for (int i = 1; i <= max_elem; ++i)
            r2r4v[i] = std::sqrt(0.5 * r2r4v[i] * std::sqrt((double)i));
        setr0ab(max_elem, a0, r0ab);
        copyc6(maxc, max_elem, c6ab, mxc);
        bool D3H4 = (keywrd.find("D3H4") != std::string::npos ||
                     keywrd.find("D3(H4)") != std::string::npos);
        if (D3H4) { s6=0.88; alp6=22.0; rs6=1.18; s18=0.0; }
        else      { s6=1.0;  alp6=14.0; rs6=1.560; s18=1.009; }
        hbscale = 1.301;
        rs8 = 1.0; alp8 = alp6 + 2.0;
    }
    // coordinates: Angstrom -> au
    std::vector<std::vector<double>> xyz(4, std::vector<double>(numat + 1, 0.0));
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= numat; ++j) xyz[i][j] = coord[i-1][j] / a0;
    double e6 = 0.0, e8 = 0.0;
    edisp(max_elem, maxc, numat, xyz, nat, c6ab, mxc, r2r4v, r0ab, rcov,
          rs6, rs8, alp6, alp8, e6, e8);
    e6 *= s6; e8 *= s18;
    E_disp = (-e6 - e8) * au_to_kcal;
    double ehb = 0.0;
    std::vector<std::vector<double>> dxyz_temp(4, std::vector<double>(numat + 1, 0.0));
    bool D3H4 = (keywrd.find("D3H4") != std::string::npos ||
                 keywrd.find("D3(H4)") != std::string::npos);
    if (!D3H4)
        hbsimple(numat, nat, xyz, hbscale, ehb, l_grad, dxyz_temp);
    E_hb = ehb * au_to_kcal;
    double dftd3 = E_disp + E_hb;
    if (l_grad) {
        gdisp(xyz, r0ab, rs6, alp6, c6ab, s6, mxc, rcov, dxyz_temp);
        for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= numat; ++j)
                dxyz[i][j] += 2.0 * dxyz_temp[i][j] * au_to_kcal;
        if (keywrd.find(" DERIV") != std::string::npos) {
            std::printf("                  GRIMME'S D3 CORRECTIONS\n");
            std::printf("   NUMBER  ATOM        X             Y             Z\n");
            for (int j = 1; j <= numat; ++j)
                std::printf("%6d%4s%13.6f%13.6f%13.6f\n", j,
                    esym(nat[j]).c_str(),
                    dxyz_temp[1][j] * 2.0 * au_to_kcal,
                    dxyz_temp[2][j] * 2.0 * au_to_kcal,
                    dxyz_temp[3][j] * 2.0 * au_to_kcal);
        }
    }
    return dftd3;
}