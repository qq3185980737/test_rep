// geoutg.cpp — C++ translation of MOPAC 2016 "geoutg.F90".
// Writes the geometry in Gaussian-style format.
#include "geoutg.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "chanel_C.h"
#include "elemts_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace common_arrays_C {
extern std::vector<int> na, nb, nc, labels;
extern std::vector<std::vector<double>> coord, geo;
extern std::vector<std::vector<int>> loc;
extern std::vector<std::string> simbol, txtatm;
}
namespace molkst_C {
extern int natoms, nvar, ndep, ltxt;
}
namespace symmetry_C {
extern std::vector<int> locpar, idepfn, locdep;
}
namespace chanel_C {
extern int iscr;
}
namespace elemts_C {
extern std::vector<std::string> elemnt;
}

extern void xyzint(double* xyz, int numat, int* na, int* nb, int* nc, double degree, double* geo);
extern void xxx(char type, int i, int j, int k, int l, std::string& r);
extern const char* elemnt(int);

void geoutg(int iprt) {
    static const char type_c[4] = {' ', 'r', 'a', 'd'};
    const double degree = 57.29577951308232;
    int n = molkst_C::natoms;
    // Do any atoms need Cartesian conversion?
    int i = 2;
    for (; i <= molkst_C::natoms; ++i)
        if (common_arrays_C::na[i] == 0) break;
    if (i <= molkst_C::natoms) {
        // Convert everything to Cartesian coordinates.
        for (int k = 1; k <= molkst_C::natoms; ++k) {
            common_arrays_C::coord[0][k] = common_arrays_C::geo[1][k];
            common_arrays_C::coord[1][k] = common_arrays_C::geo[2][k];
            common_arrays_C::coord[2][k] = common_arrays_C::geo[3][k];
        }
        std::vector<double> xyz(3 * molkst_C::natoms + 3, 0.0), go(3 * molkst_C::natoms + 3, 0.0);
        std::vector<int> a(molkst_C::natoms + 2, 0), b(molkst_C::natoms + 2, 0), cc(molkst_C::natoms + 2, 0);
        for (int k = 1; k <= molkst_C::natoms; ++k) {
            for (int q = 1; q <= 3; ++q) xyz[(k - 1) * 3 + (q - 1)] = common_arrays_C::coord[q-1][k];
            a[k] = common_arrays_C::na[k];
            b[k] = common_arrays_C::nb[k];
            cc[k] = common_arrays_C::nc[k];
        }
        xyzint(xyz.data(), molkst_C::natoms, a.data(), b.data(), cc.data(), 1.0, go.data());
        for (int k = 1; k <= molkst_C::natoms; ++k)
            for (int q = 1; q <= 3; ++q) common_arrays_C::geo[q][k] = go[(k - 1) * 3 + (q - 1)];
        molkst_C::nvar = 0;
        for (int k = 1; k <= molkst_C::natoms; ++k)
            for (int j = 1; j <= std::min(3, k - 1); ++j) {
                molkst_C::nvar++;
                common_arrays_C::loc[1][molkst_C::nvar] = k;
                common_arrays_C::loc[2][molkst_C::nvar] = j;
            }
    } else {
        // Constrain angles to the range 0..180, dihedrals to -180..180 degrees.
        for (int k = 4; k <= molkst_C::natoms; ++k) {
            double w = common_arrays_C::geo[2][k] * degree;
            double x = common_arrays_C::geo[3][k] * degree;
            w = w - std::floor(w / 360.0) * 360.0;
            if (w < 0) w = w + 360.0;
            if (w > 180.0) {
                x = x + 180.0;
                w = 360.0 - w;
            }
            double sgn = (x >= 0 ? 1.0 : -1.0);
            x = x - std::floor(x / 360.0 + sgn * (0.5 - 1.e-9) - 1.e-9) * 360.0;
            common_arrays_C::geo[2][k] = w / degree;
            common_arrays_C::geo[3][k] = x / degree;
        }
    }
    std::vector<std::vector<int>> igeo(4, std::vector<int>(molkst_C::natoms + 1, -1));
    for (int k = 1; k <= molkst_C::natoms; ++k) { igeo[1][k] = -1; igeo[2][k] = -1; igeo[3][k] = -1; }
    for (int k = 1; k <= molkst_C::nvar; ++k)
        igeo[common_arrays_C::loc[2][k]][common_arrays_C::loc[1][k]] = -2;
    for (int k = 1; k <= molkst_C::ndep; ++k) {
        if (symmetry_C::idepfn[k] == 14) {
            igeo[3][symmetry_C::locdep[k]] = -symmetry_C::locpar[k];
        } else {
            if (symmetry_C::idepfn[k] > 3) continue;
            igeo[symmetry_C::idepfn[k]][symmetry_C::locdep[k]] = symmetry_C::locpar[k];
        }
    }
    int maxtxt = molkst_C::ltxt;
    int nopt = 0;
    std::vector<std::vector<std::string>> line(4, std::vector<std::string>(molkst_C::natoms + 1, " "));
    std::vector<std::string> optdat(3 * molkst_C::natoms + 1, " ");
    for (int k = 1; k <= molkst_C::natoms; ++k) {
        for (int j = 1; j <= 3; ++j) {
            line[j][k] = " ";
            if (igeo[j][k] == -1) {
                char buf[16];
                if (j != 1) std::snprintf(buf, sizeof(buf), "%12.6f", common_arrays_C::geo[j][k] * degree);
                else std::snprintf(buf, sizeof(buf), "%12.6f", common_arrays_C::geo[j][k]);
                line[j][k] = buf;
            } else if (igeo[j][k] == -2) {
                nopt++;
                if (common_arrays_C::simbol[nopt] != "---------") {
                    std::string sb = common_arrays_C::simbol[nopt];
                    if (!sb.empty() && sb[0] == '-') line[j][k] = "   " + sb.substr(1);
                    else line[j][k] = "   " + sb;
                } else {
                    int nbi = common_arrays_C::nb[k];
                    int nci = common_arrays_C::nc[k];
                    if (j != 3) nci = 0;
                    if (j == 1) nbi = 0;
                    std::string r;
                    xxx(type_c[j], k, common_arrays_C::na[k], nbi, nci, r);
                    line[j][k] = "   " + r;
                }
                optdat[nopt] = line[j][k];
            } else if (igeo[j][k] < 0) {
                line[3][k] = line[3][-igeo[j][k]];
                if (line[3][k].size() >= 2) line[3][k][2] = '-';
            } else {
                line[j][k] = line[j][igeo[j][k]];
            }
        }
        std::string blank = std::string(elemnt(common_arrays_C::labels[k])) + common_arrays_C::txtatm[k] + "  ";
        if (common_arrays_C::labels[k] == 99 && !blank.empty()) blank[0] = ' ';
        int j = std::max(4, maxtxt + 2);
        std::string bl = blank.substr(0, j);
        switch (k) {
            case 1:
                std::fprintf(stdout, "%s\n", bl.c_str());
                break;
            case 2:
                std::fprintf(stdout, "%s %4d %s\n", bl.c_str(), common_arrays_C::na[k], line[1][k].substr(3).c_str());
                break;
            case 3:
                std::fprintf(stdout, "%s %4d %s %4d %s\n", bl.c_str(), common_arrays_C::na[k],
                             line[1][k].substr(3).c_str(), common_arrays_C::nb[k], line[2][k].substr(3).c_str());
                break;
            default:
                std::fprintf(stdout, "%s %4d %s %4d %s %4d %s %4d\n", bl.c_str(), common_arrays_C::na[k],
                             line[1][k].substr(3).c_str(), common_arrays_C::nb[k], line[2][k].substr(3).c_str(),
                             common_arrays_C::nc[k], line[3][k].substr(3).c_str(), 0);
                break;
        }
    }
    std::fprintf(stdout, "\n");
    for (int l = 1; l <= 3; ++l) {
        for (int k = 1; k <= nopt; ++k) {
            if (common_arrays_C::loc[2][k] != l) continue;
            if (common_arrays_C::loc[2][k] != 1)
                std::fprintf(stdout, "%s%12.6f\n", optdat[k].substr(3).c_str(),
                             common_arrays_C::geo[common_arrays_C::loc[2][k]][common_arrays_C::loc[1][k]] * degree);
            else
                std::fprintf(stdout, "%s%12.6f\n", optdat[k].substr(3).c_str(),
                             common_arrays_C::geo[common_arrays_C::loc[2][k]][common_arrays_C::loc[1][k]]);
        }
    }
}
