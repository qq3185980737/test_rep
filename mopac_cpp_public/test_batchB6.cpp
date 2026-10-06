// test_batchB6.cpp — tests: getgeo (MOPAC geometry reader, internal coords).
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "getgeo.h"

namespace molkst_C {
extern int natoms, numat, maxtxt, numcal, id;
extern bool moperr, units, Angstroms;
extern std::string keywrd, line;
extern double arc_hof_1, arc_hof_2;
}
namespace common_arrays_C {
extern std::vector<double> atmass;
extern std::vector<std::string> simbol, txtatm;
extern std::vector<int> na_store, nat;
extern std::vector<char> l_atom;
}
namespace parameters_C {
extern std::vector<double> ams;
}
int main() {
    bool ok = true;
    molkst_C::keywrd = "AM1";
    molkst_C::natoms = 12;
    molkst_C::numat = 0;
    molkst_C::numcal = 1;
    molkst_C::moperr = false;
    common_arrays_C::simbol.assign(40, "---------");
    common_arrays_C::txtatm.assign(12, " ");
    common_arrays_C::atmass.assign(12, 0.0);
    common_arrays_C::nat.assign(12, 0);
    common_arrays_C::na_store.assign(12, 0);
    common_arrays_C::l_atom.assign(12, false);
    parameters_C::ams.resize(108, 0.0);
    parameters_C::ams[8] = 15.9994; parameters_C::ams[1] = 1.0079; parameters_C::ams[6] = 12.011;

    std::string in =
        "  O   0.0      0.0      0.0\n"
        "  H   0.957    0        0.0     0    0.0    0\n"
        "  H   0.957    0        104.5   1    0.0    0    1    2    0\n";
    std::istringstream is(in);
    std::cin.rdbuf(is.rdbuf());

    std::vector<int> labels(12, 0);
    std::vector<std::vector<double>> geo(4, std::vector<double>(12, 0.0));
    std::vector<std::vector<double>> xyz(4, std::vector<double>(12, 0.0));
    std::vector<std::vector<int>> lopt(4, std::vector<int>(12, 0));
    std::vector<int> na(12, 0), nb(12, 0), nc(12, 0);
    getgeo(5, labels, geo, xyz, lopt, na, nb, nc, false);

    if (molkst_C::natoms != 3) { std::fprintf(stderr, "FAIL natoms=%d\n", molkst_C::natoms); ok = false; }
    if (molkst_C::numat != 3) { std::fprintf(stderr, "FAIL numat=%d\n", molkst_C::numat); ok = false; }
    if (labels[1] != 8 || labels[2] != 1 || labels[3] != 1) {
        std::fprintf(stderr, "FAIL labels %d %d %d\n", labels[1], labels[2], labels[3]); ok = false;
    }
    if (std::fabs(geo[1][1]) > 1e-4 || std::fabs(geo[2][1]) > 1e-4 || std::fabs(geo[3][1]) > 1e-4) {
        std::fprintf(stderr, "FAIL geo1 %f %f %f\n", geo[1][1], geo[2][1], geo[3][1]); ok = false;
    }
    if (std::fabs(geo[1][2] - 0.957) > 1e-3) { std::fprintf(stderr, "FAIL bond=%f\n", geo[1][2]); ok = false; }
    double ang = geo[2][3] / 1.7453292519943e-2;
    if (std::fabs(ang - 104.5) > 1e-3) { std::fprintf(stderr, "FAIL angle=%f\n", ang); ok = false; }
    if (na[3] != 1 || nb[3] != 2 || nc[3] != 0) {
        std::fprintf(stderr, "FAIL conn %d %d %d\n", na[3], nb[3], nc[3]); ok = false;
    }
    if (std::fabs(common_arrays_C::atmass[1] - 15.9994) > 1e-3) {
        std::fprintf(stderr, "FAIL atmass1=%f\n", common_arrays_C::atmass[1]); ok = false;
    }
    // Cartesian check: O at origin, H1 at (0.957,0,0), H2 from 104.5 degree bend.
    double d2 = std::sqrt(xyz[1][3] * xyz[1][3] + xyz[2][3] * xyz[2][3] + xyz[3][3] * xyz[3][3]);
    if (std::fabs(d2 - 0.957) > 2e-3) { std::fprintf(stderr, "FAIL d2=%f\n", d2); ok = false; }
    double dot = xyz[1][2] * xyz[1][3] + xyz[2][2] * xyz[2][3] + xyz[3][2] * xyz[3][3];
    double costh = dot / (0.957 * d2);
    double angc = std::acos(costh) / 1.7453292519943e-2;
    if (std::fabs(angc - 104.5) > 1.0) { std::fprintf(stderr, "FAIL cart angle=%f\n", angc); ok = false; }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
