// test_batchB4.cpp — tests: getpdb (Brookhaven PDB reader).
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "getpdb.h"

namespace molkst_C {
extern int maxtxt, natoms, numat, na1, nbreaks, ncomments;
extern std::string keywrd, line;
}
namespace common_arrays_C {
extern std::vector<std::string> txtatm, txtatm1, all_comments;
extern std::vector<double> atmass;
extern std::vector<int> na, nb, nc, labels, breaks;
extern std::vector<std::vector<int>> lopt;
extern std::vector<std::vector<double>> coord, break_coords;
}
namespace parameters_C {
extern std::vector<double> ams;
}
void upcase(std::string& s, int n) {
    for (int i = 0; i < n && i < (int)s.size(); ++i)
        if (s[i] >= 'a' && s[i] <= 'z') s[i] = (char)(s[i] - 'a' + 'A');
}
extern "C" double reada_(const char*, int*, int);
double reada(const std::string& s, int i) {
    int p = i;
    return reada_(s.c_str(), &p, (int)s.size());
}
extern "C" double reada_(const char* s, int* i, int len) {
    (void)len;
    std::string t = (s == nullptr) ? "" : s;
    if (*i < 1 || *i > (int)t.size()) *i = 1;
    t = t.substr(*i - 1);
    double v = 0.0;
    size_t pos = 0;
    try { v = std::stod(t, &pos); } catch (...) { v = 0.0; }
    *i = *i + (int)pos;
    return v;
}
void mopend(const std::string&) {}
void web_message(int, const char*) {}
void xyzint(const std::vector<double>&, int, const std::vector<int>&,
            const std::vector<int>&, const std::vector<int>&, double, std::vector<double>&) {}

#include "elemts_C.h"
namespace elemts_C {
extern std::vector<std::string> elemnt;
}
extern const char* elemnt(int);

int main() {
    bool ok = true;
    // Build PDB ATOM lines field-by-field (columns 1-based).
    auto atom_line = [](int seq, const char* aname, const char* res, char chain,
                        int rnum, double x, double y, double z, const char* elem2) {
        std::string L(80, ' ');
        L.replace(0, 4, "ATOM");
        char buf[8];
        std::snprintf(buf, sizeof(buf), "%5d", seq); L.replace(6, 5, buf);
        L.replace(12, 4, aname);
        L.replace(17, 3, res);
        L[21] = chain;
        std::snprintf(buf, sizeof(buf), "%4d", rnum); L.replace(22, 4, buf);
        std::snprintf(buf, sizeof(buf), "%8.3f", x); L.replace(30, 8, buf);
        std::snprintf(buf, sizeof(buf), "%8.3f", y); L.replace(38, 8, buf);
        std::snprintf(buf, sizeof(buf), "%8.3f", z); L.replace(46, 8, buf);
        L.replace(54, 6, "  1.00");
        L.replace(60, 6, " 20.00");
        L.replace(76, 2, elem2); // element, right-aligned in 77-78
        return L;
    };
    std::string pdb =
        "HEADER    TEST PROTEIN\n" +
        atom_line(1, " N  ", "ALA", 'A', 1, 1.000, 2.000, 3.000, " N") + "\n" +
        atom_line(2, " CA ", "ALA", 'A', 1, 1.100, 2.100, 3.100, " C") + "\n" +
        atom_line(3, " HA ", "ALA", 'A', 1, 1.200, 2.200, 3.200, " H") + "\n" +
        "END\n";
    std::istringstream in(pdb);
    std::cin.rdbuf(in.rdbuf());
    molkst_C::keywrd = "";
    common_arrays_C::coord.assign(4, std::vector<double>(12, 0.0));
    common_arrays_C::break_coords.assign(4, std::vector<double>(402, 0.0));
    common_arrays_C::lopt.assign(4, std::vector<int>(12, 0));
    common_arrays_C::breaks.assign(402, 0);
    common_arrays_C::txtatm.assign(12, " ");
    common_arrays_C::txtatm1.assign(12, " ");
    common_arrays_C::atmass.assign(12, 0.0);
    common_arrays_C::na.assign(12, 0);
    common_arrays_C::nb.assign(12, 0);
    common_arrays_C::nc.assign(12, 0);
    common_arrays_C::labels.assign(12, 0);
    parameters_C::ams.resize(108, 0.0);
    parameters_C::ams[1] = 1.0079; parameters_C::ams[6] = 12.011; parameters_C::ams[7] = 14.007;
    std::vector<std::vector<double>> geo(4, std::vector<double>(12, 0.0));
    getpdb(geo);
    if (molkst_C::natoms != 3) { std::printf("FAIL natoms=%d\n", molkst_C::natoms); ok = false; }
    if (molkst_C::numat != 3) { std::printf("FAIL numat=%d\n", molkst_C::numat); ok = false; }
    if (common_arrays_C::labels[1] != 7) { std::printf("FAIL label1=%d\n", common_arrays_C::labels[1]); ok = false; }
    if (common_arrays_C::labels[2] != 6) { std::printf("FAIL label2=%d\n", common_arrays_C::labels[2]); ok = false; }
    if (common_arrays_C::labels[3] != 1) { std::printf("FAIL label3=%d\n", common_arrays_C::labels[3]); ok = false; }
    if (std::fabs(geo[1][2] - 1.1) > 1.e-4) { std::printf("FAIL x=%f\n", geo[1][2]); ok = false; }
    if (std::fabs(geo[2][2] - 2.1) > 1.e-4) { std::printf("FAIL y=%f\n", geo[2][2]); ok = false; }
    if (std::fabs(geo[3][2] - 3.1) > 1.e-4) { std::printf("FAIL z=%f\n", geo[3][2]); ok = false; }
    if (common_arrays_C::lopt[1][1] != 1) { std::printf("FAIL lopt\n"); ok = false; }
    if (common_arrays_C::na[1] != 0) { std::printf("FAIL na\n"); ok = false; }
    if (std::fabs(common_arrays_C::atmass[2] - 12.011) > 1.e-3) { std::printf("FAIL atmass=%f\n", common_arrays_C::atmass[2]); ok = false; }
    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
