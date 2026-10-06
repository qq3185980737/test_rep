// test_batchB2.cpp — tests: geout (geometry print) and getdat (data-set read).
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "symmetry_C.h"
#include "maps_C.h"
#include "parameters_C.h"
#include "getdat.h"

namespace molkst_C {
extern int numat, natoms, nvar, ndep, maxtxt, ncomments, ijulian, run;
extern bool gui, in_house_only;
extern std::string keywrd, line, title, koment, jobnam, verson;
extern double arc_hof_1, arc_hof_2;
}
namespace chanel_C {
extern int iw, iw0;
extern std::string job_fn;
}
namespace common_arrays_C {
extern std::vector<int> labels, na, nb, nc, nat;
extern std::vector<std::vector<double>> geo;
extern std::vector<std::vector<int>> loc;
extern std::vector<std::string> txtatm, all_comments;
extern std::vector<char> l_atom;
extern std::vector<double> atmass;
}
namespace parameters_C {
extern std::vector<double> ams;   // tore is double[N_PARAM+1] in parameters_C.h
}
namespace symmetry_C {
extern std::vector<int> locpar, idepfn, locdep;
extern std::vector<double> depmul;
}
namespace maps_C {
extern int lpara1, latom1, lpara2, latom2, lparam, latom;
}
namespace elemts_C {
extern std::vector<std::string> elemnt;
}
extern const char* elemnt(int);
extern std::vector<std::string> getdat_lines;

// Stubs for routines not under test.
void wrttxt(int) {}
void mopend(const std::string&) {}
void to_screen(const std::string&) {}
void init_filenames() {}
void web_message(const char*) {}
void upcase(std::string& s, int n) {
    for (int i = 0; i < n && i < (int)s.size(); ++i)
        if (s[i] >= 'a' && s[i] <= 'z') s[i] = (char)(s[i] - 'a' + 'A');
}
extern "C" double reada_(const char* s, int* i, int len);
double reada(const std::string& s, int i) {
    int p = i;
    return reada_(s.c_str(), &p, (int)s.size());
}
extern "C" double reada_(const char* s, int* i, int len) {
    (void)len;
    std::string t = (s == nullptr) ? "" : s;
    if (*i < 0 || *i > (int)t.size()) *i = 0;
    t = t.substr(*i);
    double v = 0.0;
    size_t pos = 0;
    try { v = std::stod(t, &pos); } catch (...) { v = 0.0; }
    *i = *i + (int)pos;
    return v;
}
void chrge(const std::vector<double>&, std::vector<double>&) {}

#include "geout.h"

int main() {
    bool ok = true;
    // ---------------- geout (mode=1, MOPAC output style) ----------------
    molkst_C::numat = 3; molkst_C::natoms = 3; molkst_C::nvar = 6; molkst_C::ndep = 0;
    molkst_C::maxtxt = 0; molkst_C::gui = false; molkst_C::in_house_only = false;
    molkst_C::keywrd = "";
    chanel_C::iw = 6;
    common_arrays_C::geo.assign(4, std::vector<double>(4, 0.0));
    common_arrays_C::geo[1][1] = 0.0; common_arrays_C::geo[2][1] = 0.0; common_arrays_C::geo[3][1] = 0.0;
    common_arrays_C::geo[1][2] = 0.96; common_arrays_C::geo[2][2] = 0.0; common_arrays_C::geo[3][2] = 0.0;
    common_arrays_C::geo[1][3] = 0.96; common_arrays_C::geo[2][3] = 104.5; common_arrays_C::geo[3][3] = 0.0;
    common_arrays_C::na.assign(4, 0); common_arrays_C::nb.assign(4, 0); common_arrays_C::nc.assign(4, 0);
    common_arrays_C::na[2] = 1; common_arrays_C::nb[2] = 0; common_arrays_C::nc[2] = 0;
    common_arrays_C::na[3] = 1; common_arrays_C::nb[3] = 2; common_arrays_C::nc[3] = 0;
    common_arrays_C::labels.assign(4, 1); common_arrays_C::labels[1] = 8; common_arrays_C::labels[2] = 1; common_arrays_C::labels[3] = 1;
    common_arrays_C::nat.assign(4, 1); common_arrays_C::nat[1] = 8; common_arrays_C::nat[2] = 1; common_arrays_C::nat[3] = 1;
    common_arrays_C::loc.assign(3, std::vector<int>(8, 0));
    for (int i = 1; i <= 6; ++i) { common_arrays_C::loc[1][i] = (i + 2) / 3; common_arrays_C::loc[2][i] = (i - 1) % 3 + 1; }
    common_arrays_C::txtatm.assign(4, " ");
    common_arrays_C::l_atom.assign(4, 1);
    common_arrays_C::atmass.assign(4, 0.0);
    common_arrays_C::atmass[1] = 15.999; common_arrays_C::atmass[2] = 1.0079; common_arrays_C::atmass[3] = 1.0079;
    parameters_C::ams.resize(108, 0.0);
    parameters_C::ams[1] = 1.0079; parameters_C::ams[8] = 15.999;
    geout(1);
    // (Output printed; smoke test only.)

    // ---------------- getdat (run=2 path, read test file) ----------------
    {
        std::ofstream f("test_B2.mop");
        f << " AM1 CHARGE=0\n"
          << "water test\n"
          << "\n"
          << " O  0.000000  0  0.000000  0  0.000000  0\n"
          << " H  0.960000  1  0.000000  2  0.000000  0\n"
          << " H  0.960000  1  104.500000  2  0.000000  3\n"
          << "\n";
        f.close();
        molkst_C::run = 2;
        molkst_C::jobnam = "test_B2";
        molkst_C::ijulian = 200;
        molkst_C::verson = "2016";
        chanel_C::iw0 = -1;
        getdat(0, 0);
        if (getdat_lines.empty()) { std::printf("FAIL getdat no lines\n"); ok = false; }
        if (getdat_lines[0].find("AM1") == std::string::npos) {
            std::printf("FAIL getdat kwd='%s'\n", getdat_lines[0].c_str());
            ok = false;
        }
        if (molkst_C::natoms < 5) { std::printf("FAIL getdat natoms=%d\n", molkst_C::natoms); ok = false; }
    }

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
