// test_batchB5.cpp — tests: pdbout (Brookhaven PDB writer, no HTML branch).
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "pdbout.h"

namespace molkst_C {
extern int numat, natoms, ncomments, nbreaks, maxtxt, numcal, nelecs;
extern std::string verson, line, keywrd, koment, title, formula, geo_ref_name, geo_dat_name;
extern double escf, arc_hof_1, arc_hof_2;
}
namespace chanel_C {
extern int iw;
extern std::string input_fn;
extern bool log;
}
namespace common_arrays_C {
extern std::vector<std::string> txtatm, txtatm1, all_comments;
extern std::vector<std::vector<double>> coord;
extern std::vector<int> nat, breaks;
extern std::vector<double> p;
}
namespace parameters_C {
extern double tore[];
}
namespace MOZYME_C {
extern std::vector<std::string> tyres, tyr;
}
#include "elemts_C.h"
namespace elemts_C {
extern std::vector<std::string> elemnt;
}
extern const char* elemnt(int);

int main() {
    bool ok = true;
    molkst_C::numat = 2;
    molkst_C::natoms = 2;
    molkst_C::ncomments = 0;
    molkst_C::maxtxt = 26;
    molkst_C::keywrd = "";
    molkst_C::verson = "21.2.0";
    chanel_C::input_fn = "test_b5.mop";
    common_arrays_C::txtatm.assign(5, " ");
    common_arrays_C::txtatm1.assign(5, " ");
    common_arrays_C::txtatm[1] = "ATOM      1  O   HOH A   1";
    common_arrays_C::txtatm[2] = "ATOM      2  H   HOH A   1";
    common_arrays_C::coord.assign(4, std::vector<double>(6, 0.0));
    common_arrays_C::coord[1][1] = 0.0; common_arrays_C::coord[2][1] = 0.0; common_arrays_C::coord[3][1] = 0.0;
    common_arrays_C::coord[1][2] = 0.957; common_arrays_C::coord[2][2] = 0.0; common_arrays_C::coord[3][2] = 0.0;
    common_arrays_C::nat.assign(6, 0);
    common_arrays_C::nat[1] = 8;
    common_arrays_C::nat[2] = 1;
    common_arrays_C::breaks.assign(5, 999);
    common_arrays_C::all_comments.assign(5, " ");
    common_arrays_C::p.clear();
    parameters_C::tore[8] = 6.0; parameters_C::tore[1] = 1.0;
    MOZYME_C::tyres = { "", "GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", "ASN",
        "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", "TYR", "MET", "PRO", "PRO", "PRO", "UNK" };
    MOZYME_C::tyr = { "", "G", "A", "V", "L", "I", "S", "T", "D", "N", "K", "E", "Q",
        "R", "H", "F", "C", "W", "Y", "M", "P", "P", "P", "?" };

    std::string out;
    std::FILE* tf = std::tmpfile();
    // pdbout writes to stdout; capture via freopen is awkward, so instead check
    // by writing to a temporary file through dup of stdout.
    fflush(stdout);
    FILE* old = std::fopen("run_b5_capture.tmp", "w");
    if (!old) { std::printf("FAIL cannot open capture\n"); return 1; }
    std::freopen("run_b5_capture.tmp", "w", stdout);
    pdbout(1);
    fflush(stdout);
    std::freopen("CONOUT$", "w", stdout);
    std::fclose(old);
    std::FILE* f = std::fopen("run_b5_capture.tmp", "r");
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = 0;
    std::fclose(f);
    out = buf;
    std::remove("run_b5_capture.tmp");

    if (out.find("HEADER") == std::string::npos) { std::fprintf(stderr, "FAIL no HEADER\n"); ok = false; }
    if (out.find("REMARK") == std::string::npos) { std::fprintf(stderr, "FAIL no REMARK\n"); ok = false; }
    if (out.find(" O") == std::string::npos) { std::fprintf(stderr, "FAIL no O atom line\n"); ok = false; }
    if (out.find(" H") == std::string::npos) { std::fprintf(stderr, "FAIL no H atom line\n"); ok = false; }
    if (out.find("0.957") == std::string::npos) { std::fprintf(stderr, "FAIL no coord\n"); ok = false; }
    if (out.find("END") == std::string::npos) { std::fprintf(stderr, "FAIL no END\n"); ok = false; }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
