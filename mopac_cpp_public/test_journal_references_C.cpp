// test_journal_references_C.cpp — spot-check parsed reference strings.
#include "journal_references_C.h"

#include <cstdio>
#include <string>

using namespace journal_references_C;

int main() {
    bool ok = true;
    auto has = [&](const std::string& s, const char* frag) { return s.find(frag) != std::string::npos; };
    // refrm1: H entry contains RM1 citation
    if (!has(refrm1[1], "RM1") || !has(refrm1[1], "ROCHA")) { std::printf("FAIL refrm1[1]: %s\n", refrm1[1].c_str()); ok = false; }
    // refmd: Na entry contains MNDO/d
    if (!has(refmd[11], "MNDO/d") || !has(refmd[11], "THIEL")) { std::printf("FAIL refmd[11]: %s\n", refmd[11].c_str()); ok = false; }
    // refpm6 lanthanide La
    if (!has(refpm6[57], "PM6") || !has(refpm6[57], "Freire")) { std::printf("FAIL refpm6[57]: %s\n", refpm6[57].c_str()); ok = false; }
    // refpm7 lanthanide Lu
    if (!has(refpm7[71], "PM7") || !has(refpm7[71], "Dutra")) { std::printf("FAIL refpm7[71]: %s\n", refpm7[71].c_str()); ok = false; }
    // refmn/refam/refpm3 have no data statements in Fortran -> stay empty
    if (!refmn[1].empty() || !refam[6].empty() || !refpm3[6].empty()) { std::printf("FAIL empty arrays\n"); ok = false; }
    // untouched slots empty
    if (!refrm1[2].empty()) { std::printf("FAIL refrm1[2] should be empty\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
