// test_batchB7.cpp — tests: gettxt (keyword/comment/title reader, SETUP merge).
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "molkst_C.h"
#include "chanel_C.h"
#include "gettxt.h"

namespace molkst_C {
extern std::string keywrd, koment, title, line, verson;
extern std::vector<std::string> refkey;
extern bool gui, in_house_only, moperr;
extern int numcal;
}
namespace chanel_C {
extern int ir, iw, isetup;
extern std::string input_fn;
}
void upcase(std::string& s, int n) {
    for (int i = 0; i < n && i < (int)s.size(); ++i)
        if (s[i] >= 'a' && s[i] <= 'z') s[i] = (char)(s[i] - 'a' + 'A');
}
void mopend(const std::string&) {}
void web_message(int, const char*) {}
void add_path(std::string&) {}

static void reset() {
    molkst_C::keywrd.clear();
    molkst_C::koment.clear();
    molkst_C::title.clear();
    molkst_C::line.clear();
    molkst_C::numcal = 1;
    molkst_C::gui = false;
    molkst_C::in_house_only = false;
    molkst_C::moperr = false;
    molkst_C::verson = "20.0.0";
    chanel_C::input_fn = "test.mop";
}

int main() {
    bool ok = true;
    // Test 1: plain three-line input.
    {
        reset();
        std::string in = "AM1\nWater molecule\nTest title\n";
        std::istringstream is(in);
        std::cin.rdbuf(is.rdbuf());
        gettxt();
        if (molkst_C::keywrd != "AM1") { std::fprintf(stderr, "FAIL T1 keywrd='%s'\n", molkst_C::keywrd.c_str()); ok = false; }
        if (molkst_C::koment != "Water molecule") { std::fprintf(stderr, "FAIL T1 koment='%s'\n", molkst_C::koment.c_str()); ok = false; }
        if (molkst_C::title != "Test title") { std::fprintf(stderr, "FAIL T1 title='%s'\n", molkst_C::title.c_str()); ok = false; }
        if (molkst_C::refkey[0] != "AM1") { std::fprintf(stderr, "FAIL T1 refkey0='%s'\n", molkst_C::refkey[0].c_str()); ok = false; }
    }
    // Test 2: '+' second keyword line.
    {
        reset();
        std::string in = "PM7 +\nAM1 DFP\nSome comment\nA title\n";
        std::istringstream is(in);
        std::cin.rdbuf(is.rdbuf());
        gettxt();
        if (molkst_C::keywrd.find("PM7") == std::string::npos ||
            molkst_C::keywrd.find("AM1") == std::string::npos ||
            molkst_C::keywrd.find("DFP") == std::string::npos) {
            std::fprintf(stderr, "FAIL T2 keywrd='%s'\n", molkst_C::keywrd.c_str()); ok = false;
        }
    }
    // Test 3: '&' continuation.
    {
        reset();
        std::string in = "AM1&\nPM7\nA title\n";
        std::istringstream is(in);
        std::cin.rdbuf(is.rdbuf());
        gettxt();
        if (molkst_C::keywrd.find("PM7") == std::string::npos) {
            std::fprintf(stderr, "FAIL T3 keywrd='%s'\n", molkst_C::keywrd.c_str()); ok = false;
        }
    }
    // Test 4: SETUP file merge.
    {
        reset();
        std::ofstream fo("setup_test.txt");
        fo << " CHARGE=0 GNORM=0.01\n";
        fo.close();
        std::string in = "SETUP=setup_test.txt\nSetup comment\nSetup title\n";
        std::istringstream is(in);
        std::cin.rdbuf(is.rdbuf());
        gettxt();
        if (molkst_C::keywrd.find("CHARGE") == std::string::npos ||
            molkst_C::keywrd.find("GNORM") == std::string::npos) {
            std::fprintf(stderr, "FAIL T4 keywrd='%s'\n", molkst_C::keywrd.c_str()); ok = false;
        }
        std::remove("setup_test.txt");
    }
    // Test 5: get_text with quoted label.
    {
        std::string t = "SETUP=\"my file.set\"  XYZ";
        std::string got = get_text(t, 7, 0);
        if (got != "my file.set") { std::fprintf(stderr, "FAIL T5 got='%s'\n", got.c_str()); ok = false; }
        if (t[6] != ' ' || t[17] != ' ') { std::fprintf(stderr, "FAIL T5 not blanked: '%s'\n", t.c_str()); ok = false; }
    }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
