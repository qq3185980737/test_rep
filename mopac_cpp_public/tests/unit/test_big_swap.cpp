// test_big_swap.cpp
#include <cstdio>
#include <string>

#include "big_swap.h"
#include "molkst_C.h"

int main() {
    // Add a new keyword into a keywrd that has no long space-run.
    molkst_C::keywrd = " AAA BBB";
    l_control("CCC", 1);
    bool r1 = molkst_C::keywrd.find("CCC") != std::string::npos;
    std::printf("after add CCC: '%s'  %s\n", molkst_C::keywrd.c_str(), r1 ? "PASS" : "FAIL");

    // Remove it.
    l_control("CCC", -1);
    bool r2 = molkst_C::keywrd.find("CCC") == std::string::npos;
    std::printf("after del CCC: '%s'  %s\n", molkst_C::keywrd.c_str(), r2 ? "PASS" : "FAIL");

    // Replace an existing keyword: old AAA should be gone, new AAA stays once.
    molkst_C::keywrd = " AAA BBB";
    l_control("AAA", 1);
    int cnt = 0;
    size_t pos = 0;
    while ((pos = molkst_C::keywrd.find("AAA", pos)) != std::string::npos) { ++cnt; ++pos; }
    bool r3 = (cnt == 1);
    std::printf("after replace AAA: '%s' (AAA count=%d) %s\n",
                molkst_C::keywrd.c_str(), cnt, r3 ? "PASS" : "FAIL");

    return (r1 && r2 && r3) ? 0 : 1;
}
