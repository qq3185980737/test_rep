// wrttxt.cpp
#include "wrttxt.h"
#include <string>
#include <cstdio>
#include <cctype>
#include "molkst_C.h"
using namespace molkst_C;
void wrttxt(int iprt) {
    bool l_chains = false, l_start = false;
    for (int i = 0; i < 6; ++i) {
        if (refkey[i].find(" NULL") != std::string::npos) break;
        line = " " + refkey[i];
        for (auto& ch : line) ch = (char)toupper(ch);
        if (!l_chains) l_chains = (line.find(" CHAINS") != std::string::npos);
        if (!l_start)  l_start  = (line.find(" START_RES") != std::string::npos);
    }
    size_t i = keywrd.find(" CHAINS");
    if (i != std::string::npos && !l_chains) {
        size_t j = keywrd.find(")", i + 7) + i + 7;
        refkey[0] = keywrd.substr(i, j - i) + refkey[0];
    }
    i = keywrd.find(" START_RES");
    if (i != std::string::npos && !l_start) {
        size_t j = keywrd.find(")", i + 10) + i + 10;
        refkey[0] = keywrd.substr(i, j - i) + refkey[0];
    }
    for (int k = 0; k < 6; ++k) {
        if (refkey[k].find(" NULL") != std::string::npos) break;
        std::fprintf(stdout, "%s\n", refkey[k].c_str());
    }
    if (koment.find(" NULL ") == std::string::npos) std::fprintf(stdout, "%s\n", koment.c_str());
    if (title.find(" NULL ") == std::string::npos) std::fprintf(stdout, "%s\n", title.c_str());
}
