// ionout.cpp — C++ translation of MOPAC 2016 "ionout.F90".
#include "ionout.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "elemts_C.h"
#include "common_arrays_C.h"
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

namespace {
using molkst_C::natoms;
using molkst_C::numat;
using molkst_C::maxtxt;
using chanel_C::iw;
using elemts_C::elemnt;
using common_arrays_C::na;
using common_arrays_C::nb;
using common_arrays_C::nc;
using common_arrays_C::txtatm;
using common_arrays_C::labels;
}  // namespace

void ionout(int* ions, int m, int num) {
    std::vector<int> nlabels(natoms + 1, 0);
    int j = 0;
    for (int i = 1; i <= natoms; ++i) {
        if (labels[i] != 99) ++j;
        nlabels[i] = j;
    }
    j = (std::min)(3, numat);
    bool Int = false;
    int klim = 0;
    for (int i = 1; i <= j; ++i) {
        if (na[i] != 0) {
            Int = true;
            klim = 3;
            break;
        }
    }
    std::string blank(40, ' ');
    int i = maxtxt + 5;
    if (i != 5) i = i + 5;
    std::fprintf(stdout, "\n");
    if (Int) {
        std::fprintf(stdout, "     ION   ATOM   WITH   TYPE%s CONNECTIVITY\n", blank.substr(0, (size_t)i).c_str());
        std::fprintf(stdout, "            No.  DUMMYS\n");
    } else {
        std::fprintf(stdout, "     ION   ATOM   WITH   TYPE\n");
        std::fprintf(stdout, "            No.  DUMMYS\n");
    }
    int nabc11 = na[1];
    na[1] = 0;
    if (m != 3) {
        if (maxtxt != 0) {
            for (int i1 = 1; i1 <= num; ++i1) {
                j = ions[(i1 - 1) * 4 + (m - 1)];
                int il, iu;
                if (!txtatm[j].empty() && txtatm[j][0] == '(') {
                    il = 2; iu = maxtxt - 1;
                } else {
                    il = 1; iu = maxtxt;
                }
                std::string t = (il >= 1 && iu <= (int)txtatm[j].size())
                                    ? txtatm[j].substr((size_t)il - 1, (size_t)(iu - il + 1))
                                    : "";
                std::fprintf(stdout, "%7d%7d%7d   %2s(%s)%7d%7d%7d\n", i1, nlabels[j], j,
                             elemnt[labels[j]].c_str(), t.c_str(), na[j], nb[j], nc[j]);
            }
        } else {
            for (int i1 = 1; i1 <= num; ++i1) {
                j = ions[(i1 - 1) * 4 + (m - 1)];
                std::fprintf(stdout, "%7d%7d%7d   %2s%7d%7d%7d\n", i1, nlabels[j], j,
                             elemnt[labels[j]].c_str(), na[j], nb[j], nc[j]);
            }
        }
    } else if (maxtxt != 0) {
        for (int i1 = 1; i1 <= num; ++i1) {
            j = ions[(i1 - 1) * 4 + (m - 1)];
            int il, iu;
            if (!txtatm[j].empty() && txtatm[j][0] == '(') {
                il = 2; iu = maxtxt - 1;
            } else {
                il = 1; iu = maxtxt;
            }
            std::string t = (il >= 1 && iu <= (int)txtatm[j].size())
                                ? txtatm[j].substr((size_t)il - 1, (size_t)(iu - il + 1))
                                : "";
            if (klim == 0) {
                std::fprintf(stdout, "%7d%7d%7d   %2s(%s)%s%+d\n", i1, nlabels[j], j,
                             elemnt[labels[j]].c_str(), t.c_str(), "  CHARGE: ", ions[(i1 - 1) * 4 + 3]);
            } else {
                std::fprintf(stdout, "%7d%7d%7d   %2s(%s)%s%+d%7d%7d%7d\n", i1, nlabels[j], j,
                             elemnt[labels[j]].c_str(), t.c_str(), "  CHARGE: ", ions[(i1 - 1) * 4 + 3],
                             na[j], nb[j], nc[j]);
            }
        }
    } else {
        for (int i1 = 1; i1 <= num; ++i1) {
            j = ions[(i1 - 1) * 4 + (m - 1)];
            if (klim == 0) {
                std::fprintf(stdout, "%7d%7d%7d   %2s%s%+d\n", i1, nlabels[j], j,
                             elemnt[labels[j]].c_str(), "  CHARGE: ", ions[(i1 - 1) * 4 + 3]);
            } else {
                std::fprintf(stdout, "%7d%7d%7d   %2s%7d%7d%7d%s%+d\n", i1, nlabels[j], j,
                             elemnt[labels[j]].c_str(), na[j], nb[j], nc[j], "  CHARGE: ",
                             ions[(i1 - 1) * 4 + 3]);
            }
        }
    }
    na[1] = nabc11;
}
