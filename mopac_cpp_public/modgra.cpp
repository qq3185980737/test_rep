// modgra.cpp — C++ translation of MOPAC 2016 "modgra.F90".
#include "modgra.h"
extern void build_res_start_etc();
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {
using molkst_C::numat;
using molkst_C::nvar;
using molkst_C::id;
using chanel_C::iw;
using common_arrays_C::grad;
using common_arrays_C::loc;
using common_arrays_C::txtatm;
using MOZYME_C::nres;
using MOZYME_C::at_res;
using MOZYME_C::res_start;
}  // namespace

// build_res_start_etc: work out which residue each atom belongs to and the
// location of the first atom in each residue (from modchg.F90).

void modgra() {
    std::vector<double> rgrad(numat + id + 1, 0.0), res_grad_b(numat + id + 1, 0.0),
        res_grad_s(numat + id + 1, 0.0);
    build_res_start_etc();
    for (int i = 1; i <= nvar; ++i) {
        int j = loc[1][i];
        rgrad[j] = rgrad[j] + grad[i] * grad[i];
    }
    bool l_protein = false;
    for (int i = 1; i <= numat; ++i) {
        int j = at_res[i];
        std::string t = txtatm[i];
        std::string tag = (t.size() >= 15) ? t.substr(12, 3) : "";
        if (t.rfind("ATOM  ", 0) == 0 && (tag == " CA" || tag == " N " || tag == " C ")) {
            res_grad_b[j] = res_grad_b[j] + rgrad[i];
            l_protein = true;
        } else {
            res_grad_s[j] = res_grad_s[j] + rgrad[i];
        }
    }
    double sum = 0.0;
    for (int i = 1; i <= nres; ++i) {
        res_grad_b[i] = std::sqrt(res_grad_b[i]);
        res_grad_s[i] = std::sqrt(res_grad_s[i]);
        sum = sum + res_grad_b[i] + res_grad_s[i];
    }
    std::fprintf(stdout, "\n");
    if (l_protein) {
        if (sum < 0.5) {
            std::fprintf(stdout, "%26s\n", " ALL GRADIENTS FOR ALL BACKBONES PLUS SIDE CHAINS ARE VERY SMALL");
            std::fprintf(stdout, "\n");
        } else {
            std::fprintf(stdout, "%26s\n", " GRADIENTS FOR ALL BACKBONES");
            std::fprintf(stdout, "     Residue        Backbone    Side-Chain        Total\n");
            std::fprintf(stdout, "\n");
            for (int i = 1; i <= nres; ++i) {
                std::string name = (txtatm[res_start[i]].size() >= 20) ? txtatm[res_start[i]].substr(17, 3) : "";
                if (name != "HOH" || std::fabs(res_grad_s[i]) > 1.0)
                    std::fprintf(stdout, "%4s%15.3f%12.3f%15.3f\n", name.c_str(),
                                 res_grad_b[i], res_grad_s[i], res_grad_b[i] + res_grad_s[i]);
            }
        }
    } else {
        if (sum < 0.5) {
            std::fprintf(stdout, "%26s\n", " ALL GRADIENTS FOR ALL GROUPS ARE VERY SMALL");
            std::fprintf(stdout, "\n");
        } else {
            std::fprintf(stdout, "%27s\n", " GRADIENTS FOR ALL GROUPS");
            std::fprintf(stdout, "%11s%10s\n", "     Group       Gradient", "");
            std::fprintf(stdout, "\n");
            for (int i = 1; i <= nres; ++i) {
                std::string r1 = (txtatm[res_start[i]].size() >= 20) ? txtatm[res_start[i]].substr(17, 3) : "";
                std::string r2 = (txtatm[res_start[i]].size() >= 26) ? txtatm[res_start[i]].substr(22, 4) : "";
                std::fprintf(stdout, "%14s%13.3f\n", (r1 + " " + r2).c_str(), res_grad_s[i]);
            }
        }
    }
}
