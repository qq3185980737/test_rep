// wrt_diffs.cpp — C++ translation of "wrt_diffs.F90".
// Given two geometries (geo and geoa), recompute connectivity for each, pair
// up bonds by matching element, and print bond-length differences largest first.

#include "wrt_diffs.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "set_up_dentate.h"
#include "upcase.h"
#include "big_swap.h"
#include "post_scf_corrections.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

// F90 analyze_h_bonds @109: compare hydrogen bonds between the current
// geometry (geo) and the reference geometry (geoa): common bonds, bonds only
// in the dataset, and bonds only in GEO_REF, with non-covalent energies.
namespace {
std::string hb_trim(const std::string& s) {
    size_t e = s.size();
    while (e > 0 && s[e - 1] == ' ') --e;
    return s.substr(0, e);
}
void get_H_bonds() {   // F90 get_H_bonds @257: evaluate H-bond corrections
    bool l_grad = false;
    double correction = 0.0;
    post_scf_corrections(correction, l_grad);
}
}  // namespace

void analyze_h_bonds() {
    std::vector<std::string> set_1_H_txt(numat + 1), set_2_H_txt(numat + 1);
    std::vector<std::string> set_1_H_dist(numat + 1), set_2_H_dist(numat + 1);
    std::vector<double> set_1_H_energy(numat + 1, 0.0), set_2_H_energy(numat + 1, 0.0);
    std::vector<bool> l_used_1(numat + 1, false), l_used_2(numat + 1, false);
    l_control("PRT DISP(1.0) SILENT", (int)std::string("PRT DISP(1.0) SILENT").size(), 1);
    get_H_bonds();
    double E_1_disp = E_disp, E_1_hb = E_hb, E_1_hh = E_hh;
    int set_1_P_Hbonds = P_Hbonds;
    for (int i = 1; i <= P_Hbonds; ++i) {
        size_t sz = H_txt[i].size();
        set_1_H_txt[i] = H_txt[i].substr(0, std::min<size_t>(31, sz)) +
                         hb_trim(H_txt[i].substr(std::min<size_t>(36, sz),
                                                 std::min<size_t>(32, sz > 36 ? sz - 36 : 0)));
        set_1_H_dist[i] = H_txt[i].substr(std::min<size_t>(31, sz),
                                          std::min<size_t>(5, sz > 31 ? sz - 31 : 0));
    }
    for (int i = 1; i <= P_Hbonds; ++i) set_1_H_energy[i] = H_energy[i];
    for (int i = 1; i <= numat; ++i) {
        coord[0][i] = geoa[1][i];
        coord[1][i] = geoa[2][i];
        coord[2][i] = geoa[3][i];
    }
    ++numcal;
    get_H_bonds();
    double E_2_disp = E_disp, E_2_hb = E_hb, E_2_hh = E_hh;
    int set_2_P_Hbonds = P_Hbonds;
    for (int i = 1; i <= P_Hbonds; ++i) {
        size_t sz = H_txt[i].size();
        set_2_H_txt[i] = H_txt[i].substr(0, std::min<size_t>(31, sz)) +
                         hb_trim(H_txt[i].substr(std::min<size_t>(36, sz),
                                                 std::min<size_t>(32, sz > 36 ? sz - 36 : 0)));
        set_2_H_dist[i] = H_txt[i].substr(std::min<size_t>(31, sz),
                                          std::min<size_t>(5, sz > 31 ? sz - 31 : 0));
    }
    for (int i = 1; i <= P_Hbonds; ++i) set_2_H_energy[i] = H_energy[i];

    std::printf("\n\n                    Analysis of Non-Covalent Interactions\n");

    int n_common = 0;
    for (int i = 1; i <= set_1_P_Hbonds; ++i)
        for (int j = 1; j <= set_2_P_Hbonds; ++j)
            if (set_1_H_txt[i] == set_2_H_txt[j]) {
                l_used_1[i] = true;
                l_used_2[j] = true;
                ++n_common;
            }
    if (job_fn == geo_dat_name) line = "dataset";
    else line = "GEO_DAT";
    std::printf("\n          Total non-covalent energy of %s system: %9.2f Kcal/mol\n",
                line.c_str(), E_1_disp + E_1_hb + E_1_hh);
    std::printf("          Total non-covalent energy of GEO_REF system: %9.2f Kcal/mol\n",
                E_2_disp + E_2_hb + E_2_hh);
    double sum = E_1_disp + E_1_hb + E_1_hh - E_2_disp - E_2_hb - E_2_hh;
    std::printf("                                           Difference: %9.2f Kcal/mol\n", sum);
    if (n_common > 0)
        std::printf("\n          Number of hydrogen bonds common to both systems:%5d\n", n_common);

    // Bonds in the dataset but not in GEO_REF, lowest energy first.
    int jj = 0;
    bool l_prt = true;
    for (int i = 1; i <= set_1_P_Hbonds; ++i) {
        if (!l_used_1[i]) {
            ++jj;
            set_1_H_txt[jj] = set_1_H_txt[i];
            set_1_H_energy[jj] = set_1_H_energy[i];
            set_1_H_dist[jj] = set_1_H_dist[i];
        }
    }
    double sum_a = 0.0;
    for (int l = 1; l <= jj; ++l) {
        sum = 10.0;
        int k = 0;
        for (int i = 1; i <= jj; ++i)
            if (set_1_H_energy[i] < sum) { k = i; sum = set_1_H_energy[i]; }
        if (sum > -0.8) break;
        if (set_1_H_dist[k] > "3.000") continue;
        if (l_prt) {
            if (job_fn == geo_dat_name) line = "dataset";
            else line = "GEO_DAT";
            std::printf("\n\n                            Hydrogen bonds in %s but not in GEO_REF\n",
                        line.c_str());
            std::printf("\n               Donor atom          Hydrogen atom      H-bond length(A)     Energy    Sum\n\n");
            l_prt = false;
        }
        sum_a += set_1_H_energy[k];
        std::printf("%4d   %s        %s%17.3f%9.3f\n", l, set_1_H_txt[k].c_str(),
                    set_1_H_dist[k].c_str(), set_1_H_energy[k], sum_a);
        set_1_H_energy[k] = 5.0;
    }

    // Bonds in GEO_REF but not in the dataset, lowest energy first.
    jj = 0;
    l_prt = true;
    for (int i = 1; i <= set_2_P_Hbonds; ++i) {
        if (!l_used_2[i]) {
            ++jj;
            set_2_H_txt[jj] = set_2_H_txt[i];
            set_2_H_energy[jj] = set_2_H_energy[i];
            set_2_H_dist[jj] = set_2_H_dist[i];
        }
    }
    sum_a = 0.0;
    for (int l = 1; l <= jj; ++l) {
        sum = 10.0;
        int k = 0;
        for (int i = 1; i <= jj; ++i)
            if (set_2_H_energy[i] < sum) { k = i; sum = set_2_H_energy[i]; }
        if (sum > -0.8) break;
        if (set_2_H_dist[k] > "3.000") continue;
        if (l_prt) {
            if (job_fn == geo_dat_name) line = "dataset";
            else line = "GEO_DAT";
            std::printf("\n\n                            Hydrogen bonds in GEO_REF but not in %s\n",
                        line.c_str());
            std::printf("\n               Donor atom          Hydrogen atom      H-bond length(A)     Energy    Sum\n\n");
            l_prt = false;
        }
        sum_a += set_2_H_energy[k];
        std::printf("%4d   %s        %s%17.3f%9.3f\n", l, set_2_H_txt[k].c_str(),
                    set_2_H_dist[k].c_str(), set_2_H_energy[k], sum_a);
        set_2_H_energy[k] = 5.0;
    }
}

// wrt_diffs(): F90 wrt_diffs.F90 main body not yet translated (needs geo/geoa
// connectivity comparison + bond-length-difference printout).  Stub keeps the
// geo_ref.cpp driver linkable; call is a no-op for now.
void wrt_diffs() {
    // TODO: translate wrt_diffs.F90 (bond-by-bond geo vs geoa printout).
}
