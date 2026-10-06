// findn1.cpp — C++ translation of MOPAC 2016 "findn1.F90".

#include "findn1.h"

#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;

bool ring_6(int atom_1, int atom_2, int atom_6) {
    int nii = nbonds[atom_2];
    int njj = nbonds[atom_6];
    for (int i = 1; i <= nii; ++i) {
        int atom_3 = ibonds[i][atom_2];
        if (atom_3 == atom_1) continue;
        for (int j = 1; j <= njj; ++j) {
            int atom_5 = ibonds[j][atom_6];
            if (atom_5 == atom_1) continue;
            for (int k = 1; k <= nbonds[atom_5]; ++k) {
                int kka = ibonds[k][atom_5];
                if (kka == atom_6) continue;
                for (int l = 1; l <= nbonds[atom_3]; ++l) {
                    int atom_4 = ibonds[l][atom_3];
                    if (atom_4 == atom_2) continue;
                    if (kka == atom_4) {
                        if (atom_1 == atom_3 || atom_1 == atom_4 || atom_1 == atom_5) return false;
                        if (atom_2 == atom_3 || atom_2 == atom_4 || atom_2 == atom_5) return false;
                        if (atom_3 == atom_4 || atom_3 == atom_5 || atom_3 == atom_6) return false;
                        if (atom_4 == atom_5 || atom_4 == atom_6 || atom_5 == atom_6) return false;
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void findn1(int& n1, const bool* ioptl, int& io) {
    using namespace common_arrays_C;
    static bool ok = true;
    n1 = 0;
    int iatom = 0;
    bool peptide = false;
    bool is_ring = false;
    int j = 0;
    int ii, jj, kk, ll, mm;
    int i, k, l;
    bool outer_done = false;
    for (int iatom_loop = 1; iatom_loop <= molkst_C::numat && !outer_done; ++iatom_loop) {
        iatom = iatom_loop;
        if (!ioptl[iatom]) {
            if (nat[iatom] == 8) {
                for (i = 1; i <= nbonds[iatom]; ++i) {
                    if (nat[ibonds[i][iatom]] == 6) {
                        l = ibonds[i][iatom];
                        for (k = 1; k <= nbonds[l]; ++k) {
                            if (nat[ibonds[k][l]] == 6) {
                                j = ibonds[k][l];
                                for (jj = 1; jj <= nbonds[j]; ++jj) {
                                    kk = ibonds[jj][j];
                                    if (kk != l) {
                                        if (ring_6(j, l, kk)) continue;
                                    }
                                    if (nat[kk] != 6 && nat[kk] != 1 && nat[kk] != 7) continue;
                                }
                                for (ii = 1; ii <= nbonds[j]; ++ii) {
                                    if (nat[ibonds[ii][j]] == 7 && !ioptl[ibonds[ii][j]]) {
                                        n1 = ibonds[ii][j];
                                        is_ring = false;
                                        for (jj = 1; jj <= nbonds[n1]; ++jj) {
                                            kk = ibonds[jj][n1];
                                            for (ll = 1; ll <= nbonds[n1]; ++ll) {
                                                mm = ibonds[ll][n1];
                                                if (nat[kk] > 1 && nat[mm] > 1 && kk > mm) {
                                                    if (ring_6(n1, mm, kk)) is_ring = true;
                                                }
                                            }
                                        }
                                        if (is_ring) { n1 = 0; }
                                        else { outer_done = true; break; }
                                    }
                                }
                                if (outer_done) break;
                            }
                        }
                        if (outer_done) break;
                    }
                }
            }
        }
    }
    io = iatom;
    if (n1 == 0) return;
    j = 0;
    for (i = 1; i <= nbonds[n1]; ++i) {
        if (nat[ibonds[i][n1]] != 1) ++j;
    }
    if (j == 1) return;
    peptide = false;
    int k2 = 0;
    for (int ires_loc = 1; ires_loc <= 1000000; ++ires_loc) {
        int found_c = 0;
        for (i = 1; i <= nbonds[n1]; ++i) {
            j = ibonds[i][n1];
            if (nat[j] == 6) {
                bool has_o = false;
                for (k = 1; k <= nbonds[j]; ++k) {
                    if (nat[ibonds[k][j]] == 8) { has_o = true; break; }
                }
                if (!has_o) continue;
                found_c = j;
                bool has_c = false;
                for (k = 1; k <= nbonds[j]; ++k) {
                    if (nat[ibonds[k][j]] == 6) { peptide = true; has_c = true; k2 = k; break; }
                }
                if (has_c) break;
            }
        }
        j = found_c;
        j = ibonds[k2][j];
        int k_found = 0;
        for (k = 1; k <= nbonds[j]; ++k) {
            if (nat[ibonds[k][j]] == 7) { k_found = k; break; }
        }
        if (peptide) {
            for (k = 1; k <= nbonds[j]; ++k) {
                if (nat[ibonds[k][j]] == 1) n1 = ibonds[k][j];
            }
            ok = ok && (n1 == ibonds[k][j]);
            molkst_C::isok = ok;
        }
        if (k_found == 0) return;
        if (n1 == ibonds[k_found][j]) return;
        if (!ioptl[ibonds[k_found][j]]) n1 = ibonds[k_found][j];
        j = 0;
        for (ii = 1; ii <= nbonds[n1]; ++ii) {
            if (nat[ibonds[ii][n1]] != 1) ++j;
        }
        if (j == 1) return;
    }
}
