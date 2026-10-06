// atomrs.cpp — C++ translation of MOPAC 2016 "atomrs.F90".
// Protein residue identification. Output unit iw maps to stdout; unported
// external routines (dihed, greek, mopend) are stubbed.

#include "atomrs.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

using namespace common_arrays_C;
using namespace elemts_C;

namespace {
// --- stubs for not-yet-portable external routines ---
void dihed(std::vector<std::vector<double>>&, int, int, int, int, double& sum) {
    sum = 0.0;  // torsion angle stub; not exercised in peptide_n tests
}
void greek(int) {}
[[noreturn]] void mopend(const char* msg) {
    std::printf("MOPEND: %s\n", msg);
    std::exit(1);
}

// mres(5, size_mres): element counts (C,N,O,S,H) per residue template.
std::vector<std::vector<int>> build_mres() {
    std::vector<std::vector<int>> m(24, std::vector<int>(6, 0));
    m[1]  = {0, 2, 1, 1, 0, 3};
    m[2]  = {0, 3, 1, 1, 0, 5};
    m[3]  = {0, 5, 1, 1, 0, 9};
    m[4]  = {0, 6, 1, 1, 0, 11};
    m[5]  = {0, 6, 1, 1, 0, 11};
    m[6]  = {0, 3, 1, 2, 0, 5};
    m[7]  = {0, 4, 1, 2, 0, 7};
    m[8]  = {0, 4, 1, 3, 0, 5};
    m[9]  = {0, 4, 2, 2, 0, 6};
    m[10] = {0, 6, 2, 1, 0, 13};
    m[11] = {0, 5, 1, 3, 0, 7};
    m[12] = {0, 5, 2, 2, 0, 8};
    m[13] = {0, 6, 4, 1, 0, 13};
    m[14] = {0, 6, 3, 1, 0, 7};
    m[15] = {0, 9, 1, 1, 0, 9};
    m[16] = {0, 3, 1, 1, 1, 5};
    m[17] = {0, 11, 2, 1, 0, 10};
    m[18] = {0, 9, 1, 2, 0, 9};
    m[19] = {0, 5, 1, 1, 1, 9};
    m[20] = {0, 6, 2, 1, 0, 7};
    m[21] = {0, 5, 1, 2, 0, 9};
    m[22] = {0, 5, 1, 2, 0, 9};
    m[23] = {0, 0, 0, 0, 0, 0};
    return m;
}
std::vector<std::vector<int>> mres = build_mres();
const int natomr = 10000;
std::vector<int> once(6, 0);
int n_once = 0;
}

bool peptide_n(int l) {
    using namespace common_arrays_C;
    if (nat[l] != 7) return false;
    if (nbonds[l] != 3) return false;
    int nc = 0, nh = 0, nco = 0;
    for (int i = 1; i <= 3; ++i) {
        if (nat[ibonds[i][l]] == 6) {
            ++nc;
            int k = ibonds[i][l];
            if (nbonds[k] == 3) {
                for (int j = 1; j <= 3; ++j) {
                    if (nat[ibonds[j][k]] == 8) {
                        if (nbonds[ibonds[j][k]] == 1) {
                            ++nco;
                        } else {
                            return false;  // N-C-O-R
                        }
                    }
                }
            } else if (nbonds[k] == 4) {
                for (int j = 1; j <= 4; ++j) {
                    if (nat[ibonds[j][k]] == 8) return false;
                }
            }
        }
        if (nat[ibonds[i][l]] == 1) ++nh;
    }
    if (nc != 2 || nh != 1 || nco != 1) return false;
    return true;
}

void atomrs(std::vector<int>& lused, std::vector<bool>& ioptl, int& ires,
            int n1, int io, int uni_res, bool first_res) {
    using namespace MOZYME_C;
using namespace molkst_C;
    std::vector<int> inres(natomr + 1, 0);
    std::vector<int> npack(108, 0);
    std::vector<int> live(501, 0);
    int extra_atoms1 = 0, extra_atoms2 = 0;
    (void)io;

    int ninres = 0;
    int ihcr = nbackb[0], ico = nbackb[1], jofco = nbackb[2];
    jatom = nbackb[3];

    int nlive = nbonds[iatom];
    int j3 = 0, j4 = 0;
    for (int i2 = 1; i2 <= nlive; ++i2) {
        int j2 = ibonds[i2][iatom];
        if (nat[j2] == 1 || j2 == ihcr) {
            ++j3;
            live[j3] = j2;
        } else {
            double sum = 0.0;
            for (int i = 1; i <= nlive; ++i) {
                if (nat[ibonds[i][iatom]] == 1) {
                    for (int j = 1; j <= nbonds[j2]; ++j) {
                        if (nat[ibonds[j][j2]] == 8) {
                            dihed(coord, ibonds[i][iatom], iatom, j2, ibonds[j][j2], sum);
                        }
                    }
                }
            }
            if (sum < funcon_C::pi * 0.5 && sum > -funcon_C::pi * 0.5) {
                j4 = std::min(2, j4 + 1);
                if (j4 == 1) extra_atoms1 = j2; else extra_atoms2 = j2;
            }
        }
    }
    nlive = j3;

    int l = 0, i, j, jj3 = 0, OXT = 0, k = 0, j2, j5, j6, j7, j10, i2, i3,
        j5b, j6b, j8b, j9b;
    bool okay = false;
    std::vector<bool> found(15, false);
    bool UNK = false;
    std::string loc_tyres = "   ", tmp;

    // outer_loop: breadth-first walk of residue atoms.
    int outer_state = 0;  // 0: running; 1: break out after error
    while (true) {
        l = live[1];
        if (l == 0) goto walk_done;
        if (l == jatom || ioptl[l] || l == iatom || peptide_n(l)) {
            if (nlive == 0) goto walk_done;
            live[1] = live[nlive];
            --nlive;
        } else {
            ioptl[l] = true;
            ++ninres;
            if (ninres > natomr) {
                std::printf(" There are more than %d atoms in residue %d\n", natomr, ires);
                mopend("Too many atoms in residue");
            }
            inres[ninres] = l;
            if (nbonds[l] != 0) {
                for (i2 = 2; i2 <= nbonds[l]; ++i2) {
                    j = ibonds[i2][l];
                    bool dup = false;
                    for (i3 = 1; i3 <= nlive; ++i3)
                        if (live[i3] == j) { dup = true; break; }
                    if (dup) continue;
                    ++nlive;
                    if (nlive > 500) { outer_state = 1; break; }
                    live[nlive] = j;
                }
                if (outer_state == 1) break;
                live[1] = ibonds[1][l];
            } else {
                if (nlive == 0) goto walk_done;
                live[1] = live[nlive];
                --nlive;
            }
        }
    }
    if (outer_state == 1) {
        std::printf(" Number of live atoms in residue %d greater than 500\n", ires);
        for (i = 1; i <= numat; ++i)
            if (txtatm[i].find("UNK") != std::string::npos) break;
        mopend("More than 500 atoms in residue");
    }
walk_done:

    // Check atoms not counted twice.
    for (j = 1; j <= ninres; ++j)
        for (k = j + 1; k <= ninres; ++k)
            if (inres[j] == inres[k]) inres[j] = 0;
    npack.assign(108, 0);
    if (nat[iatom] == 7) npack[7] = 1;

    for (i2 = 1; i2 <= ninres; ++i2) {
        if (inres[i2] != 0) {
            j2 = inres[i2];
            for (j = 1; j <= nbonds[j2]; ++j) {
                j5 = nat[ibonds[j][j2]];
                if (!(j5 == 1 || (j5 >= 6 && j5 <= 8) || j5 == 16)) jj3 = j5;
            }
            npack[nat[j2]] = npack[nat[j2]] + 1;
        }
        lused[inres[i2]] = ires;
    }
    lused[iatom] = -ires;
    lused[ihcr] = -ires;
    if (ico > 0) lused[ico] = -ires;
    if (jofco != 0) lused[jofco] = -ires;

    // Determine residue name.
    OXT = 0;
    if (ico > 0) {
        k = 0; i2 = 0;
        for (j = 1; j <= nbonds[ico]; ++j) {
            j2 = ibonds[j][ico];
            if (nat[j2] == 8) {
                ++k;
                if (nbonds[j2] == 2) i2 = j2;
            }
        }
    }
    if (k == 2) { OXT = i2; npack[8] = npack[8] - 1; } else { OXT = 0; }
    int nh = npack[1];
    k = 0;
    loc_tyres = "   ";
    if (jj3 == 34) npack[16] = npack[16] + 1;

    if (jj3 == 0 || jj3 == 34) {
        for (loop = 1; loop <= mxeno; ++loop) {
            npack[1] = npack[6] - nxeno[0][loop-1];
            npack[2] = npack[7] - nxeno[1][loop-1];
            npack[3] = npack[8] - nxeno[2][loop-1];
            npack[4] = npack[16] - nxeno[3][loop-1];
            if (iatom == jatom && ico > 0) {
                j = 0;
                for (i = 1; i <= nbonds[ico]; ++i)
                    if (nat[ibonds[i][ico]] == 8) ++j;
                if (j == 2) npack[3] = std::max(0, npack[3] - 1);
            }

            for (j = 1; j <= size_mres - 1; ++j) {
                bool eq = true;
                for (k = 1; k <= 4; ++k)
                    if (npack[k] != mres[k][j]) { eq = false; break; }
                if (!eq) continue;
                if (j == 21 || j == 22) {
                    j2 = 0;
                    for (j = 1; j <= nbonds[ico]; ++j)
                        if (nat[ibonds[j][ico]] == 8) ++j2;
                    if (j2 != 2) continue;
                }
                if (j == 3) {
                    j4 = 0;
                    for (k = 1; k <= nbonds[ihcr]; ++k) {
                        l = 0; j2 = ibonds[k][ihcr];
                        if (nat[j2] == 6)
                            for (j = 1; j <= nbonds[j2]; ++j)
                                if (nat[ibonds[j][j2]] == 6) ++l;
                        j4 = std::max(j4, (int)l);
                    }
                    if (j4 == 2) {
                        okay = false;
                        for (k = 1; k <= nbonds[ihcr]; ++k) {
                            if (nat[ibonds[k][ihcr]] == 7) {
                                j2 = ibonds[k][ihcr];
                                for (j = 1; j <= nbonds[j2]; ++j) {
                                    if (nat[ibonds[j][j2]] == 6 && ibonds[j][j2] != ihcr) {
                                        j4 = ibonds[j][j2];
                                        for (j5b = 1; j5b <= nbonds[j4]; ++j5b) {
                                            if (nat[ibonds[j5b][j4]] == 6) {
                                                j6b = ibonds[j5b][j4];
                                                for (j7 = 1; j7 <= nbonds[j6b]; ++j7) {
                                                    if (nat[ibonds[j7][j6b]] == 6) {
                                                        j8b = ibonds[j7][j6b];
                                                        for (j9b = 1; j9b <= nbonds[j8b]; ++j9b)
                                                            for (j10 = 1; j10 <= nbonds[ihcr]; ++j10)
                                                                if (ibonds[j9b][j8b] == ibonds[j10][ihcr]) {
                                                                    okay = true; break;
                                                                }
                                                        if (okay) break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (okay) break;
                        }
                        k = okay ? 20 : 0;
                        if (!okay) continue;
                    }
                }
                k = j;
                goto matched_1200;
            }
            j = 0;
            if (iatom == jatom) {
                if (npack[3] != 0) { k = 1; npack[3] = npack[3] + 1; }
            } else if (jatom == 0) {
                k = 1; npack[3] = npack[3] - 1;
            }
            if (k != 0) {
                for (k = 1; k <= size_mres - 1; ++k) {
                    bool eq = true;
                    for (j = 1; j <= 4; ++j)
                        if (npack[j] != mres[j][k]) { eq = false; break; }
                    if (!eq) continue;
                    if (k == 21 || k == 22) {
                        j2 = 0;
                        for (j = 1; j <= nbonds[ico]; ++j)
                            if (nat[ibonds[j][ico]] == 8) ++j2;
                        if (j2 != 2) { loc_tyres = "HYP"; continue; }
                    }
                    if (k == 3) {
                        j4 = 0;
                        for (j = 1; j <= nbonds[ihcr]; ++j) {
                            l = 0; j2 = ibonds[j][ihcr];
                            for (j = 1; j <= nbonds[j2]; ++j)
                                if (nat[ibonds[j][j2]] == 6) ++l;
                            j4 = std::max(j4, (int)l);
                        }
                        if (j4 != 3) continue;
                    }
                    goto matched_1200;
                }
            }
        }
    } else {
        for (j = 1; j <= 10; ++j)
            if (atom_names[jj3][j - 1] != ' ') break;
        for (j6 = 1; j6 <= numat; ++j6)
            if (std::fabs(coord[0][ihcr] - coorda[0][j6]) < 5.0e-4)
                if (std::fabs(coord[1][ihcr] - coorda[1][j6]) < 5.0e-4)
                    if (std::fabs(coord[2][ihcr] - coorda[2][j6]) < 5.0e-4) break;
        for (j7 = 1; j7 <= n_once; ++j7)
            if (once[j7] == j6) break;
        if (j7 > n_once) {
            n_once = j7; once[n_once] = j6;
        }
        loop = 1;
        goto matched_1200;
    }
    if (npack[6] == 0 && npack[7] == 1 && npack[8] == 0 && npack[16] == 0) {
        ioptl[jatom] = true;
        --ires;
        return;
    }
matched_1200:
    if (k == 0) k = size_mres;
    if (k == 4) {
        for (i2 = 1; i2 <= nbonds[iatom]; ++i2) {
            j = ibonds[i2][iatom];
            if (nat[j] == 6) {
                for (j2 = 1; j2 <= nbonds[j]; ++j2) {
                    j3 = ibonds[j2][j];
                    if (nat[j3] == 6) {
                        l = 0;
                        for (j = 1; j <= nbonds[j3]; ++j) {
                            j5 = ibonds[j][j3];
                            if (nat[j5] == 6) ++l;
                            if (nat[j5] != 1 && nat[j5] != 6) break;
                        }
                        if (l == 3) k = 5;
                    }
                }
            }
        }
    } else if (k == 3) {
        l = 0;
        for (j = 1; j <= nbonds[iatom]; ++j)
            if (nat[ibonds[j][iatom]] != 1) ++l;
        if ((l == 3 && iatom != n1) || (l == 2 && iatom == n1)) k = 20;
    }
    if (nat[n1] == 1) {
        ++ninres;
        inres[ninres] = n1;
    }
    if (loc_tyres == "   ") loc_tyres = tyres[k];
    if (k == 19 && jj3 == 34) loc_tyres = "MSE";
    if (k == 23) {
        okay = false;
        if (npack[6] == 2 && npack[7] == 0 && npack[8] == 1 && npack[16] == 0) {
            loc_tyres = "ACE"; line = "Acetyl group"; okay = true;
        }
        if (npack[6] == 1 && npack[7] == 1 && npack[8] == 0 && npack[16] == 0) {
            loc_tyres = "CH3"; line = "Methyl group"; okay = true;
        }
        if (okay && keywrd.find(" RESID") != std::string::npos) {
            ++ncomments;
        }
    }
    if (txeno[loop].size() >= 3) loc_tyres = txeno[loop].substr(0, 3);
    for (j = 1; j <= ninres; ++j) {
        l = inres[j];
        if (l != 0)
            txtatm[l] = "ATOM " + std::to_string(l) + "  " + cap_elemnt[nat[l]] +
                        loc_tyres + " " + std::to_string(ires);
    }
    res_start[uni_res] = iatom;
    greek(jatom);
    if (k == 20) {
        for (i = 1; i <= ninres; ++i)
            if (inres[i] == jatom) break;
        if (i == ninres + 1 && jatom > 0) txtatm[jatom] = "   ";
    }
    if (txtatm[iatom].substr(0, 15).find_first_not_of(' ') != std::string::npos) {
        allres[ires] = loc_tyres + " ";
        return;
    }
    if (ires == 1 && txtatm[ihcr].size() >= 20 &&
        txtatm[ihcr].substr(17, 3) == "PRO") {
        txtatm[n1] = std::to_string(n1) + " PRO   1";
    }
    if (OXT != 0) txtatm[OXT] = txtatm[OXT].substr(0, 14) + "XT" + txtatm[OXT].substr(16);

    // Label sanity check.
    j2 = 0;
    for (j2 = 1; j2 <= 20; ++j2)
        if (afn[j2] == txtatm[inres[1]].substr(17, 3)) break;
    std::fill(found.begin(), found.end(), false);
    UNK = false;
    if (j2 < 21) {
        for (i = 1; i <= ninres; ++i) {
            l = inres[i];
            line = elemnt[nat[l]] + txtatm[l].substr(14, 1) + txtatm[l].substr(15, 1);
            for (j = 1; j <= n_add[j2]; ++j)
                if (line.substr(0, 4) == atomname[j2][j]) found[j] = true;
        }
        for (i = 5; i <= n_add[j2]; ++i)
            if (!found[i]) break;
        if (j2 != 8) {
            l = 0;
            for (j = 1; j <= nbonds[iatom]; ++j)
                if (nat[ibonds[j][iatom]] > 1) ++l;
            if (l > 2) i = 0;
        }
        if (i <= n_add[j2]) {
            loc_tyres = "UNK";
            for (i = 1; i <= ninres; ++i)
                txtatm[inres[i]] = txtatm[inres[i]].substr(0, 17) + "UNK";
            UNK = true;
        }
    }

    // Extra atoms (proline sidechain) handling.
    int ea1 = extra_atoms1, ea2 = extra_atoms2;
    if (ea1 > 0) {
        live[1] = ea1; live[2] = ea2;
        nlive = 1;
        if (live[2] != 0) nlive = 2;
        ninres = 0;
        while (true) {
            l = live[1];
            if (l == jatom || ioptl[l] || l == iatom) {
                if (nlive == 0) break;
                live[1] = live[nlive];
                --nlive;
            } else {
                ioptl[l] = true;
                ++ninres;
                if (ninres > natomr) mopend("Too many atoms in residue");
                inres[ninres] = l;
                if (nbonds[l] != 0) {
                    for (i2 = 2; i2 <= nbonds[l]; ++i2) {
                        j = ibonds[i2][l];
                        bool dup = false;
                        for (i3 = 1; i3 <= nlive; ++i3)
                            if (live[i3] == j) { dup = true; break; }
                        if (dup) continue;
                        ++nlive;
                        if (nlive > 500) break;
                        live[nlive] = j;
                    }
                    live[1] = ibonds[1][l];
                } else {
                    if (nlive == 0) break;
                    live[1] = live[nlive];
                    --nlive;
                }
            }
        }
        tmp = "UNK";
        npack.assign(108, 0);
        for (j = 1; j <= ninres; ++j)
            ++npack[nat[inres[j]]];
        okay = false;
        if (npack[6] == 1) {
            if (npack[7] == 1 && npack[8] == 0 && npack[16] == 0) {
                tmp = "CH3"; line = "Methyl group"; okay = true;
            }
        } else if (npack[6] == 2) {
            if (npack[7] == 0 && npack[8] == 1 && npack[16] == 0) {
                tmp = "ACE"; line = "Acetyl group"; okay = true;
            }
        } else if (npack[6] == 4) {
            if (npack[7] == 0 && npack[8] == 0 && npack[16] == 0) {
                tmp = "TBU"; line = "Tertiary butyl group"; okay = true;
            }
        }
        if (okay && keywrd.find(" RESID") != std::string::npos) ++ncomments;
        for (j = 1; j <= ninres; ++j) {
            l = inres[j];
            if (l != 0)
                txtatm[l] = "ATOM " + std::to_string(l) + "  " + cap_elemnt[nat[l]] +
                            tmp + " " + std::to_string(ires);
        }
    }
    txtatm[iatom] = "ATOM " + std::to_string(iatom) + "  " + cap_elemnt[nat[iatom]] +
                    loc_tyres + " " + std::to_string(ires);
    if (keywrd.find(" ADD-H") != std::string::npos) return;

    i = std::abs(mres[5][k] - nh);
    if (ico != 0 && !UNK && i > 1 && !first_res && loop == 1) {
        if (nbonds[ico] != 4 || nbonds[jofco] != 2 || i != 3) {
            if (odd_h && k != size_mres) {
                odd_h = false;
            }
        }
    }
    if (ires > maxres) mopend("Too many residues");
    allres[ires] = loc_tyres + " ";
    for (i = 1; i <= numat; ++i)
        if (ions[i] == 1 && std::abs(lused[i]) == ires) allres[ires][3] = '+';
    j = 0; l = 0;
    for (i = 1; i <= numat; ++i) {
        if (ions[i] > 0 && std::abs(lused[i]) == ires) j += ions[i];
        if (ions[i] < 0 && std::abs(lused[i]) == ires) l += ions[i];
    }
    if (j != 0 || l != 0) {
        if (j + l == 0) {
            if (j == 1) std::printf(" Residue:  '%s%d' is Zwitterionic\n",
                                     allres[ires].c_str(), ires);
        } else if (j > 0) {
            allres[ires][3] = '+';
        } else {
            allres[ires][3] = '-';
        }
    }
}


// Pointer wrapper for the names.F90 driver (1-based Fortran-style arrays).
void atomrs(int* lused, bool* ioptl, int ires, int n1, int io,
            int uni_res, bool first_res) {
    std::vector<int> lused_v(natomr + 1, 0);
    std::vector<bool> ioptl_v(natomr + 1, false);
    for (int i = 1; i <= natomr; ++i) { lused_v[i] = lused[i]; ioptl_v[i] = ioptl[i]; }
    atomrs(lused_v, ioptl_v, ires, n1, io, uni_res, first_res);
    for (int i = 1; i <= natomr; ++i) { lused[i] = lused_v[i]; ioptl[i] = ioptl_v[i]; }
}
