// moldat_helpers.cpp — C++ translation of moldat.F90 internal subprograms:
//   setcup      (lines 897-1035):  CUTOFP and unit-cell counts for periodic runs,
//   setup_nhco  (lines 1250-1323): MM correction to -(C=O)-(NH)- linkages.
// 1-based Fortran indexing is kept for all module arrays.
#include "moldat_helpers.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "distance.h"
#include "molmec_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "reada.h"
#include "volume.h"

using namespace molkst_C;        // cutofp, id, keywrd, l1u, l2u, l3u, l123, l11, l21, l31, line
using namespace common_arrays_C; // tvec, nat
using namespace chanel_C;        // iw
using namespace molmec_C;        // nnhco, nhco, htype

void setcup() {
    cutofp = 1.e10;
    l1u = 0;
    l2u = 0;
    l3u = 0;
    l123 = 1;
    if (id == 0) return;
    std::string::size_type i = keywrd.find(" CUTOFP");
    if (i != std::string::npos) {
        cutofp = reada(keywrd, (int)i + 8);  // Fortran Index is 1-based -> C++ find + 1
    } else {
        cutofp = 30.0;
    }
    if (id == 1) {  // Polymer case
        double tv1 = std::sqrt(tvec[1][1] * tvec[1][1] + tvec[2][1] * tvec[2][1] +
                               tvec[3][1] * tvec[3][1]);
        if (tv1 < 1.0) {
            line = "  Length of translation vector is too small.";
            std::printf("%s\n", line.c_str());
            mopend(line);
            return;
        }
        l1u = (int)(cutofp * 4.0 / 3.0 / tv1) + 1;
    } else if (id == 2) {  // Layer system
        double r1 = std::sqrt(tvec[1][1] * tvec[1][1] + tvec[2][1] * tvec[2][1] +
                              tvec[3][1] * tvec[3][1]);
        double r2 = std::sqrt(tvec[1][2] * tvec[1][2] + tvec[2][2] * tvec[2][2] +
                              tvec[3][2] * tvec[3][2]);
        if (r1 < 1.0) {
            line = "  Length of first translation vector is too small.";
            std::printf("%s\n", line.c_str());
            mopend(line);
            return;
        }
        // NB: the Fortran source tests r1 twice here (source quirk, preserved).
        if (r1 < 1.0) {
            line = "  Length of second translation vector is too small.";
            std::printf("%s\n", line.c_str());
            mopend(line);
            return;
        }
        double r12 = std::sqrt((tvec[1][2] - tvec[1][1]) * (tvec[1][2] - tvec[1][1]) +
                               (tvec[2][2] - tvec[2][1]) * (tvec[2][2] - tvec[2][1]) +
                               (tvec[3][2] - tvec[3][1]) * (tvec[3][2] - tvec[3][1]));
        double tv1 = r1 * std::sin(std::acos((r1 * r1 + r2 * r2 - r12 * r12) / (2 * r1 * r2)));
        double tv2 = r2 * std::sin(std::acos((r1 * r1 + r2 * r2 - r12 * r12) / (2 * r1 * r2)));
        l1u = (int)(cutofp * 4.0 / 3.0 / tv1) + 1;
        l2u = (int)(cutofp * 4.0 / 3.0 / tv2) + 1;
    } else {  // Solid-state (three-dimensional crystal)
        double r1 = std::sqrt(tvec[1][1] * tvec[1][1] + tvec[2][1] * tvec[2][1] +
                              tvec[3][1] * tvec[3][1]);
        double r2 = std::sqrt(tvec[1][2] * tvec[1][2] + tvec[2][2] * tvec[2][2] +
                              tvec[3][2] * tvec[3][2]);
        double r3 = std::sqrt(tvec[1][3] * tvec[1][3] + tvec[2][3] * tvec[2][3] +
                              tvec[3][3] * tvec[3][3]);
        if (r1 < 1.0 || r2 < 1.0 || r3 < 1.0) {
            if (r1 < 1.0) {
                line = "  Length of first translation vector is too small.";
            } else if (r2 < 1.0) {
                line = "  Length of second translation vector is too small.";
            } else {
                line = "  Length of third translation vector is too small.";
            }
            std::printf("%s\n", line.c_str());
            mopend(line);
            return;
        }
        double r12 = std::sqrt((tvec[1][2] - tvec[1][1]) * (tvec[1][2] - tvec[1][1]) +
                               (tvec[2][2] - tvec[2][1]) * (tvec[2][2] - tvec[2][1]) +
                               (tvec[3][2] - tvec[3][1]) * (tvec[3][2] - tvec[3][1]));
        double r13 = std::sqrt((tvec[1][3] - tvec[1][1]) * (tvec[1][3] - tvec[1][1]) +
                               (tvec[2][3] - tvec[2][1]) * (tvec[2][3] - tvec[2][1]) +
                               (tvec[3][3] - tvec[3][1]) * (tvec[3][3] - tvec[3][1]));
        double r23 = std::sqrt((tvec[1][3] - tvec[1][2]) * (tvec[1][3] - tvec[1][2]) +
                               (tvec[2][3] - tvec[2][2]) * (tvec[2][3] - tvec[2][2]) +
                               (tvec[3][3] - tvec[3][2]) * (tvec[3][3] - tvec[3][2]));
        double area12 = r1 * r2 * std::sin(std::acos((r1 * r1 + r2 * r2 - r12 * r12) / (2 * r1 * r2)));
        double area13 = r1 * r3 * std::sin(std::acos((r1 * r1 + r3 * r3 - r13 * r13) / (2 * r1 * r3)));
        double area23 = r2 * r3 * std::sin(std::acos((r2 * r2 + r3 * r3 - r23 * r23) / (2 * r2 * r3)));
        double v[9] = {tvec[1][1], tvec[2][1], tvec[3][1],
                       tvec[1][2], tvec[2][2], tvec[3][2],
                       tvec[1][3], tvec[2][3], tvec[3][3]};
        double vol = volume(v, id);
        if (vol < 1.0) {
            std::printf("\n\n%10sVolume of unit cell unreasonably small: %9.6f Cubic Angstroms\n",
                        "", vol);
            mopend("Volume of unit cell unreasonably small");
            return;
        }
        double tv1 = vol / area23;
        double tv2 = vol / area13;
        double tv3 = vol / area12;
        const double sum = 4.0;
        if ((tv1 < sum || tv2 < sum || tv3 < sum) && keywrd.find(" GEO-OK") == std::string::npos) {
            std::printf("\n\n%10sTranslation vector length 1: %9.3f Angstroms\n", "", r1);
            std::printf("%10sTranslation vector length 2: %9.3f Angstroms\n", "", r2);
            std::printf("%10sTranslation vector length 3: %9.3f Angstroms\n", "", r3);
            if (tv1 < sum)
                std::printf("%10sDistance between faces 2 and 3 is unreasonably small: %6.3f Angstroms, min: %6.3f\n",
                            "", tv1, sum);
            if (tv2 < sum)
                std::printf("%10sDistance between faces 1 and 3 is unreasonably small: %6.3f Angstroms, min: %6.3f\n",
                            "", tv2, sum);
            if (tv3 < sum)
                std::printf("%10sDistance between faces 1 and 2 is unreasonably small: %6.3f Angstroms, min: %6.3f\n",
                            "", tv3, sum);
            if (tv1 < sum || tv2 < sum || tv3 < sum) {
                mopend("One or more translation vectors are unreasonably small");
                std::printf("%10sTo over-ride this safety check, add keyword 'GEO-OK'\n", "");
                return;
            }
        }
        l1u = (int)((cutofp * 4.0 / 3.0) / tv1) + 1;
        l2u = (int)((cutofp * 4.0 / 3.0) / tv2) + 1;
        l3u = (int)((cutofp * 4.0 / 3.0) / tv3) + 1;
    }
    l123 = (2 * l1u + 1) * (2 * l2u + 1) * (2 * l3u + 1);
    l11 = std::min(l1u, 1);
    l21 = std::min(l2u, 1);
    l31 = std::min(l3u, 1);
}

void setup_nhco(int& ii) {
    nnhco = 0;
    htype = 0.0;
    if (method_mndo) htype = 6.1737;
    if (method_am1) htype = 3.3191;
    if (method_pm3) htype = 7.1853;
    if (method_rm1) htype = 2.4127;
    if (method_pm7) htype = 3.1595;
    if (method_pm6) htype = 2.5000;
    ii = 0;
    if (keywrd.find("NOMM") != std::string::npos) ii = 1;
    // Identify O=C-N-H systems via the interatomic distances matrix.
    for (int j = 1; j <= numat; ++j) {
        if (nat[j] != 6) continue;
        bool skip = false;
        for (int i = 1; i <= numat && !skip; ++i) {
            if (nat[i] != 8) continue;
            if (distance(i, j) > 1.3) continue;
            for (int k = 1; k <= numat && !skip; ++k) {
                if (nat[k] != 7) continue;
                if (distance(k, j) > 1.6) continue;
                for (int l = 1; l <= numat && !skip; ++l) {
                    if (nat[l] != 1) continue;
                    if (distance(k, l) > 1.3) continue;
                    // H-N-C=O found (atoms l-k-j-i).  Now search the atom
                    // attached to nitrogen to specify the X-N-C=O system.
                    for (int m = 1; m <= numat && !skip; ++m) {
                        if (m == k || m == l || m == j) continue;
                        if (distance(m, k) > 1.7) continue;
                        bool already = false;
                        for (int jj = 1; jj <= nnhco; jj += 2) {
                            if (nhco[2][jj] == k) { already = true; break; }
                        }
                        if (already) continue;
                        ++nnhco;
                        nhco[0][nnhco] = i;
                        nhco[1][nnhco] = j;
                        nhco[2][nnhco] = k;
                        nhco[3][nnhco] = m;
                        ++nnhco;
                        nhco[0][nnhco] = i;
                        nhco[1][nnhco] = j;
                        nhco[2][nnhco] = k;
                        nhco[3][nnhco] = l;
                        if (ii != 0) {
                            ii += 2;
                            nnhco -= 2;
                        }
                        skip = true;  // Fortran: cycle l230
                    }
                }
            }
        }
    }
}
