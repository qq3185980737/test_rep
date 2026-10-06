// grid.cpp — C++ translation of MOPAC 2016 "grid.F90".
// GRID: run a 2-D grid of SCF/geometry optimizations over two internal
// coordinates and store the lowest-energy geometry for each grid point.
#include "grid.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "maps_C.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace chanel_C { extern int iw, iw0, ires, iarc, iump; extern std::string archive_fn, ump_fn, restart_fn; }
namespace common_arrays_C {
extern std::vector<std::vector<double>> geo, geoa;
extern std::vector<double> xparam, pa, pb, p;
extern std::vector<int> na, nb, nc;
}
namespace molkst_C {
extern int nvar, natoms, norbs, numat, mpack;
extern std::string keywrd, line;
extern double tleft;
extern bool gui, uhf, moperr;
}
namespace maps_C {
extern double rxn_coord1, rxn_coord2;
extern int ione, ijlp, ilp, jlp, jlp1;
extern std::vector<double> surf;
extern int lpara1, latom1, lpara2, latom2;
}

extern double reada(const std::string&, int);
extern double second(int);
extern void ef(std::vector<double>&, double&);
extern void flepo(std::vector<double>&, int, double&);
extern void geout(int);
extern void pdbout(int);
extern void wrttxt(int);
extern void bonds();
extern void to_screen(const std::string&);
extern void mopend(const std::string&);

void grid() {
    using namespace chanel_C;
    using namespace common_arrays_C;
    using namespace molkst_C;
    using namespace maps_C;

    int npts1 = 11, npts2 = 11, maxcyc = 100000;
    int i, iloop, jloop, j, k, ij, l, iw00, percent = 0, max_count = 0, big_loop, loop;
    double step1 = 0.0, step2 = 0.0, degree = 57.29577951308232, c1, c2, cputot = 0.0;
    double escf = 0.0, cpu1, cpu2, cpu3, geo11, geo22, sum = 0.0;
    std::string formt = "14.3";
    bool restrt, useef, opend, minimize_energy_in_grid = false, use_p = false;
    std::vector<double> all_points1, all_points2, all_geo, all_nabc_flat;
    std::vector<std::vector<double>> xy, surfac, all_pa, all_pb;
    std::vector<int> all_nabc;
    char txt;

    geoa.assign(4, std::vector<double>((size_t)natoms + 1, 0.0));   // 1-based rows 1..3
    for (i = 1; i <= 3; ++i)
        for (j = 1; j <= natoms; ++j) geoa[i][j] = geo[i][j];
    iw00 = iw0;
    iw0 = -1;  // While in "grid", do not print working in other subroutines.
    useef = (keywrd.find(" BFGS") == std::string::npos || nvar < 2);
    i = (int)keywrd.find("STEP1") + 6;
    step1 = reada(keywrd, i);
    i = (int)keywrd.find("STEP2") + 6;
    step2 = reada(keywrd, i);
    minimize_energy_in_grid = (keywrd.find(" SMOOTH") != std::string::npos);
    if (keywrd.find(" BIGCYCLES") != std::string::npos) {
        maxcyc = (int)std::lround(reada(keywrd, (int)keywrd.find(" BIGCYCLES")));
    }
    if (keywrd.find("POINT1") != std::string::npos) {
        npts1 = (int)std::lround(std::fabs(reada(keywrd, (int)keywrd.find("POINT1") + 7)));
    }
    if (keywrd.find("POINT2") != std::string::npos) {
        npts2 = (int)std::lround(std::fabs(reada(keywrd, (int)keywrd.find("POINT2") + 7)));
    }
    restrt = (keywrd.find(" RESTART") != std::string::npos);
    surf.clear();
    surf.resize((size_t)npts1 * npts2 + 1, 1.0e9);
    surfac.assign((size_t)npts2 + 1, std::vector<double>((size_t)npts1 + 1, 0.0));
    xy.assign(3, std::vector<double>(4 * (size_t)npts1 * npts2 + 1, 0.0));

    // First pass: left-to-right, snake over rows.
    k = 0;
    for (i = 1; i <= npts1; ++i)
        for (j = 1; j <= npts2; ++j) {
            ++k;
            xy[1][k] = i - 1;
            if (i % 2 == 1)
                xy[2][k] = j - 1;
            else
                xy[2][k] = npts2 - j;
        }
    ij = npts1 * npts2;
    // Second pass: reverse of first.
    for (i = 1; i <= ij; ++i) {
        ++k;
        xy[1][k] = xy[1][ij - i + 1];
        xy[2][k] = xy[2][ij - i + 1];
    }
    // Third pass: top-to-bottom, snake over columns.
    for (j = 1; j <= npts2; ++j)
        for (i = 1; i <= npts1; ++i) {
            ++k;
            xy[2][k] = j - 1;
            if (j % 2 == 1)
                xy[1][k] = i - 1;
            else
                xy[1][k] = npts1 - i;
        }
    // Fourth pass: reverse of third.
    ij = npts1 * npts2;
    l = 3 * ij + 1;
    for (i = 1; i <= ij; ++i) {
        ++k;
        xy[1][k] = xy[1][l - i];
        xy[2][k] = xy[2][l - i];
    }

    if (lpara1 != 1) step1 = step1 / degree;
    if (lpara2 != 1) step2 = step2 / degree;
    geo11 = geo[lpara1][latom1];
    geo22 = geo[lpara2][latom2];
    max_count = 0;
    i = npts1 * npts2;
    j = (uhf ? mpack : 1);
    // all_* grow dynamically with max_count (vector).
    all_points1.assign((size_t)i + 1, 0.0);
    all_points2.assign((size_t)i + 1, 0.0);
    all_geo.assign((size_t)i * 3 * (size_t)natoms + 1, 0.0);
    all_nabc.assign((size_t)i * 3 * (size_t)natoms + 1, 0);
    all_pa.assign((size_t)i + 1, std::vector<double>((size_t)mpack + 1, 0.0));
    all_pb.assign((size_t)i + 1, std::vector<double>((size_t)j + 1, 0.0));
    use_p = (i == 0);
    surf.assign((size_t)npts1 * npts2 + 1, 1.0e9);
    if (lpara1 != 1 && na[latom1] > 0)
        c1 = degree;
    else
        c1 = 1.0;
    if (lpara2 != 1 && na[latom2] > 0)
        c2 = degree;
    else
        c2 = 1.0;

    if (restrt) {
        std::ifstream rf(restart_fn.c_str(), std::ios::binary);
        if (!rf.good()) {
            mopend("Restart file either does not exist or is not available for reading");
            return;
        }
        rf.seekg(0);
        if (nvar > 0) {
            int ii, jj;
            rf.read((char*)&ii, sizeof(int)); rf.read((char*)&jj, sizeof(int));
            if ((norbs != jj || numat != ii) && (norbs != ii || numat != jj)) {
                mopend("Restart file read in does not match current data set");
                return;
            }
        }
        if (useef && nvar > 1) {
            for (i = 1; i <= 7; ++i) rf.read((char*)&txt, sizeof(char));
        } else if (nvar > 0) {
            for (i = 1; i <= 3; ++i) rf.read((char*)&txt, sizeof(char));
        }
        int i6[6];
        rf.read((char*)i6, 6 * sizeof(int));
        max_count = i6[0]; ijlp = i6[1]; ilp = i6[2]; jlp = i6[3]; jlp1 = i6[4]; ione = i6[5];
        rf.read((char*)&rxn_coord1, sizeof(double));
        rf.read((char*)&rxn_coord2, sizeof(double));
        all_geo.resize((size_t)max_count * 3 * (size_t)natoms + 1);
        all_nabc.resize((size_t)max_count * 3 * (size_t)natoms + 1);
        all_points1.resize((size_t)max_count + 1);
        all_points2.resize((size_t)max_count + 1);
        surf.resize((size_t)max_count + 1);
        for (i = 1; i <= max_count; ++i) {
            int base = (i - 1) * 3 * natoms;
            for (j = 1; j <= 3; ++j)
                for (int a = 1; a <= natoms; ++a)
                    rf.read((char*)&all_geo[(size_t)base + (j - 1) * natoms + a], sizeof(double));
            for (j = 1; j <= 3; ++j)
                for (int a = 1; a <= natoms; ++a)
                    rf.read((char*)&all_nabc[(size_t)base + (j - 1) * natoms + a], sizeof(int));
            rf.read((char*)&surf[i], sizeof(double));
            rf.read((char*)&all_points1[i], sizeof(double));
            rf.read((char*)&all_points2[i], sizeof(double));
        }
    } else {
        ijlp = 1;
        cputot = 0.0;
    }

    if (minimize_energy_in_grid)
        big_loop = 4 * npts1 * npts2;
    else
        big_loop = npts1 * npts2;
    if (iw00 > -1) to_screen("       FIRST VARIABLE   SECOND VARIABLE        FUNCTION      DONE LEFT");

    for (loop = ijlp; loop <= big_loop; ++loop) {
        if (loop >= maxcyc) tleft = -100.0;
        geo[lpara1][latom1] = geo11 + xy[1][loop] * step1;
        geo[lpara2][latom2] = geo22 + xy[2][loop] * step2;
        cpu1 = second(2);
        if (useef && nvar > 1)
            ef(xparam, escf);
        else
            flepo(xparam, nvar, escf);
        cpu2 = second(2);
        cpu3 = cpu2 - cpu1;
        cputot = cputot + cpu3;
        ++jlp;
        ++ijlp;

        // Find the point defined by the two variable coordinates.
        k = 0;
        for (i = 1; i <= max_count; ++i) {
            if (std::fabs(all_points1[i] - geo[lpara1][latom1]) < 1.0e-6 &&
                std::fabs(all_points2[i] - geo[lpara2][latom2]) < 1.0e-6) {
                k = i;
                break;
            }
        }
        if (k == 0) {
            ++max_count;
            k = max_count;
            if ((size_t)max_count >= all_points1.size()) {
                all_points1.resize((size_t)max_count + 16, 0.0);
                all_points2.resize((size_t)max_count + 16, 0.0);
                all_geo.resize((size_t)max_count * 3 * (size_t)natoms + 16 * 3 * (size_t)natoms, 0.0);
                all_nabc.resize((size_t)max_count * 3 * (size_t)natoms + 16 * 3 * (size_t)natoms, 0);
                all_pa.resize((size_t)max_count + 16, std::vector<double>((size_t)mpack + 1, 0.0));
                all_pb.resize((size_t)max_count + 16, std::vector<double>((size_t)(uhf ? mpack : 1) + 1, 0.0));
            }
            all_points1[k] = geo[lpara1][latom1];
            all_points2[k] = geo[lpara2][latom2];
        }

        if (escf + 1.0e-10 < surf[k]) {
            surf[k] = escf;
            int base = (k - 1) * 3 * natoms;
            for (i = 1; i <= natoms; ++i)
                for (j = 1; j <= 3; ++j) all_geo[(size_t)base + (j - 1) * natoms + i] = geo[j][i];
            for (i = 1; i <= natoms; ++i) {
                all_nabc[(size_t)base + (1 - 1) * natoms + i] = na[i];
                all_nabc[(size_t)base + (2 - 1) * natoms + i] = nb[i];
                all_nabc[(size_t)base + (3 - 1) * natoms + i] = nc[i];
            }
            if (gui) {
                if (use_p) {
                    for (ij = 1; ij <= mpack; ++ij) all_pa[k][ij] = pa[ij];
                    if (uhf)
                        for (ij = 1; ij <= mpack; ++ij) all_pb[k][ij] = pb[ij];
                }
            }
            char lbuf[160];
            std::snprintf(lbuf, sizeof(lbuf), " :%16.5f%16.5f%21.6f%10d%5d",
                          geo[lpara1][latom1] * c1, geo[lpara2][latom2] * c2, escf, loop, big_loop - loop);
            line = lbuf;
            if (iw00 > -1) {
                i = (int)std::lround((100.0 * loop) / big_loop);
                if (i != percent) {
                    percent = i;
                    std::snprintf(lbuf, sizeof(lbuf), "%4d%% of Grid Surface done", percent);
                    to_screen(lbuf);
                }
                to_screen(line);
            }
            std::fprintf(stdout, "%s\n", line.c_str());
            if (keywrd.find(" DEBUG") != std::string::npos) {
                if (keywrd.find(" PDBOUT") != std::string::npos)
                    pdbout(iw);
                else
                    geout(iw);
            }
        }
        if (tleft < 0.0 || moperr) {
            if (tleft < 0.0) {
                std::ofstream rf(restart_fn.c_str(), std::ios::binary | std::ios::trunc);
                if (rf.good()) {
                    int i6[6] = {max_count, ijlp, ilp, jlp, jlp1, ione};
                    rf.write((char*)i6, 6 * sizeof(int));
                    rf.write((char*)&rxn_coord1, sizeof(double));
                    rf.write((char*)&rxn_coord2, sizeof(double));
                    int base;
                    for (i = 1; i <= max_count; ++i) {
                        base = (i - 1) * 3 * natoms;
                        for (j = 1; j <= 3; ++j)
                            for (int a = 1; a <= natoms; ++a)
                                rf.write((char*)&all_geo[(size_t)base + (j - 1) * natoms + a], sizeof(double));
                        for (j = 1; j <= 3; ++j)
                            for (int a = 1; a <= natoms; ++a)
                                rf.write((char*)&all_nabc[(size_t)base + (j - 1) * natoms + a], sizeof(int));
                        rf.write((char*)&surf[i], sizeof(double));
                        rf.write((char*)&all_points1[i], sizeof(double));
                        rf.write((char*)&all_points2[i], sizeof(double));
                    }
                }
            }
            surf.clear();
            surfac.clear();
            return;
        }
    }

    if (minimize_energy_in_grid) std::fprintf(stdout, " Survey complete.  About to print points\n");
    for (iloop = 1; iloop <= npts1; ++iloop)
        for (jloop = 1; jloop <= npts2; ++jloop) {
            geo[lpara1][latom1] = geo11 + (iloop - 1) * step1;
            geo[lpara2][latom2] = geo22 + (jloop - 1) * step2;
            k = 0;
            for (i = 1; i <= max_count; ++i) {
                if (std::fabs(all_points1[i] - geo[lpara1][latom1]) < 1.0e-6 &&
                    std::fabs(all_points2[i] - geo[lpara2][latom2]) < 1.0e-6) {
                    k = i;
                    break;
                }
            }
            escf = surf[k];
            int base = (k - 1) * 3 * natoms;
            for (i = 1; i <= natoms; ++i)
                for (j = 1; j <= 3; ++j) geo[j][i] = all_geo[(size_t)base + (j - 1) * natoms + i];
            for (i = 1; i <= natoms; ++i) {
                na[i] = all_nabc[(size_t)base + (1 - 1) * natoms + i];
                nb[i] = all_nabc[(size_t)base + (2 - 1) * natoms + i];
                nc[i] = all_nabc[(size_t)base + (3 - 1) * natoms + i];
            }
            surfac[jloop][iloop] = escf;
            std::fprintf(stdout, "\n          FIRST VARIABLE   SECOND VARIABLE FUNCTION\n");
            std::fprintf(stdout, " :%16.5f%16.5f%16.6f\n",
                         geo[lpara1][latom1] * c1, geo[lpara2][latom2] * c2, escf);
            geout(iw);
            if (gui) {
                if (use_p) {
                    ij = 0;
                    for (i = 1; i <= norbs; ++i)
                        for (j = 1; j <= i; ++j) {
                            ++ij;
                            pa[ij] = all_pa[k][ij];
                        }
                    if (uhf) {
                        ij = 0;
                        for (i = 1; i <= norbs; ++i)
                            for (j = 1; j <= i; ++j) {
                                ++ij;
                                pb[ij] = all_pb[k][ij];
                            }
                    } else {
                        std::fill(pa.begin(), pa.end(), 0.0);
                        std::fill(pb.begin(), pb.end(), 0.0);
                    }
                    for (ij = 1; ij <= mpack; ++ij) p[ij] = pa[ij] + pb[ij];
                } else {
                    for (ij = 1; ij <= mpack; ++ij) p[ij] = 2 * pa[ij];
                }
                bonds();
            }
        }

    std::fprintf(stdout, "\n          HORIZONTAL: VARYING SECOND PARAMETER,\n          VERTICAL:   VARYING FIRST PARAMETER\n");
    std::fprintf(stdout, "\n          WHOLE OF GRID, SUITABLE FOR PLOTTING\n\n");

    // ARCHIVE
    std::ofstream af(archive_fn.c_str(), std::ios::out | std::ios::trunc);
    if (af.good()) {
        af << " ARCHIVE FILE FOR GRID CALCULATION\nGRID OF HEATS\n\n";
    }
    wrttxt(iarc);
    std::fprintf(stdout, "\n          TOTAL JOB TIME IN FLEPO : %10.3f\n", cputot);
    if (af.good()) af << "\n          TOTAL JOB TIME IN FLEPO : " << cputot << "\n";

    sum = 0.0;
    for (i = 1; i <= npts1; ++i)
        for (j = 1; j <= npts2; ++j)
            if (std::fabs(surfac[j][i]) > sum) sum = std::fabs(surfac[j][i]);
    if (sum > 0.99e7) formt = "14.3";
    else if (sum > 0.99e6) formt = "13.3";
    else if (sum > 0.99e5) formt = "12.3";
    else if (sum > 0.99e4) formt = "11.3";
    else if (sum > 0.99e3) formt = "10.3";
    else if (sum > 0.99e2) formt = "9.3";
    else if (sum > 0.99e1) formt = "8.3";

    std::fprintf(stdout, "\n");
    for (j = 1; j <= npts2; ++j)
        std::fprintf(stdout, ("%11." + formt + "f").c_str(), (geo22 + (j - 1) * step2) * c2);
    std::fprintf(stdout, "\n");
    if (af.good()) {
        af << "\n";
        for (j = 1; j <= npts2; ++j)
            af << " " << (geo22 + (j - 1) * step2) * c2;
        af << "\n";
    }
    std::ofstream uf(ump_fn.c_str(), std::ios::out | std::ios::trunc);
    for (i = 1; i <= npts1; ++i) {
        std::fprintf(stdout, "\n%10.3f", (geo11 + (i - 1) * step1) * c1);
        for (j = 1; j <= npts2; ++j) std::fprintf(stdout, (" %11." + formt + "f").c_str(), surfac[j][i]);
        std::fprintf(stdout, "\n");
        if (af.good()) {
            af << "\n" << (geo11 + (i - 1) * step1) * c1;
            for (j = 1; j <= npts2; ++j) af << " " << surfac[j][i];
            af << "\n";
        }
        if (uf.good()) {
            for (j = 1; j <= npts2; ++j)
                uf << (rxn_coord1 + step1 * (i - 1)) * c1 << " "
                   << (rxn_coord2 + step2 * (j - 1)) * c2 << " " << surfac[j][i] << "\n";
        }
    }
    if (uf.good()) uf.close();

    surf.clear();
    surfac.clear();
    iw0 = iw00;
}
