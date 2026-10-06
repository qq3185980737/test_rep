
// geo_ref.cpp — C++ translation of geo_ref.F90 (MOPAC 2016).
#include "geo_ref.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "big_swap.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "geo_diff.h"
#include "getgeo.h"
#include "getpdb.h"
#include "gmetry.h"
#include "molkst_C.h"
#include "upcase.h"
#include "web_message.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

// external hooks
double reada(const std::string&, int);            // 1-based position
void mopend(const std::string&);
void add_path(std::string&);
void geout(int);
void pdbout(int);
void wrt_diffs();
void analyze_h_bonds();
void dock(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, double&);
void web_message(int, const char*);

namespace {

// Fortran 1-based substring helpers
inline std::string substr1(const std::string& s, int a, int b) {
    if (a < 1) a = 1;
    if (b > (int)s.size()) b = (int)s.size();
    if (a > b) return std::string();
    return s.substr(a - 1, b - a + 1);
}
inline std::string substr1(const std::string& s, int a) {
    if (a < 1) a = 1;
    if (a > (int)s.size()) return std::string();
    return s.substr(a - 1);
}
// Fortran index(s, sub): 1-based position or 0.
inline int index1(const std::string& s, const std::string& sub) {
    std::size_t p = s.find(sub);
    return (p == std::string::npos) ? 0 : (int)p + 1;
}
inline std::string ltrim(const std::string& s) {
    std::size_t p = s.find_first_not_of(' ');
    return (p == std::string::npos) ? std::string() : s.substr(p);
}
inline std::string rtrim(const std::string& s) {
    std::size_t p = s.find_last_not_of(' ');
    return (p == std::string::npos) ? std::string() : s.substr(0, p + 1);
}
inline std::string trim(const std::string& s) { return rtrim(ltrim(s)); }

}

void compare_txtatm(bool& bug, bool& any_bug);

void geo_ref() {
    int i, j, k, l, ii, jj, i4, j4, k4, iquit, numat_dat, numat_ref;
    std::vector<int> map_atoms_A, atom_no;
    std::vector<std::string> tmp_txt, diffs;
    std::vector<double> tmp_atmass;
    std::vector<std::vector<double>> tmp_geoa;
    double dum1, dum2, sum, rms, rms_min, sum1, sum2, sum3, toler, xmin, sum4;
    bool intern = true, exists, bug, any_bug, swap, pdb, html, first, let;
    std::vector<bool> same, ok;
    std::string line_1, num;
    std::vector<std::string> refkey_ref(7, " ");  // 1-based 1..6

    //   For Geo-Ref to work, some very specific conditions must be satisfied.
    //   So before attempting a GEO_REF calculation, check that the data are okay
    for (j = 2; j <= natoms; ++j)
        if (na[j] != 0) break;
    if (nvar != 3 * numat + 3 * id) {
        std::fprintf(stdout, "%s\n", " GEO_REF requires all parameters to be optimized");
        std::fprintf(stdout, "%s\n", " (The easiest way to do this is to add \"OPT\" to the set of keywords.)");
        mopend("GEO_REF requires all parameters to be optimized");
    } else if (j <= natoms) {
        std::fprintf(stdout, "%s\n", " GEO_REF requires Cartesian coordinates");
        std::fprintf(stdout, "%s\n", " (The easiest way to do this is to add \"XYZ\" to the set of keywords.)");
        mopend("GEO_REF requires Cartesian coordinates");
    }
    if (moperr) return;
    geoa.assign(4, std::vector<double>(natoms + 301, 0.0));  // 1-based, generous safety factor
    c.assign(4, std::vector<double>(natoms + 301, 0.0));
    for (ii = 1; ii <= 6; ++ii) {
        line = " " + rtrim(refkey[ii]);
        upcase(line, (int)line.size());
        i = index1(line, " GEO_REF");
        if (i != 0) break;
    }
    j = index1(substr1(refkey[ii], i + 10), " ") + i + 8;
    if (index1(substr1(line, i, j), "\"") + index1(substr1(line, i, j), "'") == 0) {
        line = " File name after GEO_REF must be in quotation marks.";
        mopend(trim(line));
        return;
    }
    j = index1(substr1(refkey[ii], i + 10), "\"") + index1(substr1(refkey[ii], i + 10), "'");
    if (j == 0) {
        line = " File name after GEO_REF must end with a quotation mark.";
        mopend(trim(line));
        return;
    }
    j = j + i + 8;
    line = substr1(refkey[ii], i + 9, j);
    line_1 = trim(line);
    upcase(line_1, (int)line_1.size());
    if (line_1 == "SELF") {
        line = output_fn.substr(0, (int)output_fn.size() - 3) + "mop";
        refkey[ii] = substr1(refkey[ii], 1, i - 1) + "geo_ref=\"" + trim(line) +
                     substr1(refkey[ii], j + 1);
    }
    geo_ref_name = trim(line);
    k = index1(keywrd, " GEO_DAT") + 10;
    if (k > 10) {
        for (ii = 1; ii <= 6; ++ii) {
            line = " " + rtrim(refkey[ii]);
            upcase(line, (int)line.size());
            i = index1(line, " GEO_DAT");
            if (i != 0) break;
        }
        j = index1(substr1(refkey[ii], i + 9), "\" ") + i + 7;
        geo_dat_name = substr1(refkey[ii], i + 9, j);
    } else {
        geo_dat_name = trim(job_fn);
    }
    if (index1(keywrd, " 0SCF") != 0 && index1(keywrd, " HTML ") != 0) {
        line = trim(geo_dat_name);
        upcase(line, (int)line.size());
        line = line.substr((int)line.size() - 4);
        if (line == ".PDB") mopend("GEO_DAT file cannot be a PDB file (It would be over-written)");
        line = trim(geo_ref_name);
        upcase(line, (int)line.size());
        line = line.substr((int)line.size() - 4);
        if (line == ".PDB") mopend("GEO_REF file cannot be a PDB file (It would be over-written)");
        if (moperr) return;
    }
    line = trim(geo_ref_name);
    {   // inquire (file=line, exist=exists)
        std::ifstream probe(line.c_str());
        exists = probe.good();
    }
    if (!exists) {
        add_path(line);
        std::ifstream probe2(line.c_str());
        exists = probe2.good();
        if (!exists) {
            std::ifstream probe3((line + ".mop").c_str());
            if (probe3.good()) line = line + ".mop";
            exists = probe3.good();
        }
        if (!exists) {
            mopend("GEO_REF file '" + trim(line) + "' does not exist.");
            return;
        }
    }
    std::ifstream in99(line.c_str());
    if (!in99.good()) {
        mopend("Problem opening file " + trim(line));
        return;
    }
    for (int n = 1; n <= numat; ++n) nat[n] = labels[n];
    i = natoms;
    in99.clear();
    in99.seekg(0);
    for (i = 1; i <= 100000; ++i) {
        if (!std::getline(in99, num)) break;
    }
    in99.clear();
    in99.seekg(0);
    map_atoms_A.assign(numat + 1, 0);
    tmp_txt.assign(numat + 1, " ");
    tmp_atmass.assign(numat + 1, 0.0);
    atom_no.assign(numat + 1, 0);
    for (int n = 1; n <= numat; ++n) tmp_txt[n] = txtatm[n];
    for (int n = 1; n <= numat; ++n) tmp_atmass[n] = atmass[n];
    if (i > numat) {
        txtatm.assign(i + 1, " ");
        txtatm1.assign(i + 1, " ");
        labels.assign(i + 1, 0);
        atmass.assign(i + 1, 0.0);
        geoa.assign(4, std::vector<double>(i + 1, 0.0));
        c.assign(4, std::vector<double>(i + 1, 0.0));
        coord.assign(4, std::vector<double>(i + 1, 0.0));
        lopt.assign(4, std::vector<int>(i + 1, 0));
        na.assign(i + 1, 0);
        nb.assign(i + 1, 0);
        nc.assign(i + 1, 0);
    }
    line_1 = trim(line);
    upcase(line_1, (int)line_1.size());
    if (index1(line_1, ".ARC") != 0) {
        for (;;) {
            if (!std::getline(in99, line_1)) break;
            if (index1(line_1, "HEAT OF FORMATION") > 0) arc_hof_2 = reada(line_1, 20);
            if (index1(line_1, " FINAL GEOMETRY OBTAINED") != 0) break;
        }
        if (!in99.good()) {
            in99.clear();
            in99.seekg(0);
        }
    }
    i = 0;
    for (;;) {
        ++i;
        if (!std::getline(in99, refkey_ref[i])) goto label99;
        if (!refkey_ref[i].empty() && refkey_ref[i][0] == '*') --i;
        if (i == 1) {
            if (index1(refkey_ref[i], " +") != 0) --i;
        }
        if (i == 3) break;
    }
    refkey_ref[4] = "    NULL";
    refkey_ref[5] = "    NULL";
    refkey_ref[6] = "    NULL";
    goto label97;
label99:
    std::fprintf(stdout, " File' %s' is faulty\n", trim(line).c_str());
    return;
label97:

    i = index1(keywrd, " GEO_REF") + 11;
    for (;;) {
        if (i > (int)keywrd.size()) break;
        if (keywrd[i - 1] == '"' || keywrd[i - 1] == '\'') break;
        ++i;
    }
    if (i <= (int)keywrd.size() && keywrd[i] == ' ') {
        if (index1(keywrd, " 0SCF") == 0 || index1(keywrd, " HTML ") == 0)
            std::fprintf(stdout, "\n%10s\n", "By default, a restraining force of 3 kcal/mol/A^2 will be used");
        density = 3.0;
    } else {
        density = reada(keywrd, j + 2);
        std::fprintf(stdout, "\n%10s%8.3f%s\n", "A restraining force of", density,
                     " kcal/mol/A^2 will be used");
    }
    if (id == 3 && std::fabs(density) > 1.e-10) {
        line = "For solids, the restraining force in GEO_REF must be 0.0.";
        std::fprintf(stdout, "\n\n%10s\n", trim(line).c_str());
        mopend(trim(line));
        return;
    }
    numat_dat = numat;
    //
    //  Put atoms in the reference data set into the same order as those in the current data set
    //
    ii = numat;
    j = maxtxt;
    getgeo(99, labels, geoa, c, lopt, na, nb, nc, intern);
    //  (Fortran deallocate(c) is a no-op for the C++ global c.)
    maxtxt = j;
    if (moperr && natoms == 0) {
        mopend("An error was detected in the reference data set.");
        return;
    }
    if (natoms == -2) {
        j = ir;
        ir = 99;
        //  rewind(99): the file must be re-read from the start.
        in99.clear();
        in99.seekg(0);
        getpdb(geoa);
        ir = j;
        numat_old = numat;
        pdb = true;
    } else {
        pdb = false;
    }
    for (int n = 1; n <= ii; ++n) txtatm1[n] = tmp_txt[n];
    numat = numat - id;
    i = std::max(ii, numat);
    numat_ref = numat;
    tmp_geoa.assign(4, std::vector<double>(i + 1, 0.0));
    same.assign(i + 1, false);
    ok.assign(i + 1, false);
    diffs.assign(i + 1, " ");
    for (int n = 1; n <= i; ++n) diffs[n] = " ";
    k = 0;
    l = 0;
    if (index1(keywrd, " 0SCF") != 0 && index1(keywrd, " HTML ") != 0) {
        //
        // Put atoms in the GEO_REF file into the same order as those in the data-set
        //
        //    GEO_REF arrays:   geoa, txtatm, labels
        //    dataset arrays:   geo, txtatm1, nat
        //
        //  Align atoms in three passes...
        //
        for (j = 1; j <= ii; ++j) {
            j4 = 0;
            //  First pass
            for (i = 1; i <= numat; ++i) {
                if (same[i]) continue;
                if (substr1(txtatm1[j], 12) == substr1(txtatm[i], 12)) {
                    ++k;
                    tmp_txt[k] = txtatm[i];
                    coord[0][k] = geo[1][j];
                    coord[1][k] = geo[2][j];
                    coord[2][k] = geo[3][j];
                    tmp_geoa[1][k] = geoa[1][i];
                    tmp_geoa[2][k] = geoa[2][i];
                    tmp_geoa[3][k] = geoa[3][i];
                    labels[k] = nat[j];
                    same[i] = true;
                    ok[j] = true;
                    j4 = 1;
                    break;
                }
            }
        }
        let = (index1(keywrd, " LET") != 0);
        for (j = 1; j <= ii; ++j) {
            //  Second pass
            if (ok[j]) continue;
            for (i4 = 1; i4 <= numat; ++i4) {
                if (same[i4]) continue;
                if (substr1(txtatm1[j], 14, 14) == "H" && substr1(txtatm[i4], 14, 14) == "H" &&
                    substr1(txtatm1[j], 18) == substr1(txtatm[i4], 18)) {
                    ++k;
                    tmp_txt[k] = txtatm[i4];
                    coord[0][k] = geo[1][j];
                    coord[1][k] = geo[2][j];
                    coord[2][k] = geo[3][j];
                    tmp_geoa[1][k] = geoa[1][i4];
                    tmp_geoa[2][k] = geoa[2][i4];
                    tmp_geoa[3][k] = geoa[3][i4];
                    labels[k] = nat[j];
                    same[i4] = true;
                    ok[j] = true;
                    break;
                }
            }
            if (!ok[j] && !let) {
                ++l;
                diffs[l] = "   " + trim(txtatm1[j]);
                if (diffs[l].size() > 40) diffs[l] = diffs[l].substr(0, 40);
            }
        }
        if (let) {
            first = true;
            for (j = 1; j <= ii; ++j) {
                //  Third pass
                if (ok[j]) continue;
                for (i4 = 1; i4 <= numat; ++i4) {
                    if (same[i4]) continue;
                    if (substr1(txtatm1[j], 14, 14) == "H" && substr1(txtatm[i4], 14, 14) == "H") {
                        if (first) {
                            std::fprintf(stdout, "\n%16s\n", "Hydrogen atoms that are on different residues");
                            if (trim(job_fn) == trim(geo_dat_name)) {
                                line = "dataset";
                            } else {
                                line = "GEO_DAT";
                            }
                            std::fprintf(stdout, "\n%6s%14s\n", ("Atoms in " + trim(line) + " only").c_str(),
                                         "     Atoms in GEO_REF only");
                            first = false;
                        }
                        line = txtatm1[j];
                        if ((int)line.size() < 41) line.resize(41, ' ');
                        line.replace(40, std::string::npos, txtatm[i4]);
                        std::fprintf(stdout, "%4s\n", trim(line).c_str());
                        ++k;
                        tmp_txt[k] = txtatm[i4];
                        coord[0][k] = geo[1][j];
                        coord[1][k] = geo[2][j];
                        coord[2][k] = geo[3][j];
                        tmp_geoa[1][k] = geoa[1][i4];
                        tmp_geoa[2][k] = geoa[2][i4];
                        tmp_geoa[3][k] = geoa[3][i4];
                        labels[k] = nat[j];
                        same[i4] = true;
                        ok[j] = true;
                        break;
                    }
                }
                if (!ok[j]) {
                    ++l;
                    diffs[l] = "   " + trim(txtatm1[j]);
                    if (diffs[l].size() > 40) diffs[l] = diffs[l].substr(0, 40);
                }
            }
            if (!first) std::fprintf(stdout, "\n%18s\n", "(Atom labels from GEO_REF will be used)");
        }
        for (int n = 1; n <= k; ++n) {
            geo[1][n] = coord[0][n];
            geo[2][n] = coord[1][n];
            geo[3][n] = coord[2][n];
        }
        ii = k;
        jj = 0;
        for (i = 1; i <= numat; ++i) {
            if (!same[i]) {
                ++jj;
                diffs[jj] = (jj <= (int)diffs.size() - 1 ? diffs[jj] : std::string("")) +
                            "   " + trim(txtatm[i]);
            }
        }
        for (int n = 1; n <= k; ++n) nat[n] = labels[n];
        for (int n = 1; n <= k; ++n) txtatm[n] = tmp_txt[n];
        l = std::max(l, jj);
        if (l > 0) {
            std::fprintf(stdout, "\n%28s\n", "Differences in atoms sets");
            std::fprintf(stdout, "\n%10s%14s\n", "Atoms in data-set only", "    Atoms in GEO_REF only");
            for (i = 1; i <= l; ++i) std::fprintf(stdout, "%4s\n", trim(diffs[i]).c_str());
        }
        for (int n = 1; n <= k; ++n) {
            geoa[1][n] = tmp_geoa[1][n];
            geoa[2][n] = tmp_geoa[2][n];
            geoa[3][n] = tmp_geoa[3][n];
        }
        for (int n = 1; n <= k; ++n) tmp_txt[n] = txtatm[n];
        if (numat_dat != numat_ref) {
            if (index1(keywrd, " GEO-OK") == 0) {
                char buf[160];
                std::snprintf(buf, sizeof(buf), "Number of atoms in data-set:%5d, in GEO_REF:%5d", numat_dat, numat_ref);
                line = buf;
                mopend(trim(line));
                i = index1(keywrd, " GEO_DAT") + 9;
                if (i > 9) {
                    std::fprintf(stdout, "\n%10s\n%10s\n", "Docking can only be done when the number of atoms in ",
                                 "the files defined by GEO_DAT and GEO_REF are the same.");
                } else {
                    std::fprintf(stdout, "\n%5s\n",
                                 "Docking can only be done when the number of atoms in the data-set and GEO_REF are the same.");
                }
                i = (k * 100) / std::max(numat, ii);
                std::fprintf(stdout, "%10s%3d%s\n", "(The two systems have", i,
                             "% of atoms in common. To continue, but using");
                std::fprintf(stdout, "%10s\n", "only those atoms that are common to both systems, add \"GEO-OK\")");
                i = index1(keywrd, " GEO_DAT") + 9;
                if (i > 9) {
                    std::fprintf(stdout, "\n%10s\n", ("GEO_DAT name: \"" + trim(geo_dat_name) + "\"").c_str());
                } else {
                    std::fprintf(stdout, "\n%10s\n", ("Data-set name: " + trim(geo_dat_name)).c_str());
                }
                std::fprintf(stdout, "\n%10s\n", ("GEO_REF name: \"" + trim(geo_ref_name) + "\"").c_str());
                return;
            }
            std::fprintf(stdout, "\n%22s%5d\n", "Number of atoms in dataset:", numat_dat);
            std::fprintf(stdout, "%22s%5d\n", "Number of atoms in GEO_REF:", numat_ref);
            std::fprintf(stdout, "%10s%5d\n", "Number of atoms common to both systems:", ii);
        }
        natoms = k + id;
        numat = k;
        if (maxtxt == 26) {
            for (int n = 1; n <= numat; ++n) txtatm1[n] = txtatm[n];
        }
        for (int n = 1; n <= numat; ++n) txtatm[n] = tmp_txt[n];
        for (int n = 1; n <= numat; ++n) atmass[n] = tmp_atmass[n];
        for (int n = 1; n <= numat; ++n) {
            coord[0][n] = geoa[1][n];
            coord[1][n] = geoa[2][n];
            coord[2][n] = geoa[3][n];
        }
    }

    if (moperr) {
        i = index1(keywrd, " GEO_REF");
        j = index1(substr1(keywrd, i + 10), "\"") + i + 8;
        line = substr1(keywrd, i + 10, j);
        line = "Fault detected in GEO_REF data set: '" + trim(line) + "'";
        std::fprintf(stdout, "%s\n", trim(line).c_str());
        mopend(trim(line));
        return;
    }
    for (j = 2; j <= natoms; ++j)
        if (na[j] > 0) break;
    if (j <= natoms) {
        gmetry(geoa, coord);
        for (int n = 1; n <= natoms; ++n) {
            geoa[1][n] = coord[0][n];
            geoa[2][n] = coord[1][n];
            geoa[3][n] = coord[2][n];
        }
        for (int n = 1; n <= natoms; ++n) na[n] = 0;
    }
    if (maxtxt == 26 && substr1(txtatm[1], 1, 4) == "ATOM" && substr1(txtatm1[1], 1, 4) == "ATOM" &&
        index1(keywrd, "GEO-OK") == 0 && index1(keywrd, "\"SELF\"") == 0) {
        //
        //  Put atoms in the reference data set into the same order as those in the current data set
        //
        bug = false;
        any_bug = false;
        //
        //  Check for faults in the data-set labels
        //
        compare_txtatm(bug, any_bug);
        bug = false;
        //
        //  Check for faults in the geo_ref labels
        //
        for (i = 1; i <= numat; ++i) {
            if (labels[i] == 1) {
                l = 1;
                for (j = i + 1; j <= numat; ++j) {
                    if (labels[j] == 1) {
                        if (substr1(txtatm[i], 12) == substr1(txtatm[j], 12)) {
                            for (k = 13; k <= 16; ++k) {
                                if (substr1(txtatm[j], k, k) == " ") {
                                    ++l;
                                    txtatm[j][k - 1] = static_cast<char>(l + '0');
                                    txtatm1[j][k - 1] = static_cast<char>(l + '0');
                                    break;
                                }
                            }
                        }
                    }
                    if (l == 9) break;
                }
                if (l > 1) {
                    for (k = 13; k <= 16; ++k) {
                        if (substr1(txtatm[i], k, k) == " ") {
                            txtatm[i][k - 1] = '1';
                            txtatm1[i][k - 1] = '1';
                            break;
                        }
                    }
                }
            }
        }
        k = 0;
        for (i = 1; i <= numat; ++i) {
            for (j = i + 1; j <= numat; ++j) {
                if (substr1(txtatm1[i], 12) == substr1(txtatm1[j], 12)) break;
            }
            if (j <= numat) {
                if (!bug) {
                    ii = index1(keywrd, " GEO_REF");
                    jj = index1(substr1(keywrd, ii + 10), "\"") + ii + 8;
                    line = substr1(keywrd, ii + 10, jj);
                    std::fprintf(stdout, "\n%10s\n", ("Atoms in the GEO_REF file '" + trim(line) + "' with the same labels").c_str());
                }
                std::fprintf(stdout, "%10s%6d%6s%6d%6s\n", "Atoms", i, " and", j,
                             (";  Labels: (" + trim(txtatm1[i]) + ") and (" + trim(txtatm1[j]) + ")").c_str());
                bug = true;
                ++k;
            }
        }
        if (bug) {
            std::fprintf(stdout, "\n%10s\n", "Edit GEO_REF data set and re-run.");
            any_bug = true;
            bug = false;
        }
        //
        //  Check for labels in the data-set that are not in geo_ref
        //
        for (i = 1; i <= numat; ++i) {
            for (j = 1; j <= numat; ++j) {
                if (substr1(txtatm[i], 12) == substr1(txtatm1[j], 12)) break;
            }
            if (j > numat) {
                if (!bug) {
                    ii = index1(keywrd, "GEO_DAT=");
                    if (ii > 0) {
                        std::fprintf(stdout, "\n%10s\n", ("Atoms in the GEO_DAT file '" + trim(geo_dat_name) +
                                                          "' with no equivalent in the geo_ref file").c_str());
                    } else {
                        std::fprintf(stdout, "\n%10s\n", ("Atoms in the data-set file '" + trim(geo_dat_name) +
                                                          "' with no equivalent in the geo_ref file").c_str());
                    }
                }
                std::fprintf(stdout, "%10s%6d%6s\n", "Atom", i, (";  Label: (" + trim(txtatm[i]) + ")").c_str());
                bug = true;
                any_bug = true;
                ++k;
            }
        }
        bug = false;
        for (i = 1; i <= numat; ++i) {
            for (j = 1; j <= numat; ++j) {
                if (substr1(txtatm1[i], 12) == substr1(txtatm[j], 12)) break;
            }
            if (j > numat) {
                if (!bug) {
                    std::fprintf(stdout, "\n%10s\n", ("Atoms in the GEO_REF file '" + trim(geo_ref_name) +
                                                      "' with no equivalent in the data-set file").c_str());
                }
                std::fprintf(stdout, "%10s%6d%6s\n", "Atom", i, (";  Label: (" + trim(txtatm1[i]) + ")").c_str());
                bug = true;
                any_bug = true;
                ++k;
            }
        }
        if (any_bug) {
            if (k == 1) {
                mopend("Fault detected in atom labels in a GEO_REF run.");
            } else {
                mopend("Faults detected in atom labels in a GEO_REF run.");
            }
            mopend("(To continue with the current data set, use 'GEO-OK')");
            web_message(iw, "geo_ref.html");
            return;
        }
        map_atoms_A.assign(numat + 1, -1);
        for (i = 1; i <= numat; ++i) {
            for (j = 1; j <= numat; ++j) {
                if (substr1(txtatm[i], 12) == substr1(txtatm1[j], 12)) break;
            }
            if (j <= numat) map_atoms_A[j] = i;
        }
        for (i = 1; i <= numat; ++i) {
            if (map_atoms_A[i] == -1) {
                for (k = 1; k <= numat; ++k) {
                    for (j = 1; j <= numat; ++j) {
                        if (map_atoms_A[j] == k) break;
                    }
                    if (j > numat) {
                        map_atoms_A[i] = k;
                        break;
                    }
                }
            }
        }
        for (i = 1; i <= numat; ++i) {
            for (j = 1; j <= numat; ++j) {
                if (substr1(txtatm[i], 12) == substr1(txtatm1[j], 12)) break;
            }
            if (j <= numat) map_atoms_A[j] = i;
        }
        //
        //  At this point, the data set geometry is used in defining the atom order.
        //
        for (int n = 1; n <= numat; ++n) tmp_txt[n] = txtatm1[n];
        for (i = 1; i <= numat; ++i) {
            j = map_atoms_A[i];
            geoa[1][j] = coord[0][i];
            geoa[2][j] = coord[1][i];
            geoa[3][j] = coord[2][i];
            txtatm1[j] = tmp_txt[i];
            atom_no[j] = i;
        }
    } else {
        for (i = 1; i <= numat; ++i) atom_no[i] = i;
    }
    for (int n = 1; n <= numat; ++n) labels[n] = nat[n];
    for (int n = 1; n <= numat; ++n) {
        coord[0][n] = geoa[1][n];
        coord[1][n] = geoa[2][n];
        coord[2][n] = geoa[3][n];
    }
    rms_min = 1.e6;
    if ((int)same.size() <= numat) same.assign(numat + 1, false);
    geo_diff(sum, rms, false);
    toler = 3.0;
    i4 = 1;
    j4 = numat;
    k4 = 1;
    iquit = 0;
    xmin = 1.e8;
    exists = true;
    swap = (index1(keywrd, " NOSWAP") == 0);

    for (ii = 1; ii <= 30; ++ii) {
        dock(geoa, geo, sum);
        sum = 0.0;
        for (i = 1; i <= numat; ++i) {
            sum += std::sqrt((geo[1][i] - geoa[1][i]) * (geo[1][i] - geoa[1][i]) +
                             (geo[2][i] - geoa[2][i]) * (geo[2][i] - geoa[2][i]) +
                             (geo[3][i] - geoa[3][i]) * (geo[3][i] - geoa[3][i]));
        }
        if (swap && ii > 1) {
            //
            //  Check that atoms are in maximum coincidence.
            //  If they are not, then re-arrange geoa
            //
            for (i = 1; i <= numat; ++i) {
                sum = (geo[1][i] - geoa[1][i]) * (geo[1][i] - geoa[1][i]) +
                      (geo[2][i] - geoa[2][i]) * (geo[2][i] - geoa[2][i]) +
                      (geo[3][i] - geoa[3][i]) * (geo[3][i] - geoa[3][i]);
                same[i] = (sum < toler);
            }
            jj = 0;
            for (i = i4; i4 < j4 ? (i <= j4) : (i >= j4); i += k4) {
                sum = (geo[1][i] - geoa[1][i]) * (geo[1][i] - geoa[1][i]) +
                      (geo[2][i] - geoa[2][i]) * (geo[2][i] - geoa[2][i]) +
                      (geo[3][i] - geoa[3][i]) * (geo[3][i] - geoa[3][i]);
                if (sum > toler) {
                    sum = 1.e8;
                    k = i;
                    for (j = i4; i4 < j4 ? (j <= j4) : (j >= j4); j += k4) {
                        if (substr1(txtatm1[i], 1, 6) == substr1(txtatm1[j], 1, 6) &&
                            substr1(txtatm1[i], 22) == substr1(txtatm1[j], 22) &&
                            substr1(txtatm1[i], 14, 14) == substr1(txtatm1[j], 14, 14) &&
                            labels[i] == labels[j] && !(same[i] && same[j])) {
                            num = substr1(txtatm1[i], 15, 15);
                            if (num == " " || (num >= "0" && num <= "9")) {
                                num = substr1(txtatm1[j], 15, 15);
                                if (num == " " || (num >= "0" && num <= "9")) num = "+";
                            }
                            if (num != "+") continue;
                            sum2 = (geo[1][j] - geoa[1][j]) * (geo[1][j] - geoa[1][j]) +
                                   (geo[2][j] - geoa[2][j]) * (geo[2][j] - geoa[2][j]) +
                                   (geo[3][j] - geoa[3][j]) * (geo[3][j] - geoa[3][j]);
                            sum3 = (geo[1][i] - geoa[1][i]) * (geo[1][i] - geoa[1][i]) +
                                   (geo[2][i] - geoa[2][i]) * (geo[2][i] - geoa[2][i]) +
                                   (geo[3][i] - geoa[3][i]) * (geo[3][i] - geoa[3][i]);
                            dum1 = (geo[1][i] - geoa[1][j]) * (geo[1][i] - geoa[1][j]) +
                                   (geo[2][i] - geoa[2][j]) * (geo[2][i] - geoa[2][j]) +
                                   (geo[3][i] - geoa[3][j]) * (geo[3][i] - geoa[3][j]);
                            dum2 = (geoa[1][i] - geo[1][j]) * (geoa[1][i] - geo[1][j]) +
                                   (geoa[2][i] - geo[2][j]) * (geoa[2][i] - geo[2][j]) +
                                   (geoa[3][i] - geo[3][j]) * (geoa[3][i] - geo[3][j]);
                            //
                            //  If swapping the atoms around will improve the overlap, then do so.
                            //
                            if (dum1 + dum2 < sum - 0.1 && dum1 + dum2 < sum2 + sum3 - 0.1) {
                                sum = dum1 + dum2;
                                sum4 = sum2 + sum3;
                                k = j;
                            }
                        }
                    }
                    if (i != k) {
                        ++jj;
                        if (exists) {
                            exists = false;
                            std::fprintf(stdout, "%s\n", std::string(80, '-').c_str());
                            if (maxtxt == 26) {
                                std::fprintf(stdout, "\n%4s\n\n%8s\n", "  List of atoms in GEO_REF that are swapped in order to maximize overlap",
                                             "      Atom Label          and         Atom Label         Difference to GEO_REF");
                            } else {
                                std::fprintf(stdout, "\n%4s\n\n%8s\n", "  List of atoms that are swapped in order to maximize overlap",
                                             "  Atom name   Atom No. and Atom No.  Difference to GEO_REF");
                            }
                        }
                        if (maxtxt == 26) {
                            std::fprintf(stdout, "%5d %2s %2s %12.2f\n", jj,
                                         ("(" + trim(txtatm1[i]) + ")").c_str(),
                                         ("(" + trim(txtatm1[k]) + ")").c_str(),
                                         std::sqrt(sum) - std::sqrt(sum4));
                        } else {
                            std::fprintf(stdout, "%5d %1s %9d %13d %1s %12.2f\n", jj,
                                         atom_names[nat[i]].c_str(), atom_no[i], atom_no[k], " ",
                                         std::sqrt(sum) - std::sqrt(sum4));
                        }
                        for (j = 1; j <= 3; ++j) {
                            sum = geoa[j][i];
                            geoa[j][i] = geoa[j][k];
                            geoa[j][k] = sum;
                        }
                        line_1 = txtatm1[i];
                        txtatm1[i] = txtatm1[k];
                        txtatm1[k] = trim(line_1);
                        j = atom_no[i];
                        atom_no[i] = atom_no[k];
                        atom_no[k] = j;
                    }
                    same[i] = true;
                }
            }
        }
        if (!exists) std::fprintf(stdout, "\n");
        k = 0;
        for (i = 1; i <= numat; ++i) {
            for (l = 1; l <= 3; ++l) {
                ++k;
                xparam[k] = geo[l][i];
            }
        }
        geo_diff(sum, rms, false);
        sum1 = sum / numat;
        sum2 = std::sqrt(rms / numat);
        if (std::fabs(rms_min - rms) < 1.e-4 && ii > 2) break;
        if (xmin > rms) {
            xmin = rms;
            iquit = 0;
        } else {
            ++iquit;
            if (iquit > 3) break;
        }
        rms_min = rms;
        toler = std::max(1.0, 0.7 * toler);
        k4 = i4;
        i4 = j4;
        j4 = k4;
        if (i4 == 1) {
            k4 = 1;
        } else {
            k4 = -1;
        }
    }
    for (int n = 1; n <= numat; ++n) txtatm1[n] = txtatm[n];
    if (sum1 < 1.e-10 && index1(keywrd, " LOCATE-TS") != 0) {
        line = " The data-set and GEO_REF geometries are the same.  Correct the fault and re-run.";
        mopend(trim(line));
        return;
    }
    if (index1(keywrd, "\"SELF\"") == 0) std::fprintf(stdout, "\n%25s\n", "After docking");
    geo_diff(sum, rms, true);
    if (trim(job_fn) == trim(geo_dat_name)) {
        line = "dataset";
    } else {
        line = "GEO_DAT";
    }
    if (index1(keywrd, "\"SELF\"") == 0)
        std::fprintf(stdout, "\n%3s%8.2f%8s%8.4f%8s%8.4f%8s\n",
                     ("Difference between " + trim(line) + " and GEO_REF:").c_str(),
                     sum1 * numat, " = total,", sum1, " = Average,", sum2, " = RMS, in Angstroms");
    in99.close();
    if (index1(keywrd, " 0SCF") != 0 && index1(keywrd, " HTML ") != 0) {
        //
        // Write out the reference geometry file needed for comparing the two geometries using JSmol.
        //
        line = geo_ref_name.substr(0, (int)geo_ref_name.size() - 3) + "pdb";
        add_path(line);
        std::ofstream out99(line.c_str());
        l_control("HTML", (int)std::string("HTML").size(), -1);
        for (int n = 1; n <= numat; ++n) {
            coord[0][n] = geoa[1][n];
            coord[1][n] = geoa[2][n];
            coord[2][n] = geoa[3][n];
        }
        pdbout(99);
        l_control("HTML", (int)std::string("HTML").size(), 1);
        out99.close();
    }
    for (int n = 1; n <= natoms; ++n) {
        coord[0][n] = geo[1][n];
        coord[1][n] = geo[2][n];
        coord[2][n] = geo[3][n];
    }   // Store input geometry
    for (i = 1; i <= nvar; ++i) xparam[i] = geoa[loc[2][i]][loc[1][i]];
    for (int n = 1; n <= numat; ++n) {
        geo[1][n] = geoa[1][n];
        geo[2][n] = geoa[2][n];
        geo[3][n] = geoa[3][n];
    }
    //
    // The following "big_swap" saves the geometry and xparam for the rotated input geometry
    //
    big_swap(0, 2);
    for (int n = 1; n <= natoms; ++n) {
        geo[1][n] = coord[0][n];
        geo[2][n] = coord[1][n];
        geo[3][n] = coord[2][n];
    }   // Restore input geometry
    for (i = 1; i <= nvar; ++i) xparam[i] = geo[loc[2][i]][loc[1][i]];

    if (index1(keywrd, " 0SCF") != 0) {
            wrt_diffs();
            if (std::fabs(arc_hof_1) + std::fabs(arc_hof_2) > 1.e-4) std::fprintf(stdout, "\n");
        if (trim(job_fn) == trim(geo_dat_name)) {
            line = "dataset";
        } else {
            line = "GEO_DAT";
        }
        if (std::fabs(arc_hof_1) > 1.e-4)
            std::fprintf(stdout, "%10s%12.3f%s\n", ("Heat of formation of " + trim(line) + " system:").c_str(),
                         arc_hof_1, " Kcal/mol");
        if (std::fabs(arc_hof_2) > 1.e-4)
            std::fprintf(stdout, "%10s%12.3f%s\n", "Heat of formation of GEO_REF system:", arc_hof_2, " Kcal/mol");
        if (std::fabs(arc_hof_1) > 1.e-4 && std::fabs(arc_hof_2) > 1.e-4)
            std::fprintf(stdout, "%40s%12.3f%s\n", "Diff.:", arc_hof_1 - arc_hof_2, " Kcal/mol");
        analyze_h_bonds();
    }
    for (ii = 1; ii <= 6; ++ii) {
        line = " " + rtrim(refkey[ii]);
        upcase(line, (int)line.size());
        i = index1(line, " GEO_REF");
        if (i != 0) break;
    }
    if (index1(keywrd, "\"SELF\"") == 0) {
        if (index1(keywrd, " 0SCF") == 0 || index1(keywrd, " HTML ") == 0) {
            line = geo_dat_name.substr(0, (int)geo_dat_name.size() - 3) + "new";
            add_path(line);
            std::ofstream out99(line.c_str());
            if (!out99.good()) {
                //  open failed: ignore (iostat only printed in real runs)
            }
        }
    }
    if (index1(keywrd, " TS ") != 0) {
        std::fprintf(stdout, "\n%s\n", "    The average of the supplied and reference geometry will be written to:");
        std::fprintf(stdout, "%s\n", ("'" + trim(line) + "'").c_str());
        for (int n = 1; n <= natoms; ++n) {
            geo[1][n] = 0.5 * (geoa[1][n] + geo[1][n]);
            geo[2][n] = 0.5 * (geoa[2][n] + geo[2][n]);
            geo[3][n] = 0.5 * (geoa[3][n] + geo[3][n]);
        }
        geout(99);
    } else if (index1(keywrd, "\"SELF\"") == 0) {
        if (index1(keywrd, " 0SCF") == 0 || index1(keywrd, " HTML ") == 0) {
            if (index1(keywrd, " SWAP ") != 0) {
                std::fprintf(stdout, "\n%4s\n", ("The rotated-translated-re-arranged input geometry will be written to: '" + trim(line) + "'").c_str());
            } else {
                std::fprintf(stdout, "\n%s\n", ("    The rotated-translated input geometry will be written to: '" + trim(line) + "'").c_str());
            }
        } else {
            for (i = (int)output_fn.size(); i >= 1; --i) {
                if (i > (int)output_fn.size()) continue;
                if (output_fn[i - 1] == '/' || output_fn[i - 1] == '\\') break;
            }
            if (i > 0) {
                std::fprintf(stdout, "\n\n%s\n", ("    The following files will be written to \"" + output_fn.substr(0, i) + "\":").c_str());
            } else {
                std::fprintf(stdout, "\n\n%s\n", "    The following files will be written:");
            }
            line = job_fn.substr(i, (int)job_fn.size() - 3 - i) + "html";
            std::fprintf(stdout, "%s\n", ("    '" + trim(line) + "'").c_str());
            line = geo_dat_name.substr(0, (int)geo_dat_name.size() - 3) + "pdb";
            std::fprintf(stdout, "%s\n", ("    '" + trim(line) + "'").c_str());
            line = geo_ref_name.substr(0, (int)geo_ref_name.size() - 3) + "pdb";
            std::fprintf(stdout, "%s\n", ("    '" + trim(line) + "'").c_str());
        }
        line = refkey_ref[1];
        for (i = 1; i <= 3; ++i) {
            line = refkey[i];
            refkey[i] = refkey_ref[i];
            refkey_ref[i] = line;
        }
        refkey_ref[4] = koment;
        refkey_ref[5] = title;
        if (index1(keywrd, " 0SCF") != 0 && index1(keywrd, " HTML ") != 0) {
            //
            // Write out the input geometry file needed for comparing the two geometries using JSmol.
            //
            if (sum2 > 9.9995) {
                num = "6";
            } else {
                num = "5";
            }
            char buf2[80];
            std::snprintf(buf2, sizeof(buf2), "RMS_DIFF=%f", sum2);
            line = buf2;
            l_control(trim(line), (int)trim(line).size(), 1);
            line = geo_dat_name.substr(0, (int)geo_dat_name.size() - 3) + "pdb";
            for (int n = 1; n <= numat; ++n) {
                coord[0][n] = geo[1][n];
                coord[1][n] = geo[2][n];
                coord[2][n] = geo[3][n];
            }
            add_path(line);
            std::ofstream out99b(line.c_str());
            pdbout(99);
        } else {
            koment = " NULL";
            title = " NULL";
            if (pdb) {
                for (int n = 1; n <= natoms; ++n) {
                    coord[0][n] = geoa[1][n];
                    coord[1][n] = geoa[2][n];
                    coord[2][n] = geoa[3][n];
                }
                html = (index1(keywrd, " HTML") != 0);
                if (html) l_control("HTML", (int)std::string("HTML").size(), -1);
                pdbout(99);
                if (html) l_control("HTML", (int)std::string("HTML").size(), 1);
                for (int n = 1; n <= numat; ++n) {
                    coord[0][n] = geo[1][n];
                    coord[1][n] = geo[2][n];
                    coord[2][n] = geo[3][n];
                }   // Then restore the input geometry
            } else {
                for (int n = 1; n <= natoms; ++n) {
                    geo[1][n] = geoa[1][n];
                    geo[2][n] = geoa[2][n];
                    geo[3][n] = geoa[3][n];
                }   // Temporarily, put the reference geometry into geo
                        geout(99);   // Write out the re-organized reference data set.
                        for (int n = 1; n <= numat; ++n) {
                    geo[1][n] = coord[0][n];
                    geo[2][n] = coord[1][n];
                    geo[3][n] = coord[2][n];
                }   // Then restore the input geometry
            }
        }
        for (i = 1; i <= 3; ++i) {
            line = refkey_ref[i];
            refkey_ref[i] = refkey[i];
            refkey[i] = line;
        }
        koment = trim(refkey_ref[4]);
        title = trim(refkey_ref[5]);
    }
    if (index1(keywrd, " 0SCF") + index1(keywrd, " TS ") != 0) {
        mopend("GEO_REF with 0SCF: Job complete");
        return;
    }
    for (int n = 1; n <= natoms; ++n) {
        geo[1][n] = coord[0][n];
        geo[2][n] = coord[1][n];
        geo[3][n] = coord[2][n];
    }
    return;
}

void compare_txtatm(bool& bug, bool& any_bug) {
    if (index1(keywrd, "GEO-OK") != 0) return;
    if (index1(keywrd, "GEO_REF") == 0) return;
    bug = false;
    if (maxtxt != 26) return;
    int i, j, ii, jj;
    for (i = 1; i <= numat; ++i) {
        for (j = 1; j <= 2; ++j) {
            if (substr1(txtatm[i], 20, 20) != " ") break;
            txtatm[i].replace(17, 3, " " + substr1(txtatm[i], 18, 19));
        }
        for (j = 1; j <= 2; ++j) {
            if (substr1(txtatm1[i], 20, 20) != " ") break;
            txtatm1[i].replace(17, 3, " " + substr1(txtatm1[i], 18, 19));
        }
    }
    for (i = 1; i <= numat; ++i) {
        for (j = i + 1; j <= numat; ++j) {
            if (substr1(txtatm[i], 12) == substr1(txtatm[j], 12)) break;
        }
        if (j <= numat && nat[i] != 1) {
            if (!bug) {
                ii = index1(keywrd, "GEO_DAT=");
                if (ii > 0) {
                    jj = index1(substr1(keywrd, ii + 9), "\"") + ii + 7;
                    line = substr1(keywrd, ii + 9, jj);
                    std::fprintf(stdout, "\n%10s\n", ("Atoms in the GEO_DAT file '" + trim(line) + "' with the same labels").c_str());
                } else {
                    std::fprintf(stdout, "\n%10s\n", ("Atoms in the data-set file '" + trim(job_fn) + "' with the same labels").c_str());
                }
                std::fprintf(stdout, "%10s%6d%6s%6d%6s\n", "Atoms", i, " and", j,
                             (";  Labels: (" + trim(txtatm[i]) + ") and (" + trim(txtatm[j]) + ")").c_str());
            }
            bug = true;
        }
    }
    if (bug) {
        mopend("Error in data detected while using GEO_REF");
        std::fprintf(stdout, "%5s\n", "(To continue with the current data set, use 'GEO-OK')");
        any_bug = true;
    }
}
