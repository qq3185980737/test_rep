// lewis.cpp — C++ translation of lewis.F90 (MOPAC 2016).
//
//  LEWIS works out what atom is bonded to what atom for the purpose of
//  generating a Lewis structure.  It starts with the topography generated
//  by "set_up_dentate" (the raw topography), then tidies it up using a
//  series of filters - keywords such as CVB, and criteria such as a
//  hydrogen atom must bond to exactly one other atom, and that must be a
//  non-hydrogen atom.
//
//  On exit, NBONDS holds the number of bonds to each atom, and
//           IBONDS contains the atom numbers, in order of increasing
//                  atomic number.
#include "lewis.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "mod_atomradii.h"
#include "reada.h"

#include "mopend.h"
using namespace molkst_C;
using namespace common_arrays_C;


double distance(int i, int j);   // interatomic distance (provided externally)
void set_up_dentate();           // raw topography (provided externally)

namespace {

// 1-based find: position of "sub" in "s" starting at 1-based "pos"; 0 if absent.
int index1(const std::string& s, const std::string& sub, int pos = 1) {
    if (pos < 1 || pos > (int)s.size() + 1) return 0;
    size_t p = s.find(sub, (size_t)(pos - 1));
    return p == std::string::npos ? 0 : (int)p + 1;
}
std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(' ');
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(' ');
    return s.substr(a, b - a + 1);
}
void append_log(const std::string& txt) {
    if (!chanel_C::log) return;
    FILE* f = std::fopen(chanel_C::log_fn.c_str(), "a");
    if (f) { std::fprintf(f, "%s\n", txt.c_str()); std::fclose(f); }
}

std::vector<std::vector<double>> store_coord;   // saved copy of coord (lewis)
int icalcn_lw = -1;                            // lewis: numcal memory
int icalcn_cv = -1;                            // check_CVS: numcal memory
int icalcn_h = -1;                             // check_h: numcal memory

}  // namespace

void remove_bond(int i) {
    //  Remove the longest bond from atom i
    double rmax = 0.0;
    int j, k, l = 0;
    for (j = 1; j <= nbonds[i]; ++j) {
        k = ibonds[j][i];
        double sum = (coord[0][i] - coord[0][k]) * (coord[0][i] - coord[0][k]) +
                     (coord[1][i] - coord[1][k]) * (coord[1][i] - coord[1][k]) +
                     (coord[2][i] - coord[2][k]) * (coord[2][i] - coord[2][k]);
        if (sum > rmax) { rmax = sum; l = k; }
    }
    k = 0;
    for (j = 1; j <= nbonds[i]; ++j)
        if (ibonds[j][i] != l) { ++k; ibonds[k][i] = ibonds[j][i]; }
    nbonds[i] = nbonds[i] - 1;
    k = 0;
    for (j = 1; j <= nbonds[l]; ++j)
        if (ibonds[j][l] != i) { ++k; ibonds[k][l] = ibonds[j][l]; }
    nbonds[l] = nbonds[l] - 1;
}

void check_h(int& ibad) {
    if (!mozyme) return;
    bool prt = (icalcn_h != numcal);
    if (prt) icalcn_h = numcal;
    std::vector<int> im(natoms + 1);
    int j = 0;
    for (int i = 1; i <= natoms; ++i)
        if (labels[i] != 99) { ++j; im[j] = i; }
    ibad = 0;
    for (int i = 1; i <= numat; ++i) {
        if (nat[i] != 1) continue;
        if (nbonds[i] == 0) {
            //  Find nearest heavy atom to hydrogen atom i
            double rmin = 1.e12;
            int k = 0;
            for (j = 1; j <= numat; ++j) {
                if (j == i || nat[j] == 1) continue;
                double r = distance(i, j);
                if (r < rmin) { rmin = r; k = j; }
            }
            nbonds[i] = 1;
            ibonds[1][i] = k;
            nbonds[k] = nbonds[k] + 1;
            ibonds[nbonds[k]][k] = i;
            if (rmin < std::max(1.3, 2.0 * mod_atomradii::radius[k])) continue;
            if (ibad < 20 &&
                index1(keywrd, " LET") + index1(keywrd, " ADD-H") + index1(keywrd, " RESEQ") == 0) {
                std::string el = trim(elemts_C::elemnt[nat[k]]);
                char buf[256];
                std::snprintf(buf, sizeof buf,
                              "Hydrogen atom %d is%6.1f Angstroms from the nearest non-hydrogen atom (%s %d)",
                              i, rmin, el.c_str(), im[k]);
                line = buf;
                mopend(trim(line));
                std::fprintf(stdout, "\n          (Atom label: '%s')\n", trim(txtatm[i]).c_str());
                std::fprintf(stdout, "(To continue, add \"LET\" to the keyword line)\n");
                std::fprintf(stdout, "(Label for atom %d: \"%s\")\n", i, txtatm[i].c_str());
                append_log(trim(line));
                ++ibad;
            }
            continue;
        } else if (nbonds[i] == 1) {
            //  The normal case
            double r = distance(i, ibonds[1][i]);
        } else {
            //  Find the heavy atom, k, nearest to hydrogen atom i
            double rmin = 1.e12;
            int m = 0;
            int k = 0;
            int heavy[10] = {0};
            for (j = 1; j <= nbonds[i]; ++j) {
                if (nat[ibonds[j][i]] == 1) continue;
                ++m;
                heavy[m] = ibonds[j][i];
                double r = distance(i, ibonds[j][i]);
                if (r < rmin && r > 0.95) { rmin = r; k = ibonds[j][i]; }
            }
            if (k == 0) {
                for (j = 1; j <= nbonds[i]; ++j) {
                    if (nat[ibonds[j][i]] == 1) continue;
                    ++m;
                    heavy[m] = ibonds[j][i];
                    double r = distance(i, ibonds[j][i]);
                    if (r < rmin) { rmin = r; k = ibonds[j][i]; }
                }
            }
            // Write out error message
            if (m == 2 && prt) {
                int n = heavy[2];
                m = heavy[1];
                std::string s1 = trim(elemts_C::elemnt[nat[m]]);
                std::string s2 = trim(elemts_C::elemnt[nat[n]]);
                char buf[256];
                std::snprintf(buf, sizeof buf,
                              " Hydrogen atom %d has covalent bonds to %s %d and %s %d",
                              im[i], s1.c_str(), im[m], s2.c_str(), im[n]);
                std::fprintf(stdout, "%s\n", buf);
                append_log(buf);
            } else if (m > 2) {
                if (ibad < 20) {
                    char buf[256];
                    std::snprintf(buf, sizeof buf,
                                  " Hydrogen atom %d has more than two covalent bonds",
                                  im[i]);
                    std::fprintf(stdout, "%s\n", buf);
                    append_log(buf);
                }
            }
            if (m > 1) {
                //  Remove all bonds to hydrogen atom i except the shortest bond, to atom "k"
                for (j = 1; j <= nbonds[i]; ++j) {
                    int ii = ibonds[j][i];
                    if (ii == k || nat[ii] == 1) continue;
                    m = 0;
                    for (int l = 1; l <= nbonds[ii]; ++l)
                        if (ibonds[l][ii] != i) { ++m; ibonds[m][ii] = ibonds[l][ii]; }
                    nbonds[ii] = m;
                }
                nbonds[i] = 1;
                ibonds[1][i] = k;
                if (prt)
                    std::fprintf(stdout, "                         (Keeping bond to %s %d)\n",
                                 trim(elemts_C::elemnt[nat[k]]).c_str(), im[k]);
            }
        }
        if (ibad == 21) std::fprintf(stdout, " Remaining errors suppressed\n");
    }
}


void check_CVS(bool let) {
    bool error = false;
    bool PDB_input = false;
    int j = 0;
    if (!store_coord.empty() && icalcn_cv != numcal) { store_coord.clear(); icalcn_cv = numcal; }
    if (store_coord.empty()) {
        int ncol = (int)coord[0].size();
        store_coord.assign(3, std::vector<double>(ncol, 0.0));
        for (int i = 0; i <= 2; ++i)
            for (int jj = 1; jj < ncol; ++jj) store_coord[i][jj] = coord[i][jj];
    }
    //  Parse CVB(2:3,45:65)
    int i = index1(keywrd, " CVB");
    if (i > 0) {
        int j = index1(keywrd.substr(i - 1), ") ") + i;
        line = keywrd.substr(i, (j - 1) - (i + 1) + 1);
        int ii = (int)trim(line).size();
        int k = 0;
        for (;;) {
            ++k;
            if (k >= ii) break;
            if (line[k - 1] == '"') {
                for (;;) { ++k; if (line[k - 1] == '"') break; }
            }
            if (line[k - 1] >= '0' && line[k - 1] <= '9') break;
        }
        PDB_input = (k > ii - 3);
    }
    i = index1(keywrd, " CVB");
    int k = 0;
    if (i != 0) k = index1(keywrd.substr(i - 1), ")") + i;
    if (k == i && i != 0) {
        k = index1(keywrd.substr(i + 2), " ") + i + 3;
        std::fprintf(stdout, "\n\n          No closing parenthesis for CVB keyword\n");
        std::fprintf(stdout, " CVB keyword =\"%s\"\n", keywrd.substr(i - 1, k - i + 1).c_str());
        mopend("ERROR IN CVB KEYWORD");
        return;
    }
    if (i != 0) {
        int m, l, jj;
        for (;;) {
            k = index1(keywrd.substr(i - 1), ") ") + i;
            jj = index1(keywrd.substr(i - 1, k - i + 1), "\"") + i;
            if (jj == i) break;
            txt_to_atom_no(keywrd, jj - 1, let, m);
            if (moperr) return;
            if (m > numat) {
                if (let) {
                    //  Jump over the faulty atom label
                    k = index1(keywrd.substr(i - 1), ")") + i;
                    i = index1(keywrd.substr(i - 1, k - i + 1), "\"") + i;
                    i = index1(keywrd.substr(i - 1, k - i + 1), "\"") + i;
                    continue;
                }
            }
        }
        k = index1(keywrd.substr(i - 1), ")") + i;
        int ncvb = 0;
        int cvb_atoms[60] = {0};
        for (;;) {
            jj = (int)std::lround((double)reada(keywrd.substr(i - 1, k - i + 1), 1));
            int ii2 = index1(keywrd.substr(i - 1, k - i + 1), ":") + i;
            if (ii2 == i) {
                std::fprintf(stdout, "\n ERROR IN CVB KEYWORD - ':' EXPECTED BUT NOT FOUND\n");
                jj = index1(keywrd, " CVB");
                std::fprintf(stdout, " CVB keyword =\"%s\"\n", keywrd.substr(jj - 1, k - jj + 1).c_str());
                mopend("ERROR IN CVB KEYWORD");
                return;
            }
            i = ii2;
            if (jj > numat || jj == 0 || jj < -numat) {
                std::fprintf(stdout, " AT LEAST ONE ATOM DEFINED BY CVB IS FAULTY\n");
                std::fprintf(stdout, " THE FAULTY ATOM NUMBER IS %d\n", jj);
                jj = index1(keywrd, " CVB");
                std::fprintf(stdout, " CVB keyword =\"%s\"\n", keywrd.substr(jj - 1, k - jj + 1).c_str());
                error = true;
            }
            l = (int)std::lround((double)reada(keywrd.substr(i - 1, k - i + 1), 1));
            if (l > numat || l == 0 || l < -numat) {
                std::fprintf(stdout, " AT LEAST ONE ATOM DEFINED BY CVB IS FAULTY\n");
                std::fprintf(stdout, " THE FAULTY ATOM NUMBER IS %d\n", l);
                int mm = index1(keywrd, " CVB");
                std::fprintf(stdout, " CVB keyword =\"%s\"\n", keywrd.substr(mm - 1, k - mm + 1).c_str());
                error = true;
            }
            if (error) goto label99;
            int n = abs(l);
            int m2 = abs(jj);
            if (!store_coord.empty() && !PDB_input) {
                //   If an atom has moved, find its new location
                int ii;
                for (ii = 1; ii <= numat; ++ii)
                    if (std::fabs(coord[0][ii] - store_coord[0][n]) < 1.e-1 &&
                        std::fabs(coord[1][ii] - store_coord[1][n]) < 1.e-1 &&
                        std::fabs(coord[2][ii] - store_coord[2][n]) < 1.e-1) break;
                l = (l >= 0) ? ii : -ii;
                for (ii = 1; ii <= numat; ++ii)
                    if (std::fabs(coord[0][ii] - store_coord[0][m2]) < 1.e-1 &&
                        std::fabs(coord[1][ii] - store_coord[1][m2]) < 1.e-1 &&
                        std::fabs(coord[2][ii] - store_coord[2][m2]) < 1.e-1) break;
                jj = (jj >= 0) ? ii : -ii;
            } else {
                jj = (jj >= 0) ? m2 : -m2;
                l = (l >= 0) ? n : -n;
            }
            int ii;
            for (ii = 1; ii <= ncvb; ++ii)
                if (cvb_atoms[ii] == m2) break;
            if (ii > ncvb) { ++ncvb; cvb_atoms[ncvb] = m2; }
            m2 = abs(l);
            for (ii = 1; ii <= ncvb; ++ii)
                if (cvb_atoms[ii] == m2) break;
            if (ii > ncvb) { ++ncvb; cvb_atoms[ncvb] = m2; }
            if (jj > 0 && l > 0) {
                for (ii = 1; ii <= nbonds[jj]; ++ii)
                    if (ibonds[ii][jj] == l && !let &&
                        index1(keywrd, " GEO-OK") == 0 && index1(keywrd, " 0SCF") == 0) {
                        std::fprintf(stdout,
                                     " The bond defined by CVB between atoms%5d and%5d already exists.  Correct error and re-submit\n",
                                     jj, l);
                        std::fprintf(stdout, "\n                              PDB label              Coordinates of the two atoms\n");
                        std::fprintf(stdout, " Atom:%5d: %s (%s) %12.6f %12.6f %12.6f\n",
                                     jj, trim(elemts_C::elemnt[nat[jj]]).c_str(), txtatm[jj].c_str(),
                                     coord[0][jj], coord[1][jj], coord[2][jj]);
                        std::fprintf(stdout, " Atom:%5d: %s (%s) %12.6f %12.6f %12.6f\n",
                                     l, trim(elemts_C::elemnt[nat[l]]).c_str(), txtatm[l].c_str(),
                                     coord[0][l], coord[1][l], coord[2][l]);
                double r = distance(jj, l);
                if (r > 3.0 * (mod_atomradii::radius[nat[jj]] + mod_atomradii::radius[nat[l]]) &&
                    !let && index1(keywrd, " GEO-OK") == 0 && index1(keywrd, " 0SCF") == 0) {
                    std::fprintf(stdout,
                                 " The bond defined by CVB between atoms%5d and%5d would be%5.1f Angstroms long.  This is unrealistic.\n",
                                 jj, l, r);
                    std::fprintf(stdout, "\n                              PDB label              Coordinates of the two atoms\n");
                    std::fprintf(stdout, " Atom:%5d: %s (%s) %12.6f %12.6f %12.6f\n",
                                 jj, trim(elemts_C::elemnt[nat[jj]]).c_str(), txtatm[jj].c_str(),
                                 coord[0][jj], coord[1][jj], coord[2][jj]);
                    std::fprintf(stdout, " Atom:%5d: %s (%s) %12.6f %12.6f %12.6f\n",
                                 l, trim(elemts_C::elemnt[nat[l]]).c_str(), txtatm[l].c_str(),
                                 coord[0][l], coord[1][l], coord[2][l]);
                    error = true;
                }
                nbonds[jj] = nbonds[jj] + 1;
                nbonds[l] = nbonds[l] + 1;
                ibonds[nbonds[jj]][jj] = l;
                ibonds[nbonds[l]][l] = jj;
            } else {
                //  make sure that the faulty bond actually exists.
                jj = abs(jj);
                l = abs(l);
                if (jj == 0 || l == 0) {
                    std::fprintf(stdout, "          Error detected while reading CVB\n");
                    jj = index1(keywrd, " CVB");
                    std::fprintf(stdout, " CVB keyword =\"%s\"\n", keywrd.substr(jj - 1, k - jj + 1).c_str());
                    error = true;
                }
                bool lbond = false;
                for (n = 1; n <= nbonds[jj]; ++n)
                    if (ibonds[n][jj] == l) {
                        //  Faulty bond exists, therefore delete it.
                        int ik = l;
                        for (m = n + 1; m <= nbonds[jj]; ++m) ibonds[m - 1][jj] = ibonds[m][jj];
                        ibonds[nbonds[jj]][jj] = ik;
                        nbonds[jj] = nbonds[jj] - 1;
                        lbond = true;
                        break;
                    }
                for (n = 1; n <= nbonds[l]; ++n)
                    if (ibonds[n][l] == jj) {
                        int ik = jj;
                        for (m = n + 1; m <= nbonds[l]; ++m) ibonds[m - 1][l] = ibonds[m][l];
                        ibonds[nbonds[l]][l] = ik;
                        nbonds[l] = nbonds[l] - 1;
                        lbond = true;
                        break;
                    }
                if (!lbond && !let) {
                    if (index1(keywrd, " GEO-OK") == 0 && index1(keywrd, " 0SCF") == 0) {
                        std::fprintf(stdout,
                                     " The bond defined by CVB between atoms%5d and%5d did not exist.  Correct error and re-submit\n",
                                     jj, l);
                        std::fprintf(stdout, "\n                              PDB label              Coordinates of the two atoms\n");
                        std::fprintf(stdout, " Atom:%5d: %s (%s) %12.6f %12.6f %12.6f\n",
                                     jj, trim(elemts_C::elemnt[nat[jj]]).c_str(), txtatm[jj].c_str(),
                                     coord[0][jj], coord[1][jj], coord[2][jj]);
                        std::fprintf(stdout, " Atom:%5d: %s (%s) %12.6f %12.6f %12.6f\n",
                                     l, trim(elemts_C::elemnt[nat[l]]).c_str(), txtatm[l].c_str(),
                                     coord[0][l], coord[1][l], coord[2][l]);
                    } else {
                        if (index1(keywrd, " LOCATE-TS") == 0 && index1(keywrd, " 0SCF") == 0) {
                            std::fprintf(stdout,
                                         " The bond defined by CVB between atoms%5d and%5d did not exist,\n but, because GEO-OK is present, the job will continue.\n",
                                         jj, l);
                        }
                    }
                    error = true;
                }
            }
        label99:
            for (j = i; j <= k; ++j)
                if (keywrd[j - 1] == ',' || keywrd[j - 1] == ';') break;
            if (j <= k) { i = j; }
            else break;
        }
        for (i = 1; i <= ncvb; ++i) {
            j = cvb_atoms[i];
            if (nbonds[j] == 0) {
                if (nat[j] == 1) {
                    //  A bond to hydrogen has been broken, so see if there is a different atom within range.
                    bool found = false;
                    for (k = 1; k <= numat; ++k) {
                        bool skip = false;
                        for (l = 1; l <= ncvb; ++l)
                            if (cvb_atoms[l] == k) { skip = true; break; }
                        if (skip) continue;
                        if (distance(j, k) < 1.2) {
                            nbonds[j] = 1;
                            ibonds[1][j] = k;
                            nbonds[k] = nbonds[k] + 1;
                            ibonds[nbonds[k]][k] = j;
                            found = true;
                            break;
                        }
                    }
                    if (nbonds[j] == 1) continue;
                }
                char buf[256];
                if (trim(txtatm[j]) != " ")
                    std::snprintf(buf, sizeof buf, "All bonds to atom %d, label: \"%s\" have been deleted",
                                  j, txtatm[j].c_str());
                else
                    std::snprintf(buf, sizeof buf,
                                  " All bonds to atom %s %5d have been deleted. Coords:%10.5f %10.5f %10.5f",
                                  trim(elemts_C::elemnt[nat[j]]).c_str(), j,
                                  coord[0][j], coord[1][j], coord[2][j]);
                mopend(trim(buf));
                error = true;
            }
        }
    }
    if (error && !let) {
        if (index1(keywrd, " GEO-OK") + index1(keywrd, " 0SCF") == 0) {
            std::fprintf(stdout, "\n          To continue add \"GEO-OK\" or correct reported faults\n");
            mopend("Fault in CVB keyword");
            return;
        } else {
            moperr = false;
            error = false;
        }
    }
}
}

void lewis(bool use_cvs) {
    // Most of the time, LET should be false, the commonest exceptions being 0SCF and RESEQ
    if (!store_coord.empty() && icalcn_lw != numcal) { store_coord.clear(); icalcn_lw = numcal; }
    std::vector<int> iz(natoms + 1);
    if (store_coord.empty()) {
        int ncol = (int)coord[0].size();
        store_coord.assign(3, std::vector<double>(ncol, 0.0));
        for (int i = 0; i <= 2; ++i)
            for (int j = 1; j < ncol; ++j) store_coord[i][j] = coord[i][j];
    }
    bool let = (index1(keywrd, " 0SCF") + index1(keywrd, " RESEQ") + index1(keywrd, " SITE=") +
                index1(keywrd, " LET") + index1(keywrd, " ADD-H") != 0);
    bool debug = (index1(keywrd, " LEWIS") != 0);

    //  Work out raw connectivity.  This depends only on the topology.
    //  Arrays nbonds and ibonds are filled here.
    set_up_dentate();

    //  Put limit on maximum number of bonds to an atom.
    //  Fe, Ni, Pd, Pt
    for (int i = 1; i <= numat; ++i) {
        if (nat[i] == 26 || nat[i] == 28 || nat[i] == 46 || nat[i] == 78) {
            while (nbonds[i] > 5) remove_bond(i);
        }
    }
    if (use_cvs) check_CVS(true);

    //  Check bonding of hydrogen atoms - check for bridges
    int ibad = 0, jbad = 0;
    check_h(ibad);
    if (ibad == 1) std::fprintf(stdout, "\n");
    if (ibad != 0 && !let) {
        std::fprintf(stdout, "\n\n          Add keyword \"LET\" to allow job to continue\n");
        mopend("A hydrogen atom is badly positioned");
    }
    jbad = ibad;
    ibad = 0;

    //  Special check for nitro groups
    for (int i = 1; i <= numat; ++i) {
        if (nat[i] == 7) {
            int k = 0, l = 0, j;
            for (j = 1; j <= nbonds[i]; ++j)
                if (nat[ibonds[j][i]] == 8 && nbonds[ibonds[j][i]] == 1) { k = ibonds[j][i]; break; }
            for (j = j + 1; j <= nbonds[i]; ++j)
                if (nat[ibonds[j][i]] == 8 && nbonds[ibonds[j][i]] == 1) { l = ibonds[j][i]; break; }
            //  Found a system with N bonded to two oxygen atoms (k and l),
            //  which are not bonded to anything else.
            if (k != 0 && l != 0) {
                //  Put a bond between the two oxygen atoms
                nbonds[k] = 2;
                nbonds[l] = 2;
                ibonds[2][k] = l;
                ibonds[2][l] = k;
            }
        }
    }

    //  Put atoms into order of atomic number
    for (int i = 1; i <= numat; ++i) {
        if (nbonds[i] > 1) {
            int n = nbonds[i];
            for (int j = 1; j <= n; ++j) iz[j] = ibonds[j][i];
            for (int j = 1; j <= n; ++j) {
                int k = 0, l = 0;
                for (int m = 1; m <= n; ++m)
                    if (iz[m] != 0 && nat[iz[m]] > k) { l = m; k = nat[iz[m]]; }
                ibonds[j][i] = iz[l];
                iz[l] = 0;
            }
        }
    }
    if (prt_topo && debug && jbad != 0 && !let) {
        int j = (maxtxt == 0) ? 1 : maxtxt / 2 + 2;
        std::string hdr = line.substr(0, j);
        std::fprintf(stdout, "\n   TOPOGRAPHY OF SYSTEM\n\n");
        std::fprintf(stdout, "  ATOM No. %s  LABEL  %sAtoms connected to this atom\n",
                     hdr.c_str(), hdr.c_str());
        for (int i = 1; i <= numat; ++i) {
            if (j == 0) {
                std::fprintf(stdout, "%7d         %s", i, trim(elemts_C::elemnt[nat[i]]).c_str());
            } else {
                if (maxtxt > 2) {
                    std::string lbl = trim(elemts_C::elemnt[nat[i]]) + " (" + txtatm[i].substr(0, maxtxt) + ") ";
                    std::fprintf(stdout, "%7d         %s", i, lbl.c_str());
                } else {
                    std::fprintf(stdout, "%7d         %s", i, trim(elemts_C::elemnt[nat[i]]).c_str());
                }
            }
            for (int jj = 1; jj <= nbonds[i]; ++jj) std::fprintf(stdout, "%7d", ibonds[jj][i]);
            std::fprintf(stdout, "\n");
        }
    }
    if (jbad != 0 && !let) mopend("Geometry is faulty");
    if (let) moperr = false;
}
