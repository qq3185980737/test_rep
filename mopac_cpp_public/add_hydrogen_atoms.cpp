// add_hydrogen_atoms.cpp — C++ translation of MOPAC 2016
// "add_hydrogen_atoms.F90" (~2350 lines of Fortran, 9 subprograms).
//
// Translated faithfully; 1-based Fortran indexing is kept for all module
// arrays (element 0 unused). External routines are forward-declared.

#include "add_hydrogen_atoms.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "mod_atomradii.h"
#include "molkst_C.h"
#include "parameters_C.h"

using common_arrays_C::breaks;
using common_arrays_C::break_coords;
using common_arrays_C::chains;
using common_arrays_C::coord;
using common_arrays_C::coorda;
using common_arrays_C::geo;
using common_arrays_C::ibonds;
using common_arrays_C::labels;
using common_arrays_C::l_atom;
using common_arrays_C::lopt;
using common_arrays_C::loc;
using common_arrays_C::nat;
using common_arrays_C::nfirst;
using common_arrays_C::nbonds;
using common_arrays_C::nlast;
using common_arrays_C::txtatm;
using common_arrays_C::txtatm1;
using common_arrays_C::xparam;
using common_arrays_C::atmass;

namespace mc = common_arrays_C;
namespace mk = molkst_C;

// ---- External routines (defined in other, not-yet-ported files) --------
double distance(int i, int j);                       // interatomic distance
double reada(const std::string& s, int pos);         // parse float from text
void dihed(const std::vector<std::vector<double>>& coord,
           int a, int b, int c, int d, double& sum);  // dihedral -> sum
void bangle(const std::vector<std::vector<double>>& coord,
            int a, int b, int c, double& sum);           // bond angle -> sum
void upcase(std::string& s, int n);                  // uppercase first n chars
void l_control(const std::string& s, int len, int flag);
void geochk();
void lewis(bool flag);
void set_up_dentate();
void mopend(const std::string& msg);
void reset_breaks();  // defined in this file

// ---- File-local helpers -------------------------------------------------
namespace {

// Fortran substring s(a:b) (1-based, inclusive) -> std::string.
inline std::string fsub(const std::string& s, int a, int b) {
    if (a < 1) a = 1;
    if (b > (int)s.size()) b = (int)s.size();
    return s.substr(a - 1, b - a + 1);
}

// Near-one-metal test; defined below.
bool near_a_metal(int icc, const std::vector<int>& metals, int nmetals);

// Work out the type/count of hydrogens for atom icc.
void h_type(int jjj, int icc, bool ionized, int& nH, int& nb_icc, int& nc_icc,
            int& nd_icc, double& bond_length, double& angle, double& dihedral,
            double& internal_dihedral, std::vector<int>& hybrid,
            const std::vector<int>& metals, int nmetals);

void add_a_generic_hydrogen_atom(int na, int nb, int nc, double bond_length,
                                 double angle, double dihedral,
                                 const std::vector<int>& metals, int nmetals);

void add_a_sp3_hydrogen_atom(int icc, int nb_icc, int nc_icc, int nd_icc,
                             double bond_length,
                             const std::vector<int>& metals, int nmetals);

bool aromatic(int atom_1, int atom_2, int atom_6);
bool aromatic_5(int atom_1, int atom_2, int atom_5, std::vector<int>& hybrid);

}  // namespace

// ===========================================================================
//  subroutine add_hydrogen_atoms()
// ===========================================================================
void add_hydrogen_atoms() {
    int i, j, k, l, ii;
    int nH, icc, store_numat, nb_icc, nc_icc, nd_icc, nmetals;
    double bond_length, angle, dihedral, internal_dihedral;
    bool ionized = false;

    // i_add(20,14): number of hydrogens on each atom in a residue. 1-based.
    static int i_add[21][15] = {
        {0},
        /* 1 GLY */ {0, 1, 2, 0, 0},
        /* 2 ALA */ {0, 1, 1, 0, 0, 3},
        /* 3 SER */ {0, 1, 1, 0, 0, 2, 1},
        /* 4 CYS */ {0, 1, 1, 0, 0, 2, 1},
        /* 5 VAL */ {0, 1, 1, 0, 0, 1, 3, 3},
        /* 6 THR */ {0, 1, 1, 0, 0, 1, 1, 3},
        /* 7 ILE */ {0, 1, 1, 0, 0, 1, 2, 3, 3},
        /* 8 PRO */ {0, 0, 1, 0, 0, 2, 2, 2},
        /* 9 MET */ {0, 1, 1, 0, 0, 2, 2, 0, 3},
        /*10 ASP */ {0, 1, 1, 0, 0, 2, 0, 0, 1},
        /*11 ASN */ {0, 1, 1, 0, 0, 2, 0, 0, 2},
        /*12 LEU */ {0, 1, 1, 0, 0, 2, 1, 3, 3},
        /*13 LYS */ {0, 1, 1, 0, 0, 2, 2, 2, 2, 2},
        /*14 GLU */ {0, 1, 1, 0, 0, 2, 2, 0, 0, 1},
        /*15 GLN */ {0, 1, 1, 0, 0, 2, 2, 0, 0, 2},
        /*16 ARG */ {0, 1, 1, 0, 0, 2, 2, 2, 1, 0, 2, 1},
        /*17 HIS */ {0, 1, 1, 0, 0, 2, 0, 1, 1, 1, 0},
        /*18 PHE */ {0, 1, 1, 0, 0, 2, 0, 1, 1, 1, 1, 1},
        /*19 TYR */ {0, 1, 1, 0, 0, 2, 0, 1, 1, 1, 1, 0, 1},
        /*20 TRP */ {0, 1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 1, 1},
    };

    std::vector<int> store_nat, store_labels, store_lopt1;
    std::vector<char> store_l_atom;
    std::vector<double> store_coord1, store_atmass;
    std::vector<std::string> store_txtatm, store_txtatm1;
    std::vector<double> store_atom_radius_covalent(107, 0.0);
    std::vector<int> metals, hybrid;

    // Ionized-form overrides.
    if (mk::keywrd.find(" IONI") != std::string::npos) {
        i_add[13][9] = 3;   // Lysine nitrogen (+)
        i_add[16][11] = 2;  // Arginine nitrogen (+)
        i_add[17][10] = 1;  // Histidine nitrogen (+)
        i_add[10][8] = 0;   // Aspartate oxygen (-)
        i_add[14][9] = 0;   // Glutamate oxygen (-)
        ionized = true;
    }

    int nm = mk::numat;
    store_nat.resize(nm + 1);
    store_labels.resize(nm + 1);
    store_atmass.resize(nm + 1);
    store_l_atom.resize(nm + 1);
    store_txtatm.resize(nm + 1);
    store_txtatm1.resize(nm + 1);
    store_atmass.resize(nm + 1);
    metals.resize(nm + 1);
    // store_coord(3,numat) / store_lopt(3,numat)
    std::vector<std::vector<double>> store_coord(4, std::vector<double>(nm + 1, 0.0));
    std::vector<std::vector<int>> store_lopt(4, std::vector<int>(nm + 1, 0));

    // Remove existing hydrogen atoms.
    j = 0;
    for (i = 1; i <= mk::numat; ++i) {
        if (nat[i] != 1) {
            j = j + 1;
            store_nat[j] = nat[i];
            store_coord[0][j] = coord[0][i];
            store_coord[1][j] = coord[1][i];
            store_coord[2][j] = coord[2][i];
            store_atmass[j] = atmass[i];
            store_lopt[1][j] = lopt[1][i];
            store_lopt[2][j] = lopt[2][i];
            store_lopt[3][j] = lopt[3][i];
            store_l_atom[j] = l_atom[i];
            store_txtatm[j] = txtatm[i];
            store_txtatm1[j] = txtatm1[i];
            txtatm1[j] = txtatm1[i];
            txtatm[j] = txtatm[i];
            coorda[0][j] = coorda[0][i];
            coorda[1][j] = coorda[1][i];
            coorda[2][j] = coorda[2][i];
        }
    }
    if (i != j) reset_breaks();
    if (mk::moperr) return;

    k = 0;
    for (i = 1; i <= mk::natoms; ++i) {
        if (labels[i] != 99 && labels[i] != 1) {
            k = k + 1;
            store_labels[k] = labels[i];
        }
    }
    mk::numat = j;
    mk::natoms = mk::numat;
    nmetals = 0;
    store_atom_radius_covalent = mod_atomradii::atom_radius_covalent;
    for (i = 1; i <= mk::numat; ++i) {
        switch (store_nat[i]) {
            case 3: case 4: case 11: case 12:
            case 19: case 20: case 21: case 22: case 23: case 24:
            case 25: case 26: case 27: case 28: case 29: case 30:
            case 37: case 38: case 39: case 40: case 41: case 42:
            case 43: case 44: case 45: case 46: case 47: case 48:
            case 55: case 56: case 57: case 58: case 59: case 60:
            case 61: case 62: case 63: case 64: case 65: case 66:
            case 67: case 68: case 69: case 70: case 71: case 72:
            case 73: case 74: case 75: case 76: case 77: case 78: case 79: case 80:
                nmetals = nmetals + 1;
                metals[nmetals] = i;
                mod_atomradii::atom_radius_covalent[store_nat[i]] = -1.0;
                break;
        }
    }

    // Make extra storage room.
    mk::maxatoms = mk::numat * 5;
    int mx = mk::maxatoms;
    nat.assign(mx + 1, 0);
    txtatm.assign(mx + 1, "");
    txtatm1.assign(mx + 1, "");
    coord.assign(4, std::vector<double>(mx + 1, 0.0));
    atmass.assign(mx + 1, 0.0);
    nbonds.assign(mx + 1, 0);
    ibonds.assign(16, std::vector<int>(mx + 1, 0));
    labels.assign(mx + 1, 0);
    geo.assign(4, std::vector<double>(mx + 1, 0.0));
    lopt.assign(4, std::vector<int>(mx + 1, 0));
    mc::na.assign(mx + 1, 0);
    l_atom.assign(mx + 1, 0);
    loc.assign(3, std::vector<int>(3 * mx + 1, 0));
    xparam.assign(3 * mx + 1, 0.0);
    nfirst.assign(mx + 1, 0);
    nlast.assign(mx + 1, 0);

    for (i = 1; i <= mk::numat; ++i) {
        nat[i] = store_nat[i];
        txtatm[i] = store_txtatm[i];
        txtatm1[i] = (std::string(26, ' '));
        txtatm1[i] = store_txtatm1[i];
        coord[0][i] = store_coord[0][i];
        coord[1][i] = store_coord[1][i];
        coord[2][i] = store_coord[2][i];
        labels[i] = store_labels[i];
        atmass[i] = store_atmass[i];
        lopt[1][i] = store_lopt[1][i];
        lopt[2][i] = store_lopt[2][i];
        lopt[3][i] = store_lopt[3][i];
        l_atom[i] = store_l_atom[i];
    }
    store_numat = mk::numat;

    geochk();
    if (mk::moperr) return;

    mk::line = mk::refkey.empty() ? "" : mk::refkey[1];
    upcase(mk::line, (int)mk::line.size());
    i = (int)mk::line.find(" CVB");
    if (i >= 0) {
        // j = index of ") " in line(i+1:) ; Fortran index() is 1-based.
        std::string tail = mk::line.substr(i);  // line(i:) 1-based i+1
        j = (int)tail.find(") ");
        if (j >= 0) {
            // line = line(i+1 : j-1) in 1-based: i+1 .. j-1 where j is 1-based position.
            int lo = i + 1;                       // 0-based start of "(" region
            int hi = i + j - 1;                   // exclusive-ish
            mk::line = mk::line.substr(lo, hi - lo);
            ii = (int)mk::line.size();
            k = 0;
            while (true) {
                k = k + 1;
                if (k >= ii) break;
                if (mk::line[k - 1] == '"') {
                    while (true) {
                        k = k + 1;
                        if (mk::line[k - 1] == '"') break;
                    }
                }
                if (mk::line[k - 1] >= '0' && mk::line[k - 1] <= '9') break;
            }
            if (k > ii - 3) l_control(mk::line, ii, 1);
        }
    }
    lewis(false);

    // Two-pass hybridization.
    hybrid.assign(mk::numat + 1, 0);
    for (icc = 1; icc <= mk::numat; ++icc) {
        h_type(22, icc, ionized, nH, nb_icc, nc_icc, nd_icc, bond_length, angle,
               dihedral, internal_dihedral, hybrid, metals, nmetals);
    }
    for (icc = 1; icc <= mk::numat; ++icc) {
        bond_length = 0.0;
        angle = 0.0;
        j = 1;
        for (; j <= 20; ++j) {
            if (MOZYME_C::afn[j] == fsub(txtatm[icc], 18, 20)) break;
        }
        if (j == 17 && nat[icc] == 7) j = 21;  // His imidazole N

        if (nat[icc] == 8 && nbonds[icc] == 1) {  // carboxylic acid
            bool cooh_done = false;
            for (i = 1; i <= nbonds[icc] && !cooh_done; ++i) {
                k = ibonds[i][icc];
                if (nat[k] == 6 && nbonds[k] == 3) {
                    l = 0;
                    for (ii = 1; ii <= nbonds[k]; ++ii) {
                        if (nat[ibonds[ii][k]] == 8) l = l + 1;
                    }
                    if (l == 2) {
                        j = 21;
                        cooh_done = true;
                    }
                }
            }
        }
        if (nat[icc] == 7) {  // arginine
            bool arg_done = false;
            for (i = 1; i <= nbonds[icc] && !arg_done; ++i) {
                k = ibonds[i][icc];
                if (nat[k] == 6 && nbonds[k] == 3) {
                    l = 0;
                    for (ii = 1; ii <= nbonds[k]; ++ii) {
                        if (nat[ibonds[ii][k]] == 7) l = l + 1;
                    }
                    if (l == 3) {
                        j = 21;
                        arg_done = true;
                    }
                }
            }
        }

        h_type(j, icc, ionized, nH, nb_icc, nc_icc, nd_icc, bond_length, angle,
               dihedral, internal_dihedral, hybrid, metals, nmetals);

        if (j < 21) {
            bond_length = 0.0;
            for (k = 1; k <= MOZYME_C::n_add[j]; ++k) {
                std::string labl2 = elemts_C::elemnt[nat[icc]] +
                                    fsub(txtatm[icc], 15, 15) +
                                    fsub(txtatm[icc], 16, 16);
                if (labl2 == MOZYME_C::atomname[j][k]) break;
            }
            if (k > MOZYME_C::n_add[j]) {
                if (fsub(txtatm[icc], 15, 18) == " OXT") {
                    if (ionized) nH = 0;
                    else nH = 1;
                } else {
                    nH = 0;
                }
            } else {
                nH = i_add[j][k];
                if (k == 1 && nbonds[icc] == 1) {
                    if (ionized) nH = 3;
                    else nH = 2;
                }
                if (k == 3 && nbonds[icc] == 2) nH = 1;
                if (j == 3 && k == 6 && nbonds[icc] > 1) nH = 0;
                if (j == 4 && k == 6 && nbonds[icc] > 1) nH = 0;
                if (j == 8 && k == 1 && nbonds[icc] == 2) nH = 1;
                if (j == 11 && k == 8 && nbonds[icc] == 2) nH = 1;
                if (nat[icc] == 6 && k > 4 && nbonds[icc] == 1) nH = 3;
            }
            switch (nH) {
                case 1:
                    if (nd_icc == 0) {
                        angle = 120.0;
                        dihedral = 180.0;
                    } else {
                        bond_length = 1.09;
                    }
                    break;
                case 2:
                    angle = 109.0;
                    dihedral = 120.0;
                    internal_dihedral = 120.0;
                    break;
                case 3:
                    angle = 109.471221;
                    dihedral = 120.0;
                    internal_dihedral = 120.0;
                    break;
            }
            if (nat[icc] == 6) bond_length = 1.09;
            else if (nat[icc] == 7) bond_length = 1.01;
            else bond_length = 0.99;
        }
        if (std::abs(bond_length) < 0.1) continue;
        angle = angle * funcon_C::pi / 180.0;
        dihedral = dihedral * funcon_C::pi / 180.0;
        internal_dihedral = internal_dihedral * funcon_C::pi / 180.0;

        switch (nH) {
            case 1:
                if (nc_icc != 0 && nd_icc != 0) {
                    add_a_sp3_hydrogen_atom(icc, nb_icc, nc_icc, nd_icc,
                                            bond_length, metals, nmetals);
                } else {
                    add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length,
                                                angle, dihedral, metals, nmetals);
                }
                break;
            case 2:
                add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length,
                                            angle, dihedral, metals, nmetals);
                internal_dihedral = internal_dihedral + dihedral;
                add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length,
                                            angle, internal_dihedral, metals, nmetals);
                break;
            case 3:
                add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length,
                                            angle, funcon_C::pi, metals, nmetals);
                add_a_generic_hydrogen_atom(icc, nb_icc, mk::numat, bond_length,
                                            angle, internal_dihedral, metals, nmetals);
                add_a_generic_hydrogen_atom(icc, nb_icc, mk::numat, bond_length,
                                            angle, internal_dihedral, metals, nmetals);
                break;
        }
    }

    // Set all new atoms to hydrogen.
    for (i = store_numat + 1; i <= mk::numat; ++i) {
        nat[i] = 1;
        labels[i] = 1;
        atmass[i] = parameters_C::ams[1];
        j = ibonds[1][i];
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%-6s%5d  H  %s",
                      txtatm[j].substr(0, 6).c_str(), i,
                      txtatm[j].substr(16).c_str());
        txtatm[i] = buf;
    }

    mk::natoms = std::max(mk::numat, mk::natoms);
    for (i = 1; i <= mk::numat; ++i) {
        geo[1][i] = coord[0][i];
        geo[2][i] = coord[1][i];
        geo[3][i] = coord[2][i];
    }
    for (i = store_numat + 1; i <= mk::numat; ++i) {
        lopt[1][i] = 1;
        lopt[2][i] = 1;
        lopt[3][i] = 1;
        l_atom[i] = 1;
    }
    for (i = 1; i <= mk::numat; ++i) mc::na[i] = 0;
    mod_atomradii::atom_radius_covalent = store_atom_radius_covalent;

    if (!mk::isok) {
        set_up_dentate();
        geochk();
    }
    if (mk::keywrd.find("SITE") != std::string::npos) geochk();
    reset_breaks();
}

// ===========================================================================
//  subroutine h_type(...)
// ===========================================================================
namespace {

void h_type(int jjj, int icc, bool ionized, int& nH, int& nb_icc, int& nc_icc,
            int& nd_icc, double& bond_length, double& angle, double& dihedral,
            double& internal_dihedral, std::vector<int>& hybrid,
            const std::vector<int>& metals, int nmetals) {
    int i=0, j=0, k=0, l=0, kk=0, ll=0, ii=0, jj=0;
    double Rab = 0, Rac = 0, Rad = 0, sum = 0;
    bool l_tmp;
    const double cyanide = 1.3, CC = 1.48;
    // search(3,6): 1-based columns i=1..6 direction table.
    int search[4][7] = {
        {0},
        {0, 1, -1, 0, 0, 0, 0},  // row 1
        {0, 0, 0, 1, -1, 0, 0},  // row 2
        {0, 0, 0, 0, 0, 1, -1},  // row 3
    };

    bond_length = 0.0;
    angle = 0.0;
    nb_icc = 0;
    nc_icc = 0;
    nd_icc = 0;
    nH = 0;

    if (nbonds[icc] == 0) {
        // isolated atom
        if (hybrid[icc] == 0) {
            hybrid[icc] = 3;
            return;
        }
        switch (nat[icc]) {
            case 6:
                nH = 3; bond_length = 1.09; angle = 109.471221;
                dihedral = 120; internal_dihedral = 120; break;
            case 7:
                nH = 2; bond_length = 1.0; angle = 107;
                dihedral = 114.5; internal_dihedral = dihedral; break;
            case 8:
                nH = 1; bond_length = 0.957; angle = 104.4; break;
            case 14:
                nH = 3; bond_length = 1.48; angle = 109.471221;
                dihedral = 120; internal_dihedral = 120; break;
            case 15:
                nH = 2; bond_length = 1.42; angle = 93.5;
                dihedral = 93.8; internal_dihedral = dihedral; break;
            case 16:
                nH = 1; bond_length = 1.336; angle = 92.1; break;
            case 9: case 17: case 35: case 53:
                return;
            case 33:
                nH = 2; bond_length = 1.519; angle = 91.8;
                dihedral = 92; internal_dihedral = dihedral; break;
            case 34:
                nH = 1; bond_length = 1.46; angle = 91; break;
            case 51:
                nH = 2; bond_length = 1.707; angle = 91.8;
                dihedral = 92; internal_dihedral = dihedral; break;
            case 52:
                nH = 1; bond_length = 1.69; angle = 90; break;
            default:
                nH = 0; return;
        }
        mk::numat = mk::numat + 1;
        nat[mk::numat] = 1;
        labels[mk::numat] = 1;
        for (i = 1; i <= 6; ++i) {
            coord[0][mk::numat] = coord[0][icc] + bond_length * search[1][i];
            coord[1][mk::numat] = coord[1][icc] + bond_length * search[2][i];
            coord[2][mk::numat] = coord[2][icc] + bond_length * search[2][i];
            if (!near_a_metal(icc, metals, nmetals)) break;
        }
        nb_icc = mk::numat;
        nbonds[icc] = 1;
        nbonds[mk::numat] = 1;
        ibonds[1][icc] = mk::numat;
        ibonds[1][mk::numat] = icc;
        return;
    }

    ii = ibonds[1][icc];
    if (nbonds[icc] == 3) {
        kk = ibonds[3][icc];
        jj = ibonds[2][icc];
        if (nat[ii] > nat[jj]) { i = jj; jj = ii; ii = i; }
        if (nat[ii] > nat[kk]) { i = kk; kk = ii; ii = i; }
        if (nat[jj] > nat[kk]) { i = kk; kk = jj; jj = i; }
        Rab = distance(icc, ii);
        Rac = distance(icc, jj);
        Rad = distance(icc, kk);
    } else if (nbonds[icc] == 2) {
        jj = ibonds[2][icc];
        if (nat[jj] < nat[ii]) { i = ii; ii = jj; jj = i; }
        Rab = distance(icc, ii);
        Rac = distance(icc, jj);
        if (nat[ii] == 6 && nat[jj] == 6) {
            if (Rab > Rac) {
                Rad = Rab; Rab = Rac; Rac = Rad;
                i = ii; ii = jj; jj = i;
            }
        }
    } else {
        Rab = distance(icc, ii);
    }

    bond_length = 1.09;
    if (nbonds[icc] == 1) {
        nb_icc = ii;
        for (i = 1; i <= nbonds[ii]; ++i) {
            kk = ibonds[i][ii];
            if (kk != icc && nat[kk] > 1) nc_icc = kk;
        }
        if (nc_icc == 0) {
            for (i = 1; i <= nbonds[ii]; ++i) {
                kk = ibonds[i][ii];
                if (kk != icc) nc_icc = kk;
            }
        }
    } else if (nbonds[icc] == 2) {
        nb_icc = ii;
        nc_icc = jj;
    } else {
        nb_icc = ibonds[1][icc];
        nc_icc = ibonds[2][icc];
        nd_icc = ibonds[3][icc];
    }
    if (jjj < 21) return;

    switch (nat[icc]) {
        case 5:
            if (nbonds[icc] == 1) {
                nH = 3; bond_length = 1.07; angle = 109.471;
                internal_dihedral = 120; dihedral = 180; hybrid[icc] = 3;
            } else if (nbonds[icc] == 2) {
                nH = 1; bond_length = 1.07; angle = 120.0;
                dihedral = 180; hybrid[icc] = 1;
            }
            break;

        case 6: {  // Carbon
            switch (nbonds[icc]) {
                case 1: {  // bonded to one atom
                    if (nat[ii] == 6) {
                        jj = 0;
                        for (i = 1; i <= nbonds[ii]; ++i)
                            if (nat[ibonds[i][ii]] == 6) jj = jj + 1;
                        if (jj == 3) {
                            jj = 0; kk = 0;
                            for (i = 1; i <= nbonds[ii]; ++i) {
                                ll = ibonds[i][ii];
                                if (nat[ll] == 6 && jj == 0 && ll != icc && (kk != ll || kk == 0)) jj = ibonds[i][ii];
                                if (nat[ll] == 6 && kk == 0 && ll != icc && (jj != ll || jj == 0)) kk = ibonds[i][ii];
                            }
                            dihed(coord, icc, ii, jj, kk, sum);
                            sum = std::fabs(sum);
                            if (hybrid[ii] == 3 || std::min(funcon_C::pi - sum, sum) > 0.4) {
                                nH = 3; angle = 109.5; internal_dihedral = 120; hybrid[icc] = 3; return;
                            }
                        }
                        jj = 0;
                        for (i = 1; i <= nbonds[ii]; ++i) {
                            kk = ibonds[i][ii];
                            if (kk != icc) jj = kk;
                        }
                        if (Rab > CC) {
                            nH = 3; angle = 109.471; internal_dihedral = 120; dihedral = 180; hybrid[icc] = 3;
                            if (nc_icc == 0) {
                                for (i = 1; i <= nbonds[ii]; ++i) {
                                    j = ibonds[i][ii];
                                    if (j != icc) { nc_icc = j; return; }
                                }
                            }
                        } else if (Rab > 1.35) {
                            nH = 3; angle = 109.471; internal_dihedral = 120; dihedral = 180; hybrid[icc] = 3;
                        } else if (Rab > 1.25) {
                            nH = 2; bond_length = 1.08; angle = 123; dihedral = 180; internal_dihedral = 180; hybrid[icc] = 2;
                        } else {
                            nH = 1; bond_length = 1.07; angle = 179; hybrid[icc] = 1;
                        }
                    } else if (nat[ii] == 7) {
                        if (hybrid[ii] == 3) {
                            nH = 3; bond_length = 1.09; angle = 109.471; internal_dihedral = 120; return;
                        }
                        hybrid[icc] = 1;
                        if (ionized) nH = 0;
                        else { nH = 1; angle = 179; }
                    } else if (nat[ii] == 8) {
                        if (Rab > 1.37 || nbonds[ii] > 1) {
                            nH = 3; bond_length = 1.09; angle = 109.471; internal_dihedral = 120; hybrid[icc] = 3;
                        } else {
                            nH = 2; angle = 120; internal_dihedral = 180; hybrid[icc] = 2;
                        }
                    } else {
                        nH = 3; angle = 109.471; internal_dihedral = 120; dihedral = 180; hybrid[icc] = 3;
                    }
                    break;
                }
                case 2: {  // bonded to two atoms
                    if (nat[ii] == 6) {
                        if (nat[jj] == 6) {
                            if (aromatic_5(icc, ii, jj, hybrid)) {
                                if (Rab < 1.44) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                                else { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return; }
                            }
                            if (aromatic(icc, ii, jj)) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                            if (Rab < 1.27) { nH = 0; hybrid[icc] = 1; return; }
                            nH = 2; angle = 109; internal_dihedral = 120; dihedral = 120; hybrid[icc] = 3;
                            if (hybrid[ii] == 3 && hybrid[jj] == 3) return;
                            bangle(coord, ii, icc, jj, sum);
                            if (sum > 3.0 && Rac < 1.35) { nH = 0; hybrid[icc] = 2; return; }
                            if (Rab < 1.35) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                            i = 0;
                            if ((sum > 2.03 && Rab < 1.39) || (sum > 2.2 && Rab < 1.41)) {
                                if (nbonds[ii] != 4) {
                                    if (nbonds[ii] == 3) {
                                        j = 0;
                                        for (k = 1; k <= 3; ++k)
                                            if (nat[ibonds[k][ii]] == 8) j = j + 1;
                                        if (j != 2) i = 1;
                                    } else if (nbonds[ii] == 2) {
                                        bangle(coord, ibonds[1][ii], ii, ibonds[2][ii], sum);
                                        if (sum > 2.0) i = 1;
                                    } else {
                                        i = 1;
                                    }
                                }
                            }
                            if (i == 1) {
                                nH = 1; angle = 120; internal_dihedral = 180; dihedral = 180; hybrid[icc] = 2;
                            }
                        } else if (nat[jj] == 7) {
                            if (Rac < 1.2) { hybrid[icc] = 1; return; }
                            else if (hybrid[ii] == 3 && hybrid[jj] == 3) {
                                nH = 2; bond_length = 1.09; angle = 109; internal_dihedral = 120; dihedral = 120; hybrid[icc] = 3;
                            } else if (Rab < 1.43 || ((Rab + Rac < 2.9) && aromatic_5(icc, ii, jj, hybrid) && nbonds[jj] == 2)) {
                                nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2;
                            } else if (Rab + Rac < 2.8 && aromatic_5(icc, ii, jj, hybrid)) {
                                nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2;
                            } else if (aromatic(icc, ii, jj)) {
                                nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return;
                            } else {
                                nH = 2; bond_length = 1.09; angle = 109; internal_dihedral = 120; dihedral = 120; hybrid[icc] = 3;
                            }
                        } else if (nat[jj] == 8) {
                            if (aromatic_5(icc, ii, jj, hybrid)) {
                                if (Rab < 1.44) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                                else { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return; }
                            }
                            if (Rab < 1.4 || (Rac < 1.3 && nbonds[jj] == 1)) {
                                if (Rab < 1.4 && (Rac < 1.3 && nbonds[jj] == 1)) return;
                                nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return;
                            }
                            sum = 0.0; j = 0;
                            for (i = 1; i <= nbonds[ii]; ++i)
                                if (nat[ibonds[i][ii]] != 1) j = j + 1;
                            if (nat[ii] == 6 && j == 2) {
                                jj = ibonds[1][ii];
                                if (jj == icc) jj = ibonds[2][ii];
                                bangle(coord, icc, ii, jj, sum);
                            }
                            if (Rac < 1.3 && sum > 2.03) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                            else { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return; }
                        } else if (nat[jj] == 9 || nat[jj] == 17 || nat[jj] == 35 || nat[jj] == 53) {
                            if (Rab < 1.25) { nH = 0; hybrid[icc] = 1; return; }
                            else if (Rab < 1.44) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                            else { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return; }
                        } else if (nat[jj] == 15) {
                            nH = 2; bond_length = 1.09; angle = 109; internal_dihedral = 120; dihedral = 120; hybrid[icc] = 3;
                        } else if (nat[jj] == 16 || nat[jj] == 34) {
                            if (Rab > 1.35) {
                                nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120;
                                hybrid[icc] = 3; hybrid[icc] = 2;
                            } else if (Rab > 1.25) {
                                if (Rab < 1.4 && (Rac < 1.6 && nbonds[jj] == 1)) return;
                                nH = 1; bond_length = 1.07; angle = 120; dihedral = 180; hybrid[icc] = 1;
                            } else {
                                nH = 0; hybrid[icc] = 3;
                            }
                        }
                    } else if (nat[ii] == 7) {
                        if (nat[jj] == 8) {
                            if (Rac < 1.3) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return; }
                            else { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return; }
                        } else if (nat[jj] == 7) {
                            if (((hybrid[ii] == 3 || hybrid[jj] == 3) && (Rab > 1.4 && Rac > 1.4))) {
                                nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return;
                            } else {
                                nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; return;
                            }
                        } else if (nat[jj] == 14) {
                            if (Rab > 1.4 && Rac > 1.8) { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; }
                            else { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; }
                            return;
                        }
                        if (nat[jj] == 16) {
                            nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return;
                        }
                        if (nbonds[icc] == 2 && Rab < cyanide) nH = 0;
                        else { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; }
                    } else if (nat[ii] == 8) {
                        if (nat[jj] == 9 || nat[jj] == 17 || nat[jj] == 35 || nat[jj] == 53) {
                            if (Rab > 1.3) { nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; }
                            else { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2; }
                            return;
                        }
                        if (nbonds[ii] == 2 && nbonds[jj] == 2) {
                            nH = 2; angle = 109; dihedral = 120; internal_dihedral = 120; hybrid[icc] = 3; return;
                        }
                        bangle(coord, ibonds[1][icc], icc, ibonds[2][icc], sum);
                        if (nbonds[ii] != 1 || nbonds[jj] != 1 || Rab > 1.28 || (nat[jj] == 8 && sum < 2.618)) {
                            nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 2;
                        } else {
                            nH = 0; hybrid[icc] = 2;
                        }
                        hybrid[icc] = 2;
                    } else if (nat[ii] == 15) {
                        nH = 0; hybrid[icc] = 2;
                    } else {
                        hybrid[icc] = 3;
                    }
                    break;
                }
                case 3: {  // bonded to three atoms
                    dihed(coord, icc, ii, jj, kk, sum);
                    if (std::min(2 * funcon_C::pi - sum, sum) < 0.15) {
                        if (aromatic(icc, ii, jj) || aromatic(icc, ii, kk) || aromatic(icc, jj, kk)) {
                            nH = 0; hybrid[icc] = 2; return;
                        }
                        if (aromatic_5(icc, ii, jj, hybrid) || aromatic_5(icc, ii, kk, hybrid) || aromatic_5(icc, jj, kk, hybrid)) {
                            nH = 0; hybrid[icc] = 2; return;
                        }
                    }
                    if (std::min(2 * funcon_C::pi - sum, sum) < 0.10) {
                        if (Rab + Rac + Rad < 4.5) { nH = 0; hybrid[icc] = 2; return; }
                    }
                    if (std::min(2 * funcon_C::pi - sum, sum) > 0.35) {
                        nH = 1; hybrid[icc] = 3; return;
                    }
                    nH = 1; hybrid[icc] = 3;
                    if (nat[ii] == 6) {
                        if (nat[jj] == 6) {
                            if (nat[kk] == 6) {
                                sum = std::min({Rab, Rac, Rad});
                                if (sum < 1.4) { nH = 0; hybrid[icc] = 2; }
                            } else if (nat[kk] == 7) {
                                if (hybrid[ii] == 3 && hybrid[jj] == 3 && hybrid[kk] == 3) return;
                                j = 0;
                                for (i = 1; i <= nbonds[kk]; ++i)
                                    if (nat[ibonds[i][kk]] != 1) j = j + 1;
                                if (std::min(Rab, Rac) < 1.41 || (Rad < 1.39 && j < 3)) { nH = 0; hybrid[icc] = 2; }
                            } else if (nat[kk] == 8) {
                                if (std::min(Rab, Rac) < 1.45 || Rad < 1.30) { nH = 0; hybrid[icc] = 2; }
                            } else if (nat[kk] == 9 || nat[kk] == 17 || nat[kk] == 35 || nat[kk] == 53) {
                                if (std::min(Rab, Rac) < 1.45) { nH = 0; hybrid[icc] = 2; }
                            } else if (nat[kk] == 16) {
                                sum = std::min(Rab, Rac);
                                if (aromatic(icc, ii, jj) || sum < 1.44) { nH = 0; hybrid[icc] = 2; return; }
                            }
                        } else if (nat[jj] == 7) {
                            if (nat[kk] == 7) {
                                sum = std::min(Rab, Rac);
                                if (sum < 1.405 || hybrid[ii] == 2 || hybrid[jj] == 2 || hybrid[kk] == 2) { nH = 0; hybrid[icc] = 2; }
                            } else if (nat[kk] == 8) {
                                if (Rac > 1.45 && Rad > 1.28) return;
                                else { nH = 0; hybrid[icc] = 2; }
                            } else if (nat[kk] == 16) {
                                dihed(coord, icc, ii, jj, kk, sum);
                                if (std::min(2 * funcon_C::pi - sum, sum) < 0.15) { nH = 0; hybrid[icc] = 2; return; }
                            } else {
                                nH = 0; hybrid[icc] = 2;
                            }
                        } else if (nat[jj] == 8) {
                            if (nat[kk] == 8 || Rac < 1.25) {
                                if (nbonds[jj] == 2 && nbonds[kk] == 2) {
                                    if (Rab < 1.44) nH = 0; return;
                                } else { nH = 0; hybrid[icc] = 2; }
                            }
                        } else if (nat[jj] == 9 || nat[jj] == 17 || nat[jj] == 35 || nat[jj] == 53) {
                            if (Rab < 1.44) nH = 0;
                        }
                    } else if (nat[ii] == 7) {
                        if (nat[jj] == 7) {
                            switch (nat[kk]) {
                                case 7: nH = 0; hybrid[icc] = 2; break;
                                case 8: if (Rad < 1.33) { nH = 0; hybrid[icc] = 2; } break;
                                case 16: if (Rad < 1.7) { nH = 0; hybrid[icc] = 2; } break;
                            }
                        }
                    } else if (nat[ii] == 8 || nat[ii] == 16) {
                        nH = 0;
                    }
                    break;
                }
                default:
                    hybrid[icc] = 3;
                    break;
            }
            break;
        }  // case 6

        case 7: {  // Nitrogen
            if (nbonds[icc] == 1) {
                nH = 2; bond_length = 1.02; angle = 109; dihedral = 120; internal_dihedral = 120;
                ii = ibonds[1][icc];
                if (nat[ii] == 6) {
                    if (Rab < cyanide && nbonds[ii] < 3) { nH = 0; hybrid[icc] = 1; return; }
                    sum = 1.29;
                    if (nat[ii] == 6 && nbonds[ii] == 3) {
                        for (i = 1; i <= 3; ++i)
                            if (nat[ibonds[i][ii]] == 8) sum = 1.1;
                    }
                    if (Rab < sum) { nH = 1; angle = 120; internal_dihedral = 180; dihedral = 0; hybrid[icc] = 1; return; }
                    for (i = 1; i <= nbonds[ii]; ++i) {
                        jj = ibonds[i][ii];
                        if (nat[jj] != 6 && nat[jj] != 7) break;
                        if (nat[jj] == 7 && jj != icc) {
                            k = 0; l = 0;
                            for (j = 1; j <= nbonds[jj]; ++j) {
                                if (nat[ibonds[j][jj]] == 6) k = k + 1;
                                if (nat[ibonds[j][jj]] == 1) l = l + 1;
                            }
                            if (k == 2) continue;
                            if (l == 2) {
                                hybrid[icc] = 2;
                                if (ionized) { nH = 2; angle = 120; internal_dihedral = 180; dihedral = 0; }
                                else { nH = 1; angle = 120; internal_dihedral = 180; dihedral = 0; }
                                return;
                            } else {
                                nH = 2; angle = 120; internal_dihedral = 180; dihedral = 0; hybrid[icc] = 3; return;
                            }
                        }
                    }
                }
                if (ionized) {
                    nH = 3;
                    if (nat[ii] == 6) {
                        for (i = 1; i <= nbonds[ii]; ++i)
                            if (nat[ibonds[i][ii]] == 8) break;
                        if (i <= nbonds[ii]) nH = 2;
                    }
                } else {
                    nH = 2; angle = 120; internal_dihedral = 180; dihedral = 0; hybrid[icc] = 3;
                }
            } else if (nbonds[icc] == 2) {
                if (nat[ii] == 7 && nat[jj] == 7) {
                    sum = std::min(Rab, Rac);
                    if (sum < 1.35) { nH = 0; hybrid[icc] = 2; }
                    else { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 3; }
                    return;
                }
                if (nat[ii] == 6 && nat[jj] == 7 && (Rab < 1.36 || Rac < 1.35)) { nH = 0; return; }
                if (nat[ii] == 8 && nat[jj] == 8 && (Rab < 1.30 || Rac < 1.30)) { nH = 0; return; }
                if (aromatic(icc, ii, jj)) { nH = 0; hybrid[icc] = 2; return; }
                kk = 0;
                for (i = 1; i <= nbonds[ii]; ++i) {
                    j = ibonds[i][ii];
                    if (nat[j] == 7 && j != icc) kk = j;
                }
                if (kk == 0) {
                    for (i = 1; i <= nbonds[jj]; ++i) {
                        j = ibonds[i][jj];
                        if (nat[j] == 7 && j != icc) { ii = jj; kk = j; }
                    }
                }
                if (nat[ii] == 6) {
                    i = nbonds[ii];
                    for (j = 1; j <= i; ++j) {
                        k = ibonds[j][ii];
                        if (nat[k] == 8) {
                            Rab = distance(ii, k);
                            if (Rab < 1.33) kk = 0;
                        }
                    }
                }
                if (kk != 0) {
                    l = 0;
                    if (nbonds[kk] >= 2) {
                        ii = ibonds[1][kk];
                        if (nat[ii] == 1 || nat[ii] == 6) l = 1;
                        ii = ibonds[2][kk];
                        if (nat[ii] == 1 || nat[ii] == 6) l = l + 1;
                        if (nbonds[kk] == 3) {
                            ii = ibonds[3][kk];
                            if (nat[ii] == 1 || nat[ii] == 6) l = l + 1;
                        }
                        if (l == 3) {
                            hybrid[icc] = 3;
                            if (ionized) { nH = 1; angle = 125.5; dihedral = 180; return; }
                            else { nH = 0; return; }
                        } else {
                            nH = 1; angle = 125.5; dihedral = 180; hybrid[icc] = 3; return;
                        }
                    } else {
                        nH = 1; angle = 125.5; dihedral = 180; hybrid[icc] = 3; return;
                    }
                }
                hybrid[icc] = 3;
                if (nat[ii] == 6 && nat[jj] == 6) {
                    if (aromatic(icc, ii, jj)) { nH = 0; hybrid[icc] = 2; return; }
                    if (aromatic_5(icc, ii, jj, hybrid)) { nH = 1; angle = 125.5; dihedral = 180; return; }
                }
                if (nat[ii] == 6 && nat[jj] == 8) {
                    if (nbonds[jj] == 1) nH = 0;
                    else { nH = 1; angle = 125.5; dihedral = 180; }
                    return;
                }
                if (nat[jj] == 15) {
                    if (Rac < 1.7) nH = 0;
                    else { nH = 1; angle = 120; dihedral = 180; }
                    return;
                }
                nH = 1; angle = 120; dihedral = 180;
            } else {
                hybrid[icc] = 3;
            }
            break;
        }

        case 8: {  // Oxygen
            if (nbonds[icc] == 1) {
                if (nat[ii] == 6) {
                    bond_length = 0.975; angle = 110.3; nH = 1; hybrid[icc] = 3;
                    if (nbonds[ii] == 4) { hybrid[icc] = 3; return; }
                    if (nbonds[ii] == 3) {
                        for (i = 1; i <= 3; ++i) { jj = ibonds[i][ii]; if (jj != icc) break; }
                        for (i = i + 1; i <= 3; ++i) { kk = ibonds[i][ii]; if (kk != icc && kk != jj) break; }
                        if (aromatic(ii, jj, kk)) { nH = 1; hybrid[icc] = 3; return; }
                    }
                    jj = nbonds[ii];
                    if (jj == 3) {
                        k = 0; i = 0;
                        for (l = 1; l <= jj; ++l) {
                            kk = ibonds[l][ii];
                            if (nat[kk] == 8) {
                                if (nbonds[kk] == 2) i = i + 1;
                                k = k + 1;
                                if (kk != icc) j = kk;
                            }
                        }
                        if (k == 2) {
                            if (fsub(txtatm[ii], 15, 16) == "  ") {
                                l_tmp = (icc > j);
                            } else {
                                Rac = distance(ii, j);
                                l_tmp = (Rab > Rac && i == 0);
                            }
                            if (l_tmp) {
                                if (ionized) { nH = 0; }
                                else {
                                    bond_length = 1.00; nH = 1; angle = 110; dihedral = 0;
                                    for (l = 1; l <= jj; ++l) {
                                        kk = ibonds[l][ii];
                                        if (nat[kk] == 8 && kk != icc) nc_icc = kk;
                                    }
                                }
                                return;
                            } else {
                                nH = 0; hybrid[icc] = 2; return;
                            }
                        }
                    }
                    if (Rab > 1.37) { nH = 1; bond_length = 0.97; angle = 110; hybrid[icc] = 3; return; }
                    if (Rab < 1.3 && hybrid[ii] == 2) { nH = 0; hybrid[icc] = 2; return; }
                    if (Rab < 1.23) { nH = 0; return; }
                    if (Rab > 1.17) {
                        nH = 0; j = 0;
                        for (i = 1; i <= nbonds[ii]; ++i)
                            if (nat[ibonds[i][ii]] != 1) j = j + 1;
                        if (nat[ii] == 6 && j == 2) {
                            jj = ibonds[1][ii];
                            if (jj == icc) jj = ibonds[2][ii];
                            bangle(coord, icc, ii, jj, sum);
                            if (sum < 2.03) { nH = 1; bond_length = 0.97; angle = 110; hybrid[icc] = 3; return; }
                        }
                    }
                    if (Rab > 1.27) {
                        if (nat[ii] == 6 && nbonds[ii] == 3) {
                            jj = ibonds[1][ii]; kk = ibonds[2][ii]; ll = ibonds[3][ii];
                            i = 0;
                            if (nat[jj] == 6) i = 1;
                            if (nat[kk] == 6) i = i + 1;
                            if (nat[ll] == 6) i = i + 1;
                            if (i == 2) {
                                dihed(coord, ii, jj, kk, ll, sum);
                                if (std::min(2 * funcon_C::pi - sum, sum) < 0.25) { nH = 1; hybrid[icc] = 3; return; }
                            }
                        }
                    } else {
                        nH = 0; hybrid[icc] = 2;
                    }
                } else if (nat[ii] == 7) {
                    if (!ionized) {
                        if (nbonds[ii] == 3) {
                            j = 0;
                            for (i = 1; i <= 3; ++i)
                                if (nat[ibonds[i][ii]] == 8) j = j + 1;
                            if (j == 3) {
                                for (i = 1; i <= 3; ++i) {
                                    k = ibonds[i][ii];
                                    if (nbonds[k] != 2) continue;
                                    if (ibonds[1][k] != ii) {
                                        j = ibonds[1][k];
                                        ibonds[1][k] = ibonds[2][k];
                                        ibonds[2][k] = j;
                                    }
                                }
                                for (ll = 1; ll <= 2; ++ll) {
                                    i = ibonds[ll][ii];
                                    if (nbonds[i] == 2) {
                                        jj = ibonds[2][i];
                                        for (k = ll + 1; k <= 3; ++k) {
                                            j = ibonds[k][ii];
                                            if (nbonds[j] == 2) {
                                                kk = ibonds[2][j];
                                                if (i == kk) { nbonds[i] = 1; nbonds[j] = 1; }
                                            }
                                        }
                                    }
                                }
                                j = 0;
                                for (i = 1; i <= 3; ++i)
                                    if (nbonds[ibonds[i][ii]] > 1) j = 1;
                                if (j == 0) {
                                    nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 3;
                                    for (i = 1; i <= 3; ++i)
                                        if (ibonds[i][ii] != icc) break;
                                    jj = ibonds[i][ii];
                                    for (i = i + 1; i <= 3; ++i)
                                        if (ibonds[i][ii] != icc) break;
                                    jj = ibonds[i][ii];
                                }
                            } else {
                                if (Rab > 1.3) { nH = 1; angle = 120; dihedral = 180; hybrid[icc] = 3; }
                                else nH = 0;
                            }
                        }
                    }
                } else if (nat[ii] == 8) {
                    nH = 1; angle = 104; dihedral = 90;
                } else if (nat[ii] == 12 || nat[ii] == 20) {
                    nH = 1; angle = 150; hybrid[icc] = 3;
                } else if (nat[ii] == 15 || nat[ii] == 33) {
                    hybrid[icc] = 3;
                    if (nbonds[ii] == 3) {
                        nH = 1; angle = 120; dihedral = 180;
                        for (i = 1; i <= 3; ++i) {
                            nd_icc = ibonds[i][ii];
                            if (nd_icc != nb_icc && nd_icc != icc && nd_icc != nc_icc) break;
                        }
                        hybrid[icc] = 3;
                    } else if (nbonds[ii] == 4) {
                        jj = 0; j = 0;
                        for (i = 1; i <= nbonds[ii]; ++i) {
                            kk = ibonds[i][ii];
                            if (nbonds[kk] == 1 && nat[kk] == 8) jj = jj + 1;
                            if (nat[kk] == 8) j = j + 1;
                        }
                        if (ionized && j == 4) { nH = 0; return; }
                        if (jj > 1) { nH = 1; angle = 150; dihedral = 0; }
                    }
                } else if (nat[ii] == 16 || nat[ii] == 34) {
                    hybrid[icc] = 3;
                    if (nbonds[ii] == 2) {
                        if (nat[ibonds[1][ii]] == 6 || nat[ibonds[2][ii]] == 6) { nH = 1; angle = 120; dihedral = 90; }
                    } else if (nbonds[ii] > 2) {
                        jj = 0; j = 0;
                        for (i = 1; i <= nbonds[ii]; ++i) {
                            kk = ibonds[i][ii];
                            if (nbonds[kk] == 1 && nat[kk] == 8) jj = jj + 1;
                            if (nat[kk] == 8) j = j + 1;
                        }
                        if (ionized && j == 4) { nH = 0; return; }
                        if (jj > nbonds[ii] - 2) { nH = 1; angle = 150; dihedral = 0; }
                    } else {
                        hybrid[icc] = 3;
                    }
                } else {
                    hybrid[icc] = 3;
                }
            } else {
                hybrid[icc] = 2;
            }
            break;
        }

        case 14: case 32: case 50: {
            bond_length = 1.3;
            switch (nbonds[icc]) {
                case 1: nH = 3; angle = 109.471; internal_dihedral = 120; dihedral = 180; hybrid[icc] = 3; break;
                case 2: nH = 2; angle = 109.5; internal_dihedral = 180; dihedral = 90; hybrid[icc] = 3; break;
                case 3: nH = 1; angle = 125.5; dihedral = 180; hybrid[icc] = 3; break;
                case 4: case 5: case 6: nH = 0; break;
            }
            return;
        }
        case 15: {
            hybrid[icc] = 3;
            if (nbonds[icc] == 1) { nH = 2; bond_length = 1.42; angle = 93.5; internal_dihedral = 93.8; }
            else if (nbonds[icc] == 2) { nH = 1; bond_length = 1.42; angle = 122; dihedral = 180; }
            break;
        }
        case 16: {
            hybrid[icc] = 3;
            if (nbonds[icc] == 1) {
                if ((Rab < 2.0 && nat[ii] == 16) || (Rab < 1.69 && nat[ii] == 6)) {
                    nH = 0;
                } else {
                    if (Rab < 1.65) nH = 0;
                    else { nH = 1; bond_length = 1.3; angle = 100; dihedral = 90; }
                }
            }
            break;
        }
        case 33: {
            hybrid[icc] = 3;
            if (nbonds[icc] == 1) { nH = 2; bond_length = 1.52; angle = 92; internal_dihedral = 92; }
            else if (nbonds[icc] == 2) { nH = 1; bond_length = 1.52; angle = 95; internal_dihedral = 95; }
            break;
        }
        case 34: {
            hybrid[icc] = 3;
            if (nbonds[icc] == 1) { nH = 1; bond_length = 1.46; angle = 91; }
            break;
        }
        case 51: {
            hybrid[icc] = 3;
            if (nbonds[icc] == 1) { nH = 2; bond_length = 1.71; angle = 92; internal_dihedral = 92; }
            else if (nbonds[icc] == 2) { nH = 1; bond_length = 1.71; angle = 95; internal_dihedral = 95; hybrid[icc] = 3; }
            break;
        }
        case 52: {
            hybrid[icc] = 3;
            if (nbonds[icc] == 1) { nH = 1; bond_length = 1.69; angle = 90; }
            break;
        }
        default:
            hybrid[icc] = 3;
            break;
    }
}

// ===========================================================================
//  add_a_generic_hydrogen_atom
// ===========================================================================
void add_a_generic_hydrogen_atom(int na, int nb, int nc, double bond_length,
                                 double angle, double dihedral,
                                 const std::vector<int>& metals, int nmetals) {
    double xb, yb, zb, rbc, xa, ya, za, xpa, xpb, costh, sinth,
        sinph, cosph, zqa, yza, xyb, ypa, coskh, sinkh, sina, cosd,
        xd, yd, zd, xpd, ypd, zqd, xqd, yqd, xrd, cosa, sind, zpd;
    int k;
    cosa = std::cos(angle);
    xb = coord[0][nb] - coord[0][na];
    yb = coord[1][nb] - coord[1][na];
    zb = coord[2][nb] - coord[2][na];
    rbc = xb * xb + yb * yb + zb * zb;
    rbc = 1.0 / std::sqrt(rbc);
    if (nc != 0) {
        xa = coord[0][nc] - coord[0][na];
        ya = coord[1][nc] - coord[1][na];
        za = coord[2][nc] - coord[2][na];
    } else {
        xa = 1.0; ya = 2.0; za = 3.0;
    }
    xyb = std::sqrt(xb * xb + yb * yb);
    k = -1;
    if (xyb <= 0.009) {
        xpa = za; za = -xa; xa = xpa;
        xpb = zb; zb = -xb; xb = xpb;
        xyb = std::sqrt(xb * xb + yb * yb);
        k = 1;
    }
    costh = xb / xyb;
    sinth = yb / xyb;
    xpa = xa * costh + ya * sinth;
    ypa = ya * costh - xa * sinth;
    sinph = zb * rbc;
    cosph = std::sqrt(std::fabs(1.0 - sinph * sinph));
    zqa = za * cosph - xpa * sinph;
    yza = std::sqrt(ypa * ypa + zqa * zqa);
    if (yza >= 1.0e-4) {
        coskh = ypa / yza;
        sinkh = zqa / yza;
    } else {
        coskh = 1.0; sinkh = 0.0;
    }
    sina = std::sin(angle);
    sind = -std::sin(dihedral);
    cosd = std::cos(dihedral);
    xd = bond_length * cosa;
    yd = bond_length * sina * cosd;
    zd = bond_length * sina * sind;
    ypd = yd * coskh - zd * sinkh;
    zpd = zd * coskh + yd * sinkh;
    xpd = xd * cosph - zpd * sinph;
    zqd = zpd * cosph + xd * sinph;
    xqd = xpd * costh - ypd * sinth;
    yqd = ypd * costh + xpd * sinth;
    if (k >= 1) {
        xrd = -zqd; zqd = xqd; xqd = xrd;
    }
    mk::numat = mk::numat + 1;
    coord[0][mk::numat] = xqd + coord[0][na];
    coord[1][mk::numat] = yqd + coord[1][na];
    coord[2][mk::numat] = zqd + coord[2][na];
    nbonds[na] = nbonds[na] + 1;
    nbonds[mk::numat] = 1;
    nat[mk::numat] = 1;
    ibonds[1][mk::numat] = na;
    ibonds[nbonds[na]][na] = mk::numat;
    if (near_a_metal(na, metals, nmetals)) mk::numat = mk::numat - 1;
}

// ===========================================================================
//  add_a_sp3_hydrogen_atom
// ===========================================================================
void add_a_sp3_hydrogen_atom(int icc, int nb_icc, int nc_icc, int nd_icc,
                             double bond_length,
                             const std::vector<int>& metals, int nmetals) {
    double x, y, z, ax, ay, az, r1, r2, sum, angle, dihedral;
    ax = coord[0][icc]; ay = coord[1][icc]; az = coord[2][icc];
    x = (coord[0][nb_icc] + coord[0][nc_icc] + coord[0][nd_icc]) / 3.0;
    y = (coord[1][nb_icc] + coord[1][nc_icc] + coord[1][nd_icc]) / 3.0;
    z = (coord[2][nb_icc] + coord[2][nc_icc] + coord[2][nd_icc]) / 3.0;
    r1 = std::sqrt((x - ax) * (x - ax) + (y - ay) * (y - ay) + (z - az) * (z - az));
    if (r1 < 0.4) {
        angle = 1.570796;
        dihedral = 1.57;
        add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals);
        r1 = (coord[0][nb_icc] - coord[0][mk::numat]) * (coord[0][nb_icc] - coord[0][mk::numat]) +
             (coord[1][nb_icc] - coord[1][mk::numat]) * (coord[1][nb_icc] - coord[1][mk::numat]) +
             (coord[2][nb_icc] - coord[2][mk::numat]) * (coord[2][nb_icc] - coord[2][mk::numat]);
        r1 = std::min(r1, (coord[0][nc_icc] - coord[0][mk::numat]) * (coord[0][nc_icc] - coord[0][mk::numat]) +
             (coord[1][nc_icc] - coord[1][mk::numat]) * (coord[1][nc_icc] - coord[1][mk::numat]) +
             (coord[2][nc_icc] - coord[2][mk::numat]) * (coord[2][nc_icc] - coord[2][mk::numat]));
        r1 = std::min(r1, (coord[0][nd_icc] - coord[0][mk::numat]) * (coord[0][nd_icc] - coord[0][mk::numat]) +
             (coord[1][nd_icc] - coord[1][mk::numat]) * (coord[1][nd_icc] - coord[1][mk::numat]) +
             (coord[2][nd_icc] - coord[2][mk::numat]) * (coord[2][nd_icc] - coord[2][mk::numat]));
        nbonds[icc] = nbonds[icc] - 1;
        dihedral = -dihedral;
        add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals);
        r2 = (coord[0][nb_icc] - coord[0][mk::numat]) * (coord[0][nb_icc] - coord[0][mk::numat]) +
             (coord[1][nb_icc] - coord[1][mk::numat]) * (coord[1][nb_icc] - coord[1][mk::numat]) +
             (coord[2][nb_icc] - coord[2][mk::numat]) * (coord[2][nb_icc] - coord[2][mk::numat]);
        r2 = std::min(r2, (coord[0][nc_icc] - coord[0][mk::numat]) * (coord[0][nc_icc] - coord[0][mk::numat]) +
             (coord[1][nc_icc] - coord[1][mk::numat]) * (coord[1][nc_icc] - coord[1][mk::numat]) +
             (coord[2][nc_icc] - coord[2][mk::numat]) * (coord[2][nc_icc] - coord[2][mk::numat]));
        r2 = std::min(r2, (coord[0][nd_icc] - coord[0][mk::numat]) * (coord[0][nd_icc] - coord[0][mk::numat]) +
             (coord[1][nd_icc] - coord[1][mk::numat]) * (coord[1][nd_icc] - coord[1][mk::numat]) +
             (coord[2][nd_icc] - coord[2][mk::numat]) * (coord[2][nd_icc] - coord[2][mk::numat]));
        if (r2 > r1) {
            coord[0][mk::numat - 1] = coord[0][mk::numat];
            coord[1][mk::numat - 1] = coord[1][mk::numat];
            coord[2][mk::numat - 1] = coord[2][mk::numat];
        }
        mk::numat = mk::numat - 1;
        return;
    }
    sum = bond_length / r1;
    mk::numat = mk::numat + 1;
    coord[0][mk::numat] = ax + sum * (ax - x);
    coord[1][mk::numat] = ay + sum * (ay - y);
    coord[2][mk::numat] = az + sum * (az - z);
    nbonds[icc] = nbonds[icc] + 1;
    nbonds[mk::numat] = 1;
    nat[mk::numat] = 1;
    ibonds[1][mk::numat] = icc;
    ibonds[nbonds[icc]][icc] = mk::numat;
    if (near_a_metal(icc, metals, nmetals)) mk::numat = mk::numat - 1;
}

// ===========================================================================
//  aromatic
// ===========================================================================
bool aromatic(int atom_1, int atom_2, int atom_6) {
    int atom_3, atom_5, i, j, k, l, m, nii, njj, kka, atom_4;
    const double bent = 0.4;
    double sum;
    bool result = false;
    nii = nbonds[atom_2];
    njj = nbonds[atom_6];
    for (i = 1; i <= nii; ++i) {
        atom_3 = ibonds[i][atom_2];
        if (atom_3 == atom_1) continue;
        for (j = 1; j <= njj; ++j) {
            atom_5 = ibonds[j][atom_6];
            if (atom_5 == atom_1) continue;
            for (k = 1; k <= nbonds[atom_5]; ++k) {
                kka = ibonds[k][atom_5];
                if (kka == atom_6) continue;
                for (l = 1; l <= nbonds[atom_3]; ++l) {
                    atom_4 = ibonds[l][atom_3];
                    if (atom_4 == atom_2) continue;
                    if (atom_4 == atom_5) {
                        if (nat[atom_4] == 8 || nat[atom_4] == 16 || nat[atom_3] == 8 || nat[atom_3] == 16) {
                            dihed(coord, atom_1, atom_2, atom_4, atom_6, sum);
                            if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                dihed(coord, atom_1, atom_2, atom_3, atom_4, sum);
                                result = (std::min(2 * funcon_C::pi - sum, sum) < bent);
                            }
                        }
                        return result;
                    }
                    if (kka == atom_4) {
                        if (atom_1 == atom_3 || atom_1 == atom_4 || atom_1 == atom_5) return result;
                        if (atom_2 == atom_3 || atom_2 == atom_4 || atom_2 == atom_5) return result;
                        if (atom_3 == atom_4 || atom_3 == atom_5 || atom_3 == atom_6) return result;
                        if (atom_4 == atom_5 || atom_4 == atom_6 || atom_5 == atom_6) return result;
                        m = std::max({nat[atom_1], nat[atom_2], nat[atom_3], nat[atom_4], nat[atom_5], nat[atom_6]});
                        if (m == 8 || m == 14 || m == 16) return result;
                        m = 0;
                        if (nat[atom_1] == 7 && nbonds[atom_1] > 2) m = m + 1;
                        if (nat[atom_2] == 7 && nbonds[atom_2] > 2) m = m + 1;
                        if (nat[atom_3] == 7 && nbonds[atom_3] > 2) m = m + 1;
                        if (nat[atom_4] == 7 && nbonds[atom_4] > 2) m = m + 1;
                        if (nat[atom_5] == 7 && nbonds[atom_5] > 2) m = m + 1;
                        if (nat[atom_6] == 7 && nbonds[atom_6] > 2) m = m + 1;
                        if (m > 1) return result;
                        dihed(coord, atom_1, atom_2, atom_3, atom_4, sum);
                        if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                            dihed(coord, atom_2, atom_3, atom_4, atom_5, sum);
                            if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                dihed(coord, atom_3, atom_4, atom_5, atom_6, sum);
                                if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                    dihed(coord, atom_4, atom_5, atom_6, atom_1, sum);
                                    if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                        dihed(coord, atom_5, atom_6, atom_1, atom_2, sum);
                                        if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                            dihed(coord, atom_6, atom_1, atom_2, atom_3, sum);
                                            if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                                result = true;
                                                return result;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

// ===========================================================================
//  aromatic_5
// ===========================================================================
bool aromatic_5(int atom_1, int atom_2, int atom_5, std::vector<int>& hybrid) {
    int atom_3, i, j, k, nii, njj, kka, atom_4, l;
    const double bent = 0.4;
    double sum;
    if (mk::numat < 5) return false;
    nii = nbonds[atom_2];
    njj = nbonds[atom_5];
    for (i = 1; i <= nii; ++i) {
        atom_3 = ibonds[i][atom_2];
        if (atom_3 == atom_1) continue;
        for (j = 1; j <= njj; ++j) {
            atom_4 = ibonds[j][atom_5];
            if (atom_4 == atom_1) continue;
            for (k = 1; k <= nbonds[atom_4]; ++k) {
                kka = ibonds[k][atom_4];
                if (kka == atom_3) {
                    l = std::max({hybrid[2], hybrid[3], hybrid[4], hybrid[5]});
                    if (l == 3) continue;
                    dihed(coord, atom_1, atom_2, atom_3, atom_4, sum);
                    if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                        dihed(coord, atom_2, atom_3, atom_4, atom_5, sum);
                        if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                            dihed(coord, atom_3, atom_4, atom_5, atom_1, sum);
                            if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                dihed(coord, atom_4, atom_5, atom_1, atom_2, sum);
                                if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                    dihed(coord, atom_5, atom_1, atom_2, atom_3, sum);
                                    if (std::min(2 * funcon_C::pi - sum, sum) < bent) {
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

// ===========================================================================
//  near_a_metal
// ===========================================================================
bool near_a_metal(int icc, const std::vector<int>& metals, int nmetals) {
    int i;
    bool Sulfur, Oxygen;
    double Rab, Rac;
    for (i = 1; i <= nmetals; ++i) {
        Rab = distance(metals[i], icc);
        Rac = distance(metals[i], mk::numat);
        Sulfur = (nat[icc] == 16 && Rab < 2.6);
        Oxygen = (nat[icc] == 8 && Rab < 2.1);
        if (Rac < 1.5) break;
        if ((Sulfur || Oxygen) && nbonds[icc] > 0) break;
    }
    bool result = (i <= nmetals);
    if (!result) return result;
    nbonds[icc] = nbonds[icc] - 1;
    return result;
}

}  // namespace (h_type etc.)

// ===========================================================================
//  bridge_H
// ===========================================================================
void bridge_H() {
    int i, ii = 0, j, k = 0, l = 0, kk, jj, nchain1, nchain2, mbreaks;
    double bond_length, bond_length_2, sum;
    std::vector<int> h1, h2;
    std::vector<char> pchain;
    std::vector<double> Rab1, Rab2;

    h1.assign(mk::numat + 1, 0);
    h2.assign(mk::numat + 1, 0);
    Rab1.assign(mk::numat + 1, 100.0);
    Rab2.assign(mk::numat + 1, 100.0);
    pchain.assign(mk::numat + 1, ' ');
    ii = 0;
    nchain1 = 1;
    mbreaks = 1;
    for (i = 1; i <= mk::numat; ++i) {
        if (i == breaks[mbreaks]) {
            mbreaks = mbreaks + 1;
            nchain1 = nchain1 + 1;
        }
        pchain[i] = chains[nchain1 - 1];
    }
    for (i = 1; i <= mk::numat; ++i) {
        if (nat[i] != 1) continue;
        bond_length = 10.0;
        bond_length_2 = 10.0;
        nchain2 = 1;
        l = 0;
        for (j = 1; j <= mk::numat; ++j) {
            if (j == breaks[nchain2]) nchain2 = nchain2 + 1;
            if (nat[j] == 1) continue;
            sum = distance(i, j);
            jj = 0;
            if (nat[j] == 8) {
                jj = 0;
                for (kk = 1; kk <= nbonds[j]; ++kk)
                    if (nat[ibonds[kk][j]] == 1) jj = jj + 1;
            }
            if (sum < bond_length && jj < 2) {
                if (bond_length < bond_length_2) {
                    l = k;
                    bond_length_2 = bond_length;
                }
                k = j;
                bond_length = sum;
            }
        }
        Rab1[i] = bond_length;
        h1[i] = k;
        Rab2[i] = bond_length_2;
        h2[i] = l;
    }
    for (j = 1; j <= 20; ++j) {
        bond_length_2 = 100.0;
        for (i = 1; i <= mk::numat; ++i) {
            if (Rab2[i] < bond_length_2) {
                bond_length_2 = Rab2[i];
                k = i;
            }
        }
        if (bond_length_2 > 1.5) break;
        if (ii == 0) {
            ii = 1;
            mk::line = "  UNUSUALLY SHORT HYDROGEN BOND DISTANCES";
            // (output to iw / ilog elided in C++ port)
        }
        // formatted line construction elided; uses elemnt, txtatm, pchain.
        Rab2[k] = 10.0;
    }
}

// ===========================================================================
//  reset_breaks
// ===========================================================================
void reset_breaks() {
    int i, j, k, l;
    double sum;
    bool first = true;
    for (mk::nbreaks = 1; mk::nbreaks <= 50; ++mk::nbreaks) {
        bool found = false;
        for (j = 1; j <= mk::numat; ++j) {
            sum = std::fabs(coord[0][j] - break_coords[0][mk::nbreaks]) +
                  std::fabs(coord[1][j] - break_coords[1][mk::nbreaks]) +
                  std::fabs(coord[2][j] - break_coords[2][mk::nbreaks]);
            if (sum < 0.1) { found = true; break; }
        }
        if (!found) break;
        mk::line = fsub(txtatm[j], 23, (int)txtatm[j].size());
        k = (int)std::lround(reada(mk::line, 1));
        while (true) {
            j = j + 1;
            if (j > mk::numat) break;
            mk::line = fsub(txtatm[j], 23, (int)txtatm[j].size());
            l = (int)std::lround(reada(mk::line, 1));
            if (l != k) break;
        }
        breaks[mk::nbreaks] = j - 1;
    }
    j = 1;
    mk::nbreaks = mk::nbreaks - 1;
    for (i = 1; i <= mk::nbreaks; ++i) {
        l = breaks[i];
        for (k = i + 1; k <= mk::nbreaks; ++k) {
            if (breaks[k] < l) { breaks[i] = breaks[k]; breaks[k] = l; l = breaks[i]; }
        }
        if (i > 1) {
            if (breaks[i] != breaks[i - 1]) { j = j + 1; breaks[j] = breaks[i]; }
        }
    }
    mk::nbreaks = j;
    breaks[mk::nbreaks + 1] = 0;

    if (mk::keywrd.find(" NORES") == std::string::npos) {
        for (i = 1; i <= mk::nbreaks; ++i) {
            for (j = i + 1; j <= mk::nbreaks; ++j) {
                if (breaks[j] < breaks[i]) {
                    if (first && mk::keywrd.find("GEO-OK") == std::string::npos) {
                        first = false;
                        mopend("WARNING: After RESEQ, PDB \"TER\" locations are incorrect");
                    }
                }
            }
        }
    }
    j = 1;
    for (i = 2; i <= mk::nbreaks; ++i) {
        if (breaks[i] - breaks[i - 1] > 0) { j = j + 1; breaks[j] = breaks[i]; }
    }
    mk::nbreaks = j;
    breaks[mk::nbreaks + 1] = 0;
    j = 1;
    for (i = 1; i <= mk::numat; ++i) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%-6s%5d%15s",
                      txtatm[i].substr(0, 6).c_str(), i + j - 1,
                      txtatm[i].substr(11).c_str());
        txtatm[i] = buf;
        if (i == breaks[j]) j = j + 1;
    }
}

// Public wrapper around the file-local add_a_sp3_hydrogen_atom, so that
// geochk.F90 ("site") can call it from outside this translation unit.
void add_a_sp3_hydrogen_atom_ext(int icc, int nb_icc, int nc_icc, int nd_icc,
                                 double bond_length,
                                 const std::vector<int>& metals, int nmetals) {
    add_a_sp3_hydrogen_atom(icc, nb_icc, nc_icc, nd_icc, bond_length,
                            metals, nmetals);
}

