// hcore.cpp — C++ translation of MOPAC 2016 "hcore.F90".
// Generates the one-electron matrix H and the two-electron integrals W for a
// molecule in Cartesian coordinates. Also accumulates the nuclear energy
// ENUCLR, handles an applied electric field (FIELD keyword) and the QM/MM
// embedding charges (mol.in).
#include "hcore.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "funcon_C.h"
#include "parameters_C.h"
#include "MOZYME_C.h"
#include "overlaps_C.h"
#include "chanel_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
using namespace MOZYME_C;

namespace molkst_C {
extern int numcal, numat, norbs, id, l1u, l2u, l3u, n2elec, mpack;
extern std::string keywrd, line;
extern double enuclr, efield[4];
}
namespace common_arrays_C {
extern std::vector<int> nfirst, nlast, nat;
extern std::vector<double> uspd, h, w, wk;
extern std::vector<std::vector<double>> coord;
}
namespace cosmo_C { extern bool useps; }
namespace funcon_C { extern double a0, ev, fpc_9; }
namespace parameters_C { extern double tore[], dd[]; }
namespace MOZYME_C { extern double cutofs; }
namespace overlaps_C { extern double cutof1, cutof2; }
namespace chanel_C { extern int iw; }

extern double reada(const std::string&, int);
namespace mndod_C {
extern std::vector<int> intij, intkl, intrep;
extern std::vector<std::vector<double>> repd;
}
extern void add_path(std::string&);
extern void mopend(const std::string&);
extern void h1elec(int, int, const double*, const double*, double*);
extern void rotate(int, int, const double*, const double*, double*, int&, double*, double*, double&);
extern void solrot(int, int, const double*, const double*, double*, double*, int&, double*, double*, double&);
extern void addhcr();
extern void addnuc();
extern void vecprt(double*, int);

void hcore() {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace funcon_C;
    using namespace parameters_C;
    using namespace overlaps_C;
    using namespace cosmo_C;

    static bool first = true, debug = false, lmolaris_qmmm = false;
    static int icalcn = 0, ione = 1;
    static std::vector<double> vqc;

    first = (icalcn != numcal);
    icalcn = numcal;
    if (first) {
        ione = 1;
        MOZYME_C::cutofs = 1.0e4;
        cutof2 = 1.0e10;
        cutof1 = 225.0;  // (Cutoff distance)^2 for overlap integrals
        if (id != 0) {
            wk.clear();
            wk.resize((size_t)n2elec + 1, 0.0);
        }
        if (id != 0) ione = 0;
        debug = (keywrd.find("HCORE") != std::string::npos);
        lmolaris_qmmm = (keywrd.find("QMMM") != std::string::npos);
        if (lmolaris_qmmm) {
            // Read in the energy (kcal/mol) that an electron on each atom would
            // have, arising from the partial charges on all atoms.
            line = "mol.in";
            add_path(line);
            std::ifstream f(line.c_str());
            if (!f.good()) {
                mopend("A file named '" + line + "' was expected, but was not found.");
                return;
            }
            int q_at = 0, link_at = 0;
            if (!(f >> q_at >> link_at)) {
                if (f.eof()) {
                    mopend("A file named '" + line + "' exists, but is empty.");
                    return;
                }
                mopend("A file named '" + line + "' exists, but is corrupt.");
                return;
            }
            int mol_at = q_at + link_at;
            if (mol_at != numat) {
                char buf[160];
                std::snprintf(buf, sizeof(buf),
                              " The number of atoms in 'mol.in':%5d and in the MOPAC data set:%5d are different.",
                              mol_at, numat);
                mopend(buf);
                std::fprintf(stdout, "          Correct fault and resubmit.\n");
                return;
            }
            vqc.clear();
            vqc.resize((size_t)numat + 1, 0.0);
            for (int i = 1; i <= numat; ++i) {
                std::string s1, s2, s3, s4;
                double v;
                if (!(f >> s1 >> s2 >> s3 >> s4 >> v)) break;
                vqc[i] = v;
                if (debug) std::fprintf(stdout, "ATOM No.%4d VQC(I)%9.3f\n", i, vqc[i]);
            }
        }
        double xf = 0.0, yf = 0.0, zf = 0.0;
        std::string tmpkey = keywrd;
        int i = (int)tmpkey.find(" FIELD(") + 1;
        size_t pos_eq = tmpkey.find(" FIELD=(");
        if (pos_eq != std::string::npos) i = (int)pos_eq + 1;
        if (i != 0) {
            // Erase all text from TMPKEY except FIELD data.
            std::fill(tmpkey.begin(), tmpkey.begin() + i, ' ');
            size_t cp = tmpkey.find(')');
            if (cp != std::string::npos) {
                std::fill(tmpkey.begin() + cp, tmpkey.end(), ' ');
                tmpkey.resize(cp + 1);
            }
            xf = reada(tmpkey, i);
            size_t cm = tmpkey.find(',');
            if (cm != std::string::npos) {
                tmpkey[cm] = ' ';
                yf = reada(tmpkey, (int)cm + 1);
                cm = tmpkey.find(',');
                if (cm != std::string::npos) {
                    tmpkey[cm] = ' ';
                    zf = reada(tmpkey, (int)cm + 1);
                }
            }
            std::fprintf(stdout, "\n          THE ELECTRIC FIELD IS%10.5f%10.5f%10.5f VOLTS/ANGSTROM\n\n",
                         xf, yf, zf);
        }
        double constv = a0 / ev;
        efield[1] = xf * constv;
        efield[2] = yf * constv;
        efield[3] = zf * constv;
    }
    bool fldon = false;
    if (efield[1] != 0.0 || efield[2] != 0.0 || efield[3] != 0.0) {
        fldon = true;
    }
    enuclr = 0.0;
    std::fill(h.begin(), h.begin() + mpack, 0.0);
    int kr = 1;
    for (int i = 1; i <= numat; ++i) {
        int ia = nfirst[i], ib = nlast[i], ni = nat[i];
        // First fill the diagonals and off-diagonals on the same atom.
        if (!fldon) {
            for (int i1 = ia; i1 <= ib; ++i1) {
                int i2 = i1 * (i1 - 1) / 2 + ia - 1;
                if (i1 - ia + 1 > 0) {
                    for (int k = 1; k <= i1 - ia + 1; ++k) h[i2 + k] = 0.0;
                    i2 = i1 - ia + 1 + i2;
                }
                h[i2] = uspd[i1];
                if (lmolaris_qmmm) {
                    if (debug) std::fprintf(stdout, "OLD 1e MATRIX ELEMENT %5d%12.5f i%5d vqc(i)%12.5f\n",
                                            i2, h[i2], i, -vqc[i] / fpc_9);
                    h[i2] = h[i2] - vqc[i] / fpc_9;
                    if (debug) std::fprintf(stdout, "UPD 1e MATRIX ELEMENT %5d%12.5f\n", i2, h[i2]);
                }
            }
        } else {
            double fnuc = 0.0;
            double fldcon = ev / a0;
            for (int i1 = ia; i1 <= ib; ++i1) {
                int i2 = i1 * (i1 - 1) / 2 + ia - 1;
                for (int j1 = ia; j1 <= i1; ++j1) {
                    ++i2;
                    h[i2] = 0.0;
                    int io1 = i1 - ia, jo1 = j1 - ia;
                    if (jo1 == 0 && io1 == 1) h[i2] = -a0 * dd[ni] * efield[1] * fldcon;
                    if (jo1 == 0 && io1 == 2) h[i2] = -a0 * dd[ni] * efield[2] * fldcon;
                    if (jo1 != 0 || io1 != 3) continue;
                    h[i2] = -a0 * dd[ni] * efield[3] * fldcon;
                }
                h[i2] = uspd[i1];
                fnuc = -(efield[1] * coord[0][i] + efield[2] * coord[1][i] + efield[3] * coord[2][i]) * fldcon;
                h[i2] = h[i2] + fnuc;
            }
            enuclr = enuclr - fnuc * tore[nat[i]];
        }
        if (lmolaris_qmmm) enuclr = enuclr + vqc[i] / fpc_9 * tore[nat[i]];
        // Fill the atom-other-atom one-electron matrix <PSI(lambda)|PSI(sigma)>.
        int im1 = i - ione;
        for (int j = 1; j <= im1; ++j) {
            double half = 1.0;
            if (i == j) half = 0.5;
            int ja = nfirst[j], jb = nlast[j], nj = nat[j];
            double di[100];  // 1-based semantics: h1elec writes flat[11..91] (9x9 col-major), flat[0..10] padding
            if (id == 0) {
                            double xi3[3] = {coord[0][i], coord[1][i], coord[2][i]};
                double xj3[3] = {coord[0][j], coord[1][j], coord[2][j]};
                h1elec(ni, nj, xi3, xj3, &di[11]);
            } else {
                for (int kk = 0; kk < 100; ++kk) di[kk] = 0.0;
                double dibits[100];
                double xj[4];
                for (int ii = -l1u; ii <= l1u; ++ii)
                    for (int jj = -l2u; jj <= l2u; ++jj)
                        for (int k = -l3u; k <= l3u; ++k) {
                            for (int a = 1; a <= 3; ++a)
                                xj[a] = coord[a-1][j] + common_arrays_C::tvec[a][1] * ii + common_arrays_C::tvec[a][2] * jj + common_arrays_C::tvec[a][3] * k;
                            double xi3b[3] = {coord[0][i], coord[1][i], coord[2][i]};
                            h1elec(ni, nj, xi3b, &xj[1], &dibits[11]);
                            for (int kk = 0; kk < 81; ++kk) di[11 + kk] += dibits[11 + kk];
                        }
            }
            int i2 = 0;
            for (int i1 = ia; i1 <= ib; ++i1) {
                int ii = i1 * (i1 - 1) / 2 + ja - 1;
                ++i2;
                int jj = std::min(i1, jb);
                for (int k = 1; k <= jj - ja + 1; ++k)
                    h[ii + k] = h[ii + k] + di[11 + (k - 1) * 9 + (i2 - 1)];  // 1-based (i2,k) col-major, 9-stride
            }
            // Calculate the two-electron integrals W, the electron-nuclear
            // terms E1B and E2A, and the nuclear-nuclear term ENUC.
            double e1b[45], e2a[45], enuc = 0.0;
            if (id == 0) {
                            double xi3c[3] = {coord[0][i], coord[1][i], coord[2][i]};
                double xj3c[3] = {coord[0][j], coord[1][j], coord[2][j]};
                rotate(ni, nj, xi3c, xj3c, &w[kr], kr, e1b, e2a, enuc);
                        } else {
                int kro = kr;
                double wjd[2025], wkd[2025];
                solrot(ni, nj, &coord[1][i], &coord[1][j], wjd, wkd, kr, e1b, e2a, enuc);
                for (int k = 0; k < kr - kro; ++k) {
                    w[kro + k] = wjd[k];
                    wk[kro + k] = wkd[k];
                }
            }
            enuclr = enuclr + enuc;  // F90 hcore:248
            // Add on the electron-nuclear attraction term for atom I.
            i2 = 0;
            for (int i1 = ia; i1 <= ib; ++i1) {
                int ii = i1 * (i1 - 1) / 2 + ia - 1;
                if (i1 - ia + 1 > 0) {
                    for (int k = 1; k <= i1 - ia + 1; ++k)
                        h[ii + k] = h[ii + k] + e1b[i2 + k - 1] * half;  // e1b is 0-based packed; F90 reads e1b(i2+1:)
                    i2 = i1 - ia + 1 + i2;
                }
            }
            // Add on the electron-nuclear attraction term for atom J.
            i2 = 0;
            for (int i1 = ja; i1 <= jb; ++i1) {
                int ii = i1 * (i1 - 1) / 2 + ja - 1;
                if (i1 - ja + 1 > 0) {
                    for (int k = 1; k <= i1 - ja + 1; ++k)
                        h[ii + k] = h[ii + k] + e2a[i2 + k - 1] * half;  // e2a is 0-based packed; F90 reads e2a(i2+1:)
                    i2 = i1 - ja + 1 + i2;
                }
            }
        }
        int ii = ib - ia + 1;
        ii = (ii * (ii + 1)) / 2;
        if (id != 0) {
            for (int i1 = kr; i1 <= kr + ii * ii - 1; ++i1) wk[i1] = 0.0;
        }
        // F90 hcore.F90:279-281: one-center block stored at end of each atom's
        // loop; fock2's kk advances per atom (pairs then fock1dorbs), aligned.
        if (ii != 0) wstore(&w[kr], kr, ni, ii);
    }

    if (useps) {
        // Dielectric correction to the core-core interaction is added to ENUCLR.
        addnuc();
        // Dielectric correction for the electron interaction is added to H.
        addhcr();
    }
    if (debug) {
        std::fprintf(stdout, "\n\n          ONE-ELECTRON MATRIX FROM HCORE\n");
        vecprt(&h[1], norbs);
        int j = kr - 1;
        if (id == 0) {
            std::fprintf(stdout, "\n\n          TWO-ELECTRON MATRIX IN HCORE\n");
        } else {
            std::fprintf(stdout, "\n\n          TWO-ELECTRON J MATRIX IN HCORE\n");
        }
        for (int i = 1; i <= j; ++i) std::fprintf(stdout, "%8.4f", w[i]);
        std::fprintf(stdout, "\n");
        if (id != 0) {
            std::fprintf(stdout, "\n\n          TWO-ELECTRON K MATRIX IN HCORE\n");
            for (int i = 1; i <= j; ++i) std::fprintf(stdout, "%8.4f", wk[i]);
            std::fprintf(stdout, "\n");
        }
    }
}

void wstore(double* w, int& kr, int ni, int ilim) {
    // F90 mndod.F90:2055: fill the one-center two-electron block and
    // unconditionally advance kr by ilim^2.
    // Official Fortran w(ij,kl) is column-major: linear = (kl-1)*ilim + ij.
    auto W=[&](int ij,int kl)->double& { return w[(kl-1)*ilim + (ij-1)]; };
    for (int i = 0; i < ilim * ilim; ++i) w[i] = 0.0;
    int ip = 1;
    W(ip,ip) = parameters_C::gss[ni];
    if (parameters_C::natorb[ni] > 2) {
        int ipx = ip + 2, ipy = ip + 5, ipz = ip + 9;
        W(ipx,ip)=parameters_C::gsp[ni]; W(ipy,ip)=parameters_C::gsp[ni]; W(ipz,ip)=parameters_C::gsp[ni];
        W(ip,ipx)=parameters_C::gsp[ni]; W(ip,ipy)=parameters_C::gsp[ni]; W(ip,ipz)=parameters_C::gsp[ni];
        W(ipx,ipx)=parameters_C::gpp[ni]; W(ipy,ipy)=parameters_C::gpp[ni]; W(ipz,ipz)=parameters_C::gpp[ni];
        W(ipy,ipx)=parameters_C::gp2[ni]; W(ipz,ipx)=parameters_C::gp2[ni]; W(ipz,ipy)=parameters_C::gp2[ni];
        W(ipx,ipy)=parameters_C::gp2[ni]; W(ipx,ipz)=parameters_C::gp2[ni]; W(ipy,ipz)=parameters_C::gp2[ni];
        W(ip+1,ip+1)=parameters_C::hsp[ni]; W(ip+3,ip+3)=parameters_C::hsp[ni]; W(ip+6,ip+6)=parameters_C::hsp[ni];
        W(ip+4,ip+4)=0.5*(parameters_C::gpp[ni]-parameters_C::gp2[ni]);
        W(ip+7,ip+7)=0.5*(parameters_C::gpp[ni]-parameters_C::gp2[ni]);
        W(ip+8,ip+8)=0.5*(parameters_C::gpp[ni]-parameters_C::gp2[ni]);
        if (ilim > 10) {
            int ij0 = ip - 1;
            for (int i = 1; i <= 243; ++i) {
                int ij = mndod_C::intij[i], kl = mndod_C::intkl[i], Int = mndod_C::intrep[i];
                W(ij+ij0, kl+ij0) = mndod_C::repd[Int][ni];
            }
        }
    }
    kr += ilim * ilim;
}
