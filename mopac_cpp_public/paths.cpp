// paths.cpp — C++ translation of "paths.F90".
// Follow a user-supplied reaction coordinate (PATHS): at each step, set the
// reaction-coordinate value and optimize with EF (default) or FLEPO/DFP.

#include "paths.h"

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "ef_C.h"
#include "elemts_C.h"
#include "maps_C.h"
#include "molkst_C.h"

// External routines (ported elsewhere or stubbed).
void ef(std::vector<double>& xparam, double& escf);
void flepo(std::vector<double>& xparam, int nvar, double& escf);
void writmo();
void pdbout(int unit);
void dfpsav(double& totime, std::vector<double>& xparam, std::vector<double>& gd,
            std::vector<double>& xlast, double& funct1, std::vector<int>& mdfp,
            std::vector<double>& xdfp);
double reada(const std::string& s, int pos);
double second(int);
void mopend(const std::string& s);
void to_screen(const std::string& s);

void paths() {
    using namespace common_arrays_C;
    using namespace maps_C;
    using namespace ef_C;

    std::vector<int> mdfp(21, 0);
    std::vector<double> xdfp(21, 0.0);
    int maxcyc = 100000;
    int percent = 0, imodel = 0;

    iloop = 1;
    int iw00 = chanel_C::iw0;
    bool debug = molkst_C::keywrd.find("DEBUG") != std::string::npos;

    std::ofstream ixyz;
    if (molkst_C::keywrd.find("PDBOUT") != std::string::npos) {
        // open(ipdb, xyz_fn minus ".xyz" + ".pdb")
        std::string fn = chanel_C::xyz_fn;
        if (fn.size() > 4) fn = fn.substr(0, fn.size() - 4) + ".pdb";
        // ipdb handled via a local stream; original used unit 14.
    }
    ixyz.open(chanel_C::xyz_fn);
    if (!debug) chanel_C::iw0 = -1;

    hesinv.assign(molkst_C::nvar * molkst_C::nvar + 1, 0.0);

    size_t ppos = molkst_C::keywrd.find("BIGCYCLES");
    if (ppos != std::string::npos)
        maxcyc = (int)reada(molkst_C::keywrd, (int)ppos);

    bool lef = (molkst_C::keywrd.find("DFP") == std::string::npos) &&
               (molkst_C::nvar > 0);

    alparm.assign(4, std::vector<double>(molkst_C::nvar + 1, 0.0));

    // Restart handling.
    if (lef) {
        if (molkst_C::keywrd.find("RESTAR") != std::string::npos) {
            std::ifstream fs(chanel_C::restart_fn, std::ios::binary);
            if (fs) {
                int i, j;
                fs.read(reinterpret_cast<char*>(&i), sizeof(i));
                fs.read(reinterpret_cast<char*>(&j), sizeof(j));
                if (molkst_C::norbs != j || molkst_C::numat != i) {
                    mopend("Restart file read in does not match current data set");
                    return;
                }
                for (int ii = 1; ii <= molkst_C::nvar; ++ii) {
                    for (int jj = 1; jj <= 3; ++jj) {
                        double v;
                        fs.read(reinterpret_cast<char*>(&v), sizeof(v));
                        alparm[jj][ii] = v;
                    }
                }
                fs.read(reinterpret_cast<char*>(&iloop), sizeof(iloop));
                fs.read(reinterpret_cast<char*>(&x0), sizeof(x0));
                fs.read(reinterpret_cast<char*>(&x1), sizeof(x1));
                fs.read(reinterpret_cast<char*>(&x2), sizeof(x2));
            } else {
                mopend("Restart file is corrupt!");
                return;
            }
        }
    } else {
        if (molkst_C::keywrd.find("RESTAR") != std::string::npos) {
            mdfp[9] = 0;
            std::vector<double> gd(3 * molkst_C::numat + 1, 0.0);
            std::vector<double> xlast(3 * molkst_C::numat + 1, 0.0);
            double totime = 0, funct1 = 0;
            dfpsav(totime, xparam, gd, xlast, funct1, mdfp, xdfp);
        }
    }

    double c1 = (lparam != 1 && na[latom] != 0) ? 57.29577951308232 : 1.0;

    // First reaction step.
    if (iloop <= 1) {
        molkst_C::time0 = second(1);
        if (maxcyc == 0) molkst_C::tleft = -100.0;
        if (lef) ef(xparam, molkst_C::escf);
        else flepo(xparam, molkst_C::nvar, molkst_C::escf);

        size_t ir = molkst_C::keywrd.find("RESTAR");
        if (ir != std::string::npos)
            molkst_C::keywrd.replace(ir, 7, "       ");
        if (iw00 > -1) {
            molkst_C::line = ": " + std::to_string(geo[lparam][latom] * c1) +
                             "  " + std::to_string(molkst_C::escf);
            to_screen(molkst_C::line);
        }
        if (molkst_C::moperr) return;
        if (molkst_C::iflepo == -1) return;
        if (molkst_C::keywrd.find("PDBOUT") != std::string::npos) {
            ++imodel;
            pdbout(14);
        }
        writmo();
        molkst_C::time0 = second(1);
    }

    // Second reaction step.
    if (iloop <= 2) {
        geo[lparam][latom] = react[2];
        if (iloop == 1) {
            x0 = react[1];
            x1 = x0;
            x2 = react[2];
            if (x2 < -100.0) { mopend("Error in PATHS"); return; }
            for (int i = 1; i <= molkst_C::nvar; ++i) alparm[2][i] = xparam[i];
            for (int i = 1; i <= molkst_C::nvar; ++i) alparm[1][i] = xparam[i];
            iloop = 2;
        }
        if (maxcyc == 1) molkst_C::tleft = -100.0;
        if (lef) ef(xparam, molkst_C::escf);
        else flepo(xparam, molkst_C::nvar, molkst_C::escf);
        if (iw00 > -1) {
            molkst_C::line = ": " + std::to_string(geo[lparam][latom] * c1) +
                             "  " + std::to_string(molkst_C::escf);
            to_screen(molkst_C::line);
        }
        if (molkst_C::iflepo == -1) return;
        rxn_coord = react[2];
        if (lparam > 1 && na[latom] > 0) rxn_coord *= 57.29577951308232;

        writmo();
        if (ixyz) {
            ixyz << molkst_C::nl_atoms << " \n";
            ixyz << "PATH \n";
            for (int i = 1; i <= molkst_C::numat; ++i)
                if (l_atom[i])
                    ixyz << elemts_C::elemnt[nat[i]]
                         << "  " << coord[0][i] << "  "
                         << coord[1][i] << "  " << coord[2][i] << "\n";
        }
        if (molkst_C::keywrd.find("PDBOUT") != std::string::npos) {
            ++imodel;
            pdbout(14);
        }
        molkst_C::time0 = second(1);
        for (int i = 1; i <= molkst_C::nvar; ++i) alparm[3][i] = xparam[i];
        if (iloop == 2) iloop = 3;
    }

    // Find number of reaction path points.
    int lpr = iloop;
    int npts = 1;
    while (npts <= 10000 && react[npts] >= -100.0) ++npts;

    for (int ii = lpr; ii <= npts - 1; ++ii) {
        iloop = ii;
        if (iloop - lpr > maxcyc - 3) molkst_C::tleft = -100.0;
        rxn_coord = react[iloop];

        double x3 = react[iloop];
        double c3 = (x0 * x0 - x1 * x1) * (x1 - x2) -
                    (x1 * x1 - x2 * x2) * (x0 - x1);
        double cc1, cc2;
        if (std::abs(c3) < 1e-8) { cc1 = 0; cc2 = 0; }
        else { cc1 = (x1 - x2) / c3; cc2 = (x0 - x1) / c3; }
        double cb1 = 1.0 / (x1 - x2);
        double cb2 = (x1 * x1 - x2 * x2) * cb1;

        for (int i = 1; i <= molkst_C::nvar; ++i) {
            double delf0 = alparm[1][i] - alparm[2][i];
            double delf1 = alparm[2][i] - alparm[3][i];
            double aconst = cc1 * delf0 - cc2 * delf1;
            double bconst = cb1 * delf1 - aconst * cb2;
            double cconst = alparm[3][i] - bconst * x2 - aconst * x2 * x2;
            xparam[i] = cconst + bconst * x3 + aconst * x3 * x3;
            alparm[1][i] = alparm[2][i];
            alparm[2][i] = alparm[3][i];
        }
        for (int i = 1; i <= molkst_C::nvar; ++i) {
            if (std::abs(xparam[i] - alparm[3][i]) > 0.2) {
                for (int k = 1; k <= molkst_C::nvar; ++k) xparam[k] = alparm[3][k];
                break;
            }
        }
        x0 = x1; x1 = x2; x2 = x3;
        geo[lparam][latom] = react[iloop];

        if (lef) ef(xparam, molkst_C::escf);
        else flepo(xparam, molkst_C::nvar, molkst_C::escf);

        if (iw00 > -1) {
            int pc = (int)(100.0 * iloop / npts);
            if (pc != percent) {
                percent = pc;
                molkst_C::line = std::to_string(percent) + "% of Reaction Coordinate done";
                to_screen(molkst_C::line);
            }
            molkst_C::line = ": " + std::to_string(geo[lparam][latom] * c1) +
                             "  " + std::to_string(molkst_C::escf);
            to_screen(molkst_C::line);
        }
        if (molkst_C::iflepo == -1) return;
        if (molkst_C::keywrd.find("PDBOUT") != std::string::npos) {
            ++imodel;
            pdbout(14);
        }
        writmo();
        molkst_C::time0 = second(1);
        for (int i = 1; i <= molkst_C::nvar; ++i) alparm[3][i] = xparam[i];
    }

    chanel_C::iw0 = iw00;
}
