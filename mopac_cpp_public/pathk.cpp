// pathk.cpp — C++ translation of "pathk.F90".
// Follow a reaction path with constant step (keywords STEP / POINTS).

#include "pathk.h"

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "maps_C.h"
#include "molkst_C.h"

void ef(std::vector<double>& xparam, double& escf);
void lbfgs(std::vector<double>& xparam, double& escf);
double reada(const std::string& s, int pos);
double second(int);
void geout(int iw);
void pdbout(int unit);
void mopend(const std::string& s);
void to_screen(const std::string& s);
void wrttxt(int unit);
void dfpsav(double& cputot, std::vector<double>& xparam, std::vector<double>& gd,
            std::vector<double>& xlast, double& escf, std::vector<int>& mdfp, std::vector<double>& xdfp);
void write_path_html();
void write_data_to_html(int unit);
void add_path(std::string& line);

void pathk() {
    using namespace common_arrays_C;
    using namespace maps_C;

    std::vector<int> mdfp(21, 0);
    std::vector<double> xdfp(21, 0.0);
    int percent = 0, imodel = 0;

    bool use_lbfgs = (molkst_C::keywrd.find("LBFGS") != std::string::npos) ||
                     (molkst_C::nvar > 2000 &&
                      molkst_C::keywrd.find("EF") == std::string::npos);
    size_t pstep = molkst_C::keywrd.find("STEP");
    double step = (pstep != std::string::npos)
                      ? reada(molkst_C::keywrd, (int)pstep + 5)
                      : 0.0;
    bool debug = molkst_C::keywrd.find("DEBUG") != std::string::npos;
    size_t ppts = molkst_C::keywrd.find("POINT");
    int npts = (ppts != std::string::npos)
                   ? (int)reada(molkst_C::keywrd, (int)ppts + 6)
                   : 0;

    double degree = 57.29577951308232;
    double c1;
    if (lparam != 1 && na[latom] != 0) {
        step /= degree;
        c1 = degree;
    } else {
        c1 = 1.0;
    }

    std::ofstream ixyz(chanel_C::xyz_fn);

    if (molkst_C::keywrd.find("MINI") != std::string::npos) {
        for (auto& row : loc) std::fill(row.begin(), row.end(), 0);
        molkst_C::nvar = 0;
        for (int i = 1; i <= molkst_C::natoms; ++i)
            for (int l = 1; l <= 3; ++l)
                if (lparam != l || latom != i) {
                    ++molkst_C::nvar;
                    loc[1][molkst_C::nvar] = i;
                    loc[2][molkst_C::nvar] = l;
                    xparam[molkst_C::nvar] = geo[l][i];
                }
        grad.assign(molkst_C::nvar + 1, 0.0);
    }

    kloop = 1;
    int maxcyc = 100000;
    size_t pbc = molkst_C::keywrd.find("BIGCYCLES");
    if (pbc != std::string::npos)
        maxcyc = (int)reada(molkst_C::keywrd, (int)pbc);

    double cputot = 0.0;
    rxn_coord = geo[lparam][latom];
    profil.assign(npts + 2, 0.0);
    react.assign(npts + 2, 0.0);
    profil[1] = 0.0;
    react[1] = geo[lparam][latom];

    if (use_lbfgs) {
        if (molkst_C::keywrd.find("RESTAR") != std::string::npos) {
            mdfp[9] = 0;
            std::vector<double> gd(3 * molkst_C::numat + 1, 0.0);
            std::vector<double> xlast(3 * molkst_C::numat + 1, 0.0);
            dfpsav(cputot, xparam, gd, xlast, molkst_C::escf, mdfp, xdfp);
        }
    } else {
        if (molkst_C::keywrd.find("RESTAR") != std::string::npos) {
            std::ifstream fs(chanel_C::restart_fn, std::ios::binary);
            if (fs) {
                int i, l;
                fs.read(reinterpret_cast<char*>(&i), sizeof(i));
                fs.read(reinterpret_cast<char*>(&l), sizeof(l));
                if (molkst_C::norbs != l || molkst_C::numat != i) {
                    mopend("Restart file read in does not match current data set");
                    return;
                }
                fs.read(reinterpret_cast<char*>(&kloop), sizeof(kloop));
                fs.read(reinterpret_cast<char*>(&rxn_coord), sizeof(rxn_coord));
                for (int i = 1; i <= kloop; ++i)
                    fs.read(reinterpret_cast<char*>(&profil[i]), sizeof(double));
            } else {
                mopend("ERROR DETECTED DURING READ IN SUBROUTINE PATHK");
                return;
            }
        }
    }

    geo[lparam][latom] = rxn_coord;
    int lloop = kloop;
    bool scale = false;
    if (molkst_C::id == 1) {
        int i;
        for (i = 1; i <= molkst_C::natoms; ++i)
            if (na[i] != 0) break;
        scale = (i > molkst_C::natoms && latom == molkst_C::natoms);
    }
    int iw00 = chanel_C::iw0;
    if (!debug) chanel_C::iw0 = -1;
    if (molkst_C::keywrd.find("HTML") != std::string::npos) write_path_html();

    for (int iloop = kloop; iloop <= npts; ++iloop) {
        if (iloop - lloop >= maxcyc) molkst_C::tleft = -100.0;
        molkst_C::time0 = second(1);
        double cpu1 = second(2);
        if (iloop > 1 && scale) {
            double factor = geo[lparam][latom] / rxn_coord;
            for (int i = 1; i <= latom - 1; ++i)
                xparam[(i - 1) * 3 + lparam] *= factor;
        }
        rxn_coord = geo[lparam][latom];
        ++molkst_C::numcal;
        if (use_lbfgs) lbfgs(xparam, molkst_C::escf);
        else ef(xparam, molkst_C::escf);

        size_t ir = molkst_C::keywrd.find("RESTAR");
        if (ir != std::string::npos) molkst_C::keywrd.replace(ir, 7, "       ");
        size_t io = molkst_C::keywrd.find("OLDENS");
        if (io != std::string::npos) molkst_C::keywrd.replace(io, 6, "      ");
        if (molkst_C::iflepo == -1 || molkst_C::tleft < 0) return;

        ++kloop;
        double cpu2 = second(2);
        double cpu3 = cpu2 - cpu1;
        cputot += cpu3;
        profil[iloop] = molkst_C::escf;
        if (iw00 > -1) {
            if (!debug) chanel_C::iw0 = 0;
            int pc = (int)(100.0 * iloop / npts);
            if (pc != percent) {
                percent = pc;
                molkst_C::line = ": " + std::to_string(geo[lparam][latom] * c1) +
                                 " " + std::to_string(molkst_C::escf) + " " +
                                 std::to_string(percent) + "% of Reaction Coordinate done";
            } else {
                molkst_C::line = ": " + std::to_string(geo[lparam][latom] * c1) +
                                 " " + std::to_string(molkst_C::escf);
            }
            to_screen(molkst_C::line);
        }
        to_screen("To_file: Reaction path");
        geout(chanel_C::iw);
        if (molkst_C::keywrd.find("PDBOUT") != std::string::npos) {
            ++imodel;
            if (molkst_C::ncomments == 0) molkst_C::ncomments = 1;
            all_comments[1] = " HEADER Heat of Formation = " +
                              std::to_string(molkst_C::escf) + " Kcal/mol";
            pdbout(14);
        }
        geo[lparam][latom] += step;

        if (ixyz) {
            ixyz << molkst_C::nl_atoms << " \n";
            ixyz << "Profile. " << iloop
                 << " HEAT OF FORMATION =" << molkst_C::escf
                 << " KCAL =" << molkst_C::escf * 4.184
                 << " KJ\n";
            for (int i = 1; i <= molkst_C::numat; ++i)
                if (l_atom[i])
                    ixyz << elemts_C::elemnt[nat[i]]
                         << "  " << coord[0][i] << "  "
                         << coord[1][i] << "  " << coord[2][i] << "\n";
        }
        if (!debug) chanel_C::iw0 = -1;
    }

    if (cputot > 1e7) cputot -= 1e7;
    react[1] *= c1;
    double stepc1 = step * c1;
    for (int i = 2; i <= npts; ++i) react[i] = react[i - 1] + stepc1;

    // Archive file.
    std::ofstream iarc(chanel_C::archive_fn, std::ios::app);
    if (iarc) {
        iarc << " ARCHIVE FILE FOR PATH CALCULATION\n";
        iarc << "A PROFILE OF COORDINATES - HEATS\n\n";
        wrttxt(chanel_C::iw);  // text block to main iw; archive handled in Fortran
        iarc << "\n TOTAL JOB TIME : " << cputot << "\n";
        for (int i = 1; i <= npts; ++i)
            iarc << react[i] << " " << profil[i] << "\n";
    }
    if (!debug) chanel_C::iw0 = iw00;
}

// write_path_html: pathk.F90 internal subroutine; M07 batch will provide full translation.
void write_path_html() {}
