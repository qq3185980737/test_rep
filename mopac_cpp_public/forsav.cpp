// forsav.cpp — C++ translation of MOPAC 2016 "forsav.F90".
// Saves/restores data used in the FORCE calculation (unformatted RESTART file).
#include "forsav.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace molkst_C {
extern int numat, norbs;
}
namespace chanel_C {
extern int iw, ires;
extern std::string restart_fn;
}

extern void mopend(const std::string&);
extern void den_in_out(int);

void forsav(double& time, std::vector<std::vector<double>>& deldip,
            int& ipt, std::vector<double>& fmatrx, std::vector<double>& coord,
            int nvar, double& refh, std::vector<double>& evecs,
            int& jstart, std::vector<double>& fconst) {
    int i99 = 0, j = 0;
retry:
    j = i99;
    // Fortran unit ires may be open; close if so.
    if (ipt == 0) {
        std::ifstream chk(chanel_C::restart_fn, std::ios::binary);
        if (!chk.is_open()) {
            std::fprintf(stdout, " RESTART file does not exist\n");
            mopend("RESTART file does not exist");
            return;
        }
    }
    std::fstream fs;
    {
        std::ios::openmode om = std::ios::binary;
        if (ipt == 0) om |= std::ios::in;
        else om |= std::ios::in | std::ios::out | std::ios::trunc;
        fs.open(chanel_C::restart_fn, om);
    }
    if (!fs.is_open()) {
        if (j != 0) {
            std::fprintf(stdout, " Fatal error in trying to open RESTART file\n");
            mopend("Fatal error in trying to open RESTART file");
            return;
        }
        // A restart file exists but cannot be mounted: delete and retry once.
        std::remove(chanel_C::restart_fn.c_str());
        goto retry;
    }
    fs.seekg(0);
    fs.seekp(0);
    if (ipt == 0) {
        // READ IN FORCE DATA
        int old_numat = 0, old_norbs = 0;
        fs.read((char*)&time, sizeof(double));
        fs.read((char*)&ipt, sizeof(int));
        fs.read((char*)&refh, sizeof(double));
        fs.read((char*)&old_numat, sizeof(int));
        fs.read((char*)&old_norbs, sizeof(int));
        if (molkst_C::norbs != old_norbs || molkst_C::numat != old_numat) {
            mopend("Restart file read in does not match current data set");
            return;
        }
        for (int i = 1; i <= nvar; ++i) fs.read((char*)&coord[i], sizeof(double));
        int linear = (nvar * (nvar + 1)) / 2;
        for (int i = 1; i <= linear; ++i) fs.read((char*)&fmatrx[i], sizeof(double));
        for (int i = 1; i <= ipt; ++i)
            for (int jj = 1; jj <= 3; ++jj) fs.read((char*)&deldip[jj][i], sizeof(double));
        int n33 = nvar * nvar;
        for (int i = 1; i <= n33; ++i) fs.read((char*)&evecs[i], sizeof(double));
        fs.read((char*)&jstart, sizeof(int));
        for (int i = 1; i <= nvar; ++i) fs.read((char*)&fconst[i], sizeof(double));
        if (fs.fail()) {
            mopend("INSUFFICIENT DATA ON DISK FILES FOR A FORCE CALCULATION RESTART.");
            std::fprintf(stdout, "\n          PERHAPS THIS STARTED OFF AS A FORCE CALCULATION\n");
            std::fprintf(stdout, "          BUT THE GEOMETRY HAD TO BE OPTIMIZED FIRST, IN WHICH CASE\n");
            std::fprintf(stdout, "          REMOVE THE KEY-WORD \"FORCE\".\n");
        }
        return;
    }
    // WRITE FORCE DATA
    if (time > 1.e7) time = time - 1.e7;
    fs.write((const char*)&time, sizeof(double));
    fs.write((const char*)&ipt, sizeof(int));
    fs.write((const char*)&refh, sizeof(double));
    fs.write((const char*)&molkst_C::numat, sizeof(int));
    fs.write((const char*)&molkst_C::norbs, sizeof(int));
    int linear = (nvar * (nvar + 1)) / 2;
    for (int i = 1; i <= nvar; ++i) fs.write((const char*)&coord[i], sizeof(double));
    for (int i = 1; i <= linear; ++i) fs.write((const char*)&fmatrx[i], sizeof(double));
    for (int i = 1; i <= ipt; ++i)
        for (int jj = 1; jj <= 3; ++jj) fs.write((const char*)&deldip[jj][i], sizeof(double));
    int n33 = nvar * nvar;
    for (int i = 1; i <= n33; ++i) fs.write((const char*)&evecs[i], sizeof(double));
    fs.write((const char*)&jstart, sizeof(int));
    for (int i = 1; i <= nvar; ++i) fs.write((const char*)&fconst[i], sizeof(double));
    den_in_out(1);
    fs.close();
}
