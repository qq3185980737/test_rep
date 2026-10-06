// parsav.cpp — C++ translation of MOPAC 2016 "parsav.F90".
// Saves/restores NLLSQ gradient-minimization restart data (unformatted).
#include "parsav.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace molkst_C {
extern int numat, norbs;
extern std::string keywrd;
}
namespace chanel_C {
extern int iw, ires;
extern std::string restart_fn;
}
namespace common_arrays_C {
extern std::vector<double> aicorr, errfn;
}

extern void mopend(const std::string&);
extern void geout(int);
extern void den_in_out(int);

void parsav(int mode, int& n, int& m, double* q, double* r, double* efslst,
            double* xlast, int* iiium) {
    std::fstream fs;
    {
        std::ios::openmode om = std::ios::binary;
        if (mode == 0) om |= std::ios::in;
        else om |= std::ios::in | std::ios::out | std::ios::trunc;
        fs.open(chanel_C::restart_fn, om);
    }
    if (!fs.is_open()) {
        std::fprintf(stdout, " Restart file either does not exist or is not available for reading\n");
        mopend("Restart file either does not exist or is not available for reading");
        return;
    }
    fs.seekg(0);
    fs.seekp(0);
    int io_stat = 0;
    if (mode == 0) {
        // MODE=0: RETRIEVE DATA FROM DISK.
        int old_numat = 0, old_norbs = 0;
        fs.read((char*)&old_numat, sizeof(int));
        fs.read((char*)&old_norbs, sizeof(int));
        for (int i = 1; i <= n; ++i) fs.read((char*)&xlast[i], sizeof(double));
        fs.read((char*)&m, sizeof(int));
        for (int i = 0; i < 6; ++i) fs.read((char*)&iiium[i], sizeof(int));
        for (int i = 1; i <= n; ++i) fs.read((char*)&efslst[i], sizeof(double));
        fs.read((char*)&n, sizeof(int));
        io_stat = fs.fail() ? 1 : 0;
        if (molkst_C::norbs != old_norbs || molkst_C::numat != old_numat) {
            mopend("Restart file read in does not match current data set");
            return;
        }
        // q(j,i), j=1..m, i=1..m (column-major m x m)
        for (int i = 1; i <= m; ++i)
            for (int j = 1; j <= m; ++j) fs.read((char*)&q[(i - 1) * m + (j - 1)], sizeof(double));
        // r(j,i), j=1..n, i=1..n
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j) fs.read((char*)&r[(i - 1) * n + (j - 1)], sizeof(double));
        if (molkst_C::keywrd.find("AIDER") != std::string::npos) {
            for (int i = 1; i <= n; ++i) fs.read((char*)&common_arrays_C::aicorr[i], sizeof(double));
            for (int i = 1; i <= n; ++i) fs.read((char*)&common_arrays_C::errfn[i], sizeof(double));
        }
        io_stat = fs.fail() ? 1 : 0;
        if (io_stat != 0) {
            std::fprintf(stdout, " Restart file is currupt\n");
            mopend("Restart file is currupt");
        }
        fs.close();
        return;
    }
    if (mode == 1) {
        std::fprintf(stdout, "\n\n          - - - - - - - TIME UP - - - - - - -\n\n");
        std::fprintf(stdout, "\n - THE CALCULATION IS BEING DUMPED TO DISK\n");
        std::fprintf(stdout, "   RESTART IT USING THE KEY-WORD \"RESTART\"\n");
        std::fprintf(stdout, "\n          CURRENT VALUE OF GEOMETRY\n");
        geout(chanel_C::iw);
    }
    fs.write((const char*)&molkst_C::numat, sizeof(int));
    fs.write((const char*)&molkst_C::norbs, sizeof(int));
    for (int i = 1; i <= n; ++i) fs.write((const char*)&xlast[i], sizeof(double));
    fs.write((const char*)&m, sizeof(int));
    for (int i = 0; i < 6; ++i) fs.write((const char*)&iiium[i], sizeof(int));
    for (int i = 1; i <= n; ++i) fs.write((const char*)&efslst[i], sizeof(double));
    fs.write((const char*)&n, sizeof(int));
    for (int i = 1; i <= m; ++i)
        for (int j = 1; j <= m; ++j) fs.write((const char*)&q[(i - 1) * m + (j - 1)], sizeof(double));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j) fs.write((const char*)&r[(i - 1) * n + (j - 1)], sizeof(double));
    if (molkst_C::keywrd.find("AIDER") != std::string::npos) {
        for (int i = 1; i <= n; ++i) fs.write((const char*)&common_arrays_C::aicorr[i], sizeof(double));
        for (int i = 1; i <= n; ++i) fs.write((const char*)&common_arrays_C::errfn[i], sizeof(double));
    }
    den_in_out(1);
    fs.close();
}
