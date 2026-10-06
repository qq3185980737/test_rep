// dfpsav.cpp — C++ translation of MOPAC 2016 "dfpsav.F90".
// Stores/restores DFP geometry-optimisation state to the binary restart
// file: norbs/numat, xparam, gd, mdfp/xdfp/totime/funct1, xlast, grad,
// hesinv (linear), and the ef_C/maps extension blocks when latom != 0.
// Fortran 1-based indexing preserved. Dump via ofstream, restore via
// ifstream (Fortran open(UNFORMATTED)+rewind semantics).
#include "dfpsav.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "ef_C.h"
#include "maps_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace maps_C;

extern void geout(int);
extern void mopend(const char* message);
extern void prtgra();
extern void den_in_out(int mode);

void dfpsav(double& totime, std::vector<double>& xparam,
            std::vector<double>& gd, std::vector<double>& xlast,
            double& funct1, std::vector<int>& mdfp,
            std::vector<double>& xdfp) {
    static int icalcn = 0;
    const bool first = icalcn != numcal;
    if (first) icalcn = numcal;

    if (mdfp[9] != 0) {
        if (mdfp[9] == 1) {
            std::printf("\n\n- - - - - - - - TIME UP - - - - - - -\n\n");
            if (keywrd.find("SADDLE") != std::string::npos) {
                std::printf("\n\n          NO RESTART EXISTS FOR SADDLE\n\n"
                            "          HERE IS A DATA-FILE FILES THAT MIGHT BE SUITABLE\n\n"
                            "          FOR RESTARTING THE CALCULATION\n\n");
                std::printf("%s\n%s\n%s\n", keywrd.c_str(), koment.c_str(),
                            title.c_str());
                for (int loop = 1; loop <= 2; ++loop) {
                    geout(-chanel_C::iw);
                    for (int ia = 1; ia <= natoms; ++ia)
                        for (int j = 1; j <= 3; ++j)
                            geo[j][ia] = geoa[j][ia];
                    na.assign(natoms + 1, 0);
                }
                std::printf("\n\n\n          CALCULATION TERMINATED HERE\n");
                mopend("NO RESTART EXISTS FOR SADDLE. ");
                return;
            }
            std::printf("\n\n          - THE CALCULATION IS BEING DUMPED TO DISK\n"
                        "            RESTART IT USING THE KEYWORD \"RESTART\"\n");
            if (keywrd.find("STEP1") == std::string::npos) {
                std::printf("\n\n          CURRENT VALUE OF HEAT OF FORMATION =%12.6f\n",
                            funct1);
                if (prt_gradients && keywrd.find(" GRADI") != std::string::npos &&
                    mozyme) {
                    std::printf("\n\n\n          CURRENT  POINT  AND  DERIVATIVES\n\n");
                    prtgra();
                }
                if (mdfp[9] == 1) geout(chanel_C::iw);
            }
        }
        {
            std::ofstream of(chanel_C::restart_fn,
                             std::ios::binary | std::ios::trunc);
            if (of) {
                of.write(reinterpret_cast<const char*>(&norbs), sizeof(int));
                of.write(reinterpret_cast<const char*>(&numat), sizeof(int));
                for (int i = 1; i <= nvar; ++i)
                    of.write(reinterpret_cast<const char*>(&xparam[i]),
                             sizeof(double));
                for (int i = 1; i <= nvar; ++i)
                    of.write(reinterpret_cast<const char*>(&gd[i]),
                             sizeof(double));
                for (int i = 1; i <= 9; ++i)
                    of.write(reinterpret_cast<const char*>(&mdfp[i]),
                             sizeof(int));
                for (int i = 1; i <= 9; ++i)
                    of.write(reinterpret_cast<const char*>(&xdfp[i]),
                             sizeof(double));
                of.write(reinterpret_cast<const char*>(&totime),
                         sizeof(double));
                of.write(reinterpret_cast<const char*>(&funct1),
                         sizeof(double));
                for (int i = 1; i <= nvar; ++i)
                    of.write(reinterpret_cast<const char*>(&xlast[i]),
                             sizeof(double));
                for (int i = 1; i <= nvar; ++i)
                    of.write(reinterpret_cast<const char*>(&grad[i]),
                             sizeof(double));
                const int linear = nvar * (nvar + 1) / 2;
                if ((int)hesinv.size() < linear + 1)
                    hesinv.resize(linear + 1, 0.0);
                for (int i = 1; i <= linear; ++i)
                    of.write(reinterpret_cast<const char*>(&hesinv[i]),
                             sizeof(double));
            }
            den_in_out(1);
            if (latom != 0) {
                if (keywrd.find(" STEP=") != std::string::npos) {
                    if (of) {
                        of.write(reinterpret_cast<const char*>(&kloop),
                                 sizeof(int));
                        of.write(reinterpret_cast<const char*>(&rxn_coord),
                                 sizeof(double));
                        for (int i = 1; i <= kloop; ++i)
                            of.write(reinterpret_cast<const char*>(&profil[i]),
                                     sizeof(double));
                    }
                } else {
                    if (of) {
                        for (int i = 1; i <= nvar; ++i)
                            for (int j = 1; j <= 3; ++j)
                                of.write(reinterpret_cast<const char*>(
                                             &ef_C::alparm[j][i]),
                                         sizeof(double));
                        of.write(reinterpret_cast<const char*>(&ef_C::iloop),
                                 sizeof(int));
                        of.write(reinterpret_cast<const char*>(&ef_C::x0),
                                 sizeof(double));
                        of.write(reinterpret_cast<const char*>(&ef_C::x1),
                                 sizeof(double));
                        of.write(reinterpret_cast<const char*>(&ef_C::x2),
                                 sizeof(double));
                    }
                }
            }
            if (keywrd.find("STEP1") != std::string::npos) return;
        }
    } else {
        if (first) std::printf("\n\n          RESTORING DATA FROM DISK\n\n");
        int old_norbs = 0, old_numat = 0;
        std::ifstream inf(chanel_C::restart_fn, std::ios::binary);
        bool io_ok = inf.good();
        if (!io_ok) {
            mopend("NO RESTART FILE EXISTS!");
            return;
        }
        if (io_ok) {
            io_ok = inf.read(reinterpret_cast<char*>(&old_norbs), sizeof(int))
                        .good() &&
                    inf.read(reinterpret_cast<char*>(&old_numat), sizeof(int))
                        .good();
            for (int i = 1; i <= nvar && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&xparam[i]),
                                 sizeof(double)).good();
            for (int i = 1; i <= nvar && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&gd[i]),
                                 sizeof(double)).good();
            if (!io_ok) {
                mopend("RESTART FILE EXISTS, BUT IS CORRUPT");
                return;
            }
        }
        if (norbs != old_norbs || numat != old_numat) {
            mopend("Restart file read in does not match current data set");
            return;
        }
        if (io_ok) {
            io_ok = true;
            for (int i = 1; i <= 9 && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&mdfp[i]),
                                 sizeof(int)).good();
            for (int i = 1; i <= 9 && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&xdfp[i]),
                                 sizeof(double)).good();
            io_ok = inf.read(reinterpret_cast<char*>(&totime), sizeof(double))
                        .good() &&
                    inf.read(reinterpret_cast<char*>(&funct1), sizeof(double))
                        .good();
            if (!io_ok) {
                mopend("NO RESTART FILE EXISTS!");
                return;
            }
        }
        if (first) std::printf("          FUNCTION =%13.6f\n\n", funct1);
        if (io_ok) {
            io_ok = true;
            for (int i = 1; i <= nvar && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&xlast[i]),
                                 sizeof(double)).good();
            for (int i = 1; i <= nvar && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&grad[i]),
                                 sizeof(double)).good();
            const int linear = nvar * (nvar + 1) / 2;
            if ((int)hesinv.size() < linear + 1)
                hesinv.resize(linear + 1, 0.0);
            for (int i = 1; i <= linear && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&hesinv[i]),
                                 sizeof(double)).good();
            if (!io_ok) {
                mopend("RESTART FILE EXISTS, BUT IS CORRUPT");
                return;
            }
        }
        den_in_out(0);
        if (latom != 0) {
            if (keywrd.find(" STEP=") != std::string::npos) {
                if (io_ok) {
                    io_ok = inf.read(reinterpret_cast<char*>(&kloop),
                                      sizeof(int)).good() &&
                            inf.read(reinterpret_cast<char*>(&rxn_coord),
                                     sizeof(double)).good();
                    for (int i = 1; i <= kloop && io_ok; ++i)
                        io_ok = inf.read(reinterpret_cast<char*>(&profil[i]),
                                         sizeof(double)).good();
                }
            } else {
                if (io_ok) {
                    io_ok = true;
                    for (int i = 1; i <= nvar; ++i)
                        for (int j = 1; j <= 3; ++j)
                            io_ok = inf.read(reinterpret_cast<char*>(
                                                 &ef_C::alparm[j][i]),
                                             sizeof(double)).good();
                    io_ok = inf.read(reinterpret_cast<char*>(&ef_C::iloop),
                                     sizeof(int)).good() &&
                            inf.read(reinterpret_cast<char*>(&ef_C::x0),
                                     sizeof(double)).good() &&
                            inf.read(reinterpret_cast<char*>(&ef_C::x1),
                                     sizeof(double)).good() &&
                            inf.read(reinterpret_cast<char*>(&ef_C::x2),
                                     sizeof(double)).good();
                }
            }
            if (!io_ok) {
                mopend("RESTART FILE EXISTS, BUT IS CORRUPT");
                return;
            }
        }
        return;
    }
}
