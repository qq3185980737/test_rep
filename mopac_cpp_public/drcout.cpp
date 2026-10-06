// drcout.cpp — C++ translation of MOPAC 2016 "drcout.F90".
// DRC output: prints geometry/energy line for a DRC at position fract,
// extrapolates xyz/vel/geo quadratically, prints large-geometry block and
// final-geometry block, and writes the trajectory. Fortran 1-based array
// indexing is preserved for all arrays passed in.
#include "drcout.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "maps_C.h"
#include "molkst_C.h"
#include "reada.h"

using namespace common_arrays_C;
using namespace maps_C;
using namespace molkst_C;
using namespace elemts_C;

extern "C" void to_screen_(const char*);
extern void local(double* c, int nocc, double* eig, int iprint,
                  const char* txt);
extern void mullik();
extern void l_control(const std::string& keywrd, int len, int n);
extern void write_trajectory(double* xyz, int mode, double* charge,
                             double escf, double ekin, double time,
                             double xtot);

void drcout(const std::vector<std::vector<double>>& xyz3,
            const std::vector<std::vector<double>>& geo3,
            const std::vector<std::vector<double>>& vel3,
            int nvar, double time,
            const std::vector<double>& escf3, const std::vector<double>& ekin3,
            const std::vector<double>& etot3, const std::vector<double>& xtot3,
            int iloop, const std::vector<double>& charge, double fract,
            const std::string& text1, const std::string& text2,
            int ii, int& jloop) {
    static int icalcn = 0, iprint = 10000;
    static bool drc = false, large = false, graph = false, run_local = false;
    static double last_point = 0.0;
    static double last_rxn_coord = 10.0;
    static int i = 0;

    if (icalcn != numcal) {
        icalcn = numcal;
        last_point = -1e8;
        graph = keywrd.find(" GRAPH") != std::string::npos;
        run_local = keywrd.find(" LOCAL") != std::string::npos;
        if (keywrd.find("RESTART") == std::string::npos ||
            keywrd.find("IRC=") != std::string::npos)
            jloop = 0;
        drc = keywrd.find(" DRC") != std::string::npos;
        i = static_cast<int>(keywrd.find("LARGE"));
        iprint = 10000;
        large = false;
        if (i != -1) {
            iprint = 1;
            i += 5;
            large = (i >= (int)keywrd.size()) || keywrd[i] == ' ' ||
                    (i + 1 < (int)keywrd.size() && keywrd[i + 1] == '-');
            if (i < (int)keywrd.size() && keywrd[i] == '=')
                iprint = static_cast<int>(std::abs(reada(keywrd, i + 1)));
        }
    }
    if (jloop == 0 || (jloop / iprint) * iprint == jloop) {
        if (drc) {
            line = " FEMTOSECONDS  POINT  POTENTIAL + KINETIC  =   TOTAL      ERROR    REF%   MOVEMENT";
        } else {
            line = "     POINT   POTENTIAL  +  ENERGY LOST   =   TOTAL      ERROR    REF%   MOVEMENT";
        }
        std::printf("\n\n%s\n", line.c_str());
    }
    if (drc) {
        if (std::abs(last_point - time) < 5e-4) {
            return;
        }
        last_point = time;
    }
    ++jloop;
    rc_escf = escf3[1] + escf3[2] * fract + escf3[3] * fract * fract;
    ekin = ekin3[1] + ekin3[2] * fract + ekin3[3] * fract * fract;
    double etot = etot3[1] + etot3[2] * fract + etot3[3] * fract * fract;
    rxn_coord = xtot3[1] + xtot3[2] * fract + xtot3[3] * fract * fract;
    if (!drc) {
        if (std::abs(rxn_coord - last_rxn_coord) < 1e-7) {
            --jloop;
            return;
        }
        last_rxn_coord = rxn_coord;
    }
    double errr =
        std::min(9999.99999, std::max(-999.99999, rc_escf + ekin - etot));
    char frmat;
    if (rc_escf > 99999.0 || rc_escf < -9999.0)
        frmat = '3';
    else if (rc_escf > 9999.0 || rc_escf < -999.0)
        frmat = '4';
    else
        frmat = '5';

    char buf[256], fmt[80];
    if (ii != 0) {
        if (drc) {
            std::snprintf(fmt, sizeof(fmt),
                          "%%10.3f%%8d%%12.%cf%%11.5f%%12.%cf%%10.5f %%5d   %%%% %%s%%s%%3d",
                          frmat, frmat);
            std::snprintf(buf, sizeof(buf), fmt, time, iloop - 2, rc_escf,
                          ekin, rc_escf + ekin, errr, jloop, text1.c_str(),
                          text2.c_str(), ii);
        } else {
            std::snprintf(fmt, sizeof(fmt),
                          "%%8d%%14.%cf%%13.5f%%17.5f%%10.5f%%6d   %%%% %%s%%s%%3d",
                          frmat);
            std::snprintf(buf, sizeof(buf), fmt, iloop - 2, rc_escf, ekin,
                          rc_escf + ekin, errr, jloop, text1.c_str(),
                          text2.c_str(), ii);
        }
    } else {
        if (drc) {
            if (text1 == " " && text2 == " ") {
                std::snprintf(fmt, sizeof(fmt),
                              "%%10.3f%%8d%%12.%cf%%11.5f%%12.%cf%%10.5f %%5d   %%%% %%8.4f",
                              frmat, frmat);
                std::snprintf(buf, sizeof(buf), fmt, time, iloop - 2, rc_escf,
                              ekin, rc_escf + ekin, errr, jloop, rxn_coord);
            } else {
                std::snprintf(fmt, sizeof(fmt),
                              "%%10.3f%%8d%%12.%cf%%11.5f%%12.%cf%%10.5f %%5d   %%%% %%s%%s%%3d",
                              frmat, frmat);
                std::snprintf(buf, sizeof(buf), fmt, time, iloop - 2, rc_escf,
                              ekin, rc_escf + ekin, errr, jloop,
                              text1.c_str(), text2.c_str(), ii);
            }
        } else {
            if (text1 == " " && text2 == " ") {
                std::snprintf(fmt, sizeof(fmt),
                              "%%8d%%14.%cf%%13.5f%%17.5f%%10.5f %%5d   %%%% %%8.4f",
                              frmat);
                std::snprintf(buf, sizeof(buf), fmt, iloop - 2, rc_escf, ekin,
                              rc_escf + ekin, errr, jloop, rxn_coord);
            } else {
                std::snprintf(fmt, sizeof(fmt),
                              "%%8d%%14.%cf%%13.5f%%17.5f%%10.5f %%5d   %%%% %%s%%s%%3d",
                              frmat);
                std::snprintf(buf, sizeof(buf), fmt, iloop - 2, rc_escf, ekin,
                              rc_escf + ekin, errr, jloop, text1.c_str(),
                              text2.c_str(), ii);
            }
        }
    }
    line = buf;
    if (keywrd.find(" LDRC_FIRST") != std::string::npos) {
        if (drc)
            line = line.substr(0, 16) + " " + line.substr(17, 36) +
                   "   0.00000     1   %  0.0000";
        else
            line = line.substr(0, 6) + " " + line.substr(7, 45) +
                   "   0.00000     1   %  0.0000";
        jloop = 0;
    }
    std::printf("%s\n", line.c_str());

    natoms = nvar / 3;
    int l = 0;
    std::vector<double> vel(3 * numat + 1, 0.0), xyz(3 * numat + 1, 0.0);
    for (int iat = 1; iat <= natoms; ++iat) {
        for (int j = 1; j <= 3; ++j) {
            vel[(iat - 1) * 3 + j] = vel3[1][l + j] + vel3[2][l + j] * fract +
                                     vel3[3][l + j] * fract * fract;
            xyz[(iat - 1) * 3 + j] = xyz3[1][l + j] + xyz3[2][l + j] * fract +
                                     xyz3[3][l + j] * fract * fract;
        }
        l += 3;
    }
    if (graph) {
        if (run_local) local(&c[1][0], i, &eigs[1], 0, "c ");
        mullik();
    }
    to_screen_("To_file: IRC-DRC");
    if (keywrd.find(" LDRC_FIRST") != std::string::npos) {
        l_control("LDRC_FIRST", 10,
                  -1);
        jloop = 1;
    }
    if (large && (jloop / iprint) * iprint == jloop) {
        std::printf("                CARTESIAN GEOMETRY           VELOCITY (IN CM/SEC)\n");
        std::printf("  ATOM        X          Y          Z                X          Y          Z\n");
        for (int iat = 1; iat <= numat; ++iat) {
            std::printf("%4d   %2s%11.5f%11.5f%11.5f  %11.1f%11.1f%11.1f\n", iat,
                        elemnt[nat[iat]].c_str(), xyz[(iat - 1) * 3 + 1],
                        xyz[(iat - 1) * 3 + 2], xyz[(iat - 1) * 3 + 3],
                        -vel[(iat - 1) * 3 + 1], -vel[(iat - 1) * 3 + 2],
                        -vel[(iat - 1) * 3 + 3]);
        }
    }
    if (drc) {
        write_trajectory(&xyz[0], 1, const_cast<double*>(&charge[0]), rc_escf,
                         ekin, time, rxn_coord);
    } else {
        write_trajectory(&xyz[0], 1, const_cast<double*>(&charge[0]), rc_escf,
                         errr, 0.0, rxn_coord);
    }
    if ((jloop / iprint) * iprint == jloop) {
        int ivar = 1;
        std::printf("\n\n          FINAL GEOMETRY OBTAINED                       CHARGE\n");
        std::printf("%s\n%s\n%s\n", keywrd.c_str(), koment.c_str(),
                    title.c_str());
        l = 0;
        char alpha[3];
        int iel1[4];
        for (int iat = 1; iat <= numat; ++iat) {
            int j = iat / 26;
            alpha[0] = static_cast<char>('A' + j);
            j = iat - j * 26;
            alpha[1] = static_cast<char>('A' + j - 1);
            alpha[2] = '\0';
            iel1[1] = iel1[2] = iel1[3] = 0;
            bool again = true;
            while (again) {
                again = false;
                if (loc[1][ivar] == iat) {
                    iel1[loc[2][ivar]] = 1;
                    ++ivar;
                    again = true;
                }
            }
            if (iat < 4) {
                iel1[3] = 0;
                if (iat < 3) {
                    iel1[2] = 0;
                    if (iat < 2) iel1[1] = 0;
                }
            }
            double gg[4];
            if (labels[iat] < 99 || (labels[iat] > 102 && labels[iat] < 107)) {
                ++l;
                gg[1] = geo3[1][iat * 3 - 2] + geo3[2][iat * 3 - 2] * fract +
                         geo3[3][iat * 3 - 2] * fract * fract;
                gg[2] = geo3[1][iat * 3 - 1] + geo3[2][iat * 3 - 1] * fract +
                         geo3[3][iat * 3 - 1] * fract * fract;
                gg[3] = geo3[1][iat * 3] + geo3[2][iat * 3] * fract +
                         geo3[3][iat * 3] * fract * fract;
                std::printf("  %2s%12.6f%3d%12.6f%3d%12.6f%3d%4d%3d%3d%10.4f%8d%s\n",
                            elemnt[labels[iat]].c_str(), gg[1], iel1[1], gg[2],
                            iel1[2], gg[3], iel1[3], na[iat], nb[iat],
                            nc[iat], charge[l], jloop, alpha);
            } else {
                gg[1] = geo3[1][iat * 3 - 2] + geo3[2][iat * 3 - 2] * fract +
                         geo3[3][iat * 3 - 2] * fract * fract;
                gg[2] = geo3[1][iat * 3 - 1] + geo3[2][iat * 3 - 1] * fract +
                         geo3[3][iat * 3 - 1] * fract * fract;
                gg[3] = geo3[1][iat * 3] + geo3[2][iat * 3] * fract +
                         geo3[3][iat * 3] * fract * fract;
                std::printf("  %2s%12.6f%3d%12.6f%3d%12.6f%3d%4d%3d%3d          %8d%s\n",
                            elemnt[labels[iat]].c_str(), gg[1], iel1[1], gg[2],
                            iel1[2], gg[3], iel1[3], na[iat], nb[iat],
                            nc[iat], jloop, alpha);
            }
        }
    }
}

// reverse_aux(): force.F90 internal — reverse the reaction path already
// written to the trajectory file.  Stub: no-op (keeps force/drc linkable).
void reverse_aux() {
    // TODO: translate force.F90 reverse_aux (trajectory reversal).
}
