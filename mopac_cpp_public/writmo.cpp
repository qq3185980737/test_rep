// writmo.cpp — C++ translation of "writmo.F90" (partial).
// WRITMO prints out most of the job results. It must not alter parameters.
//
// Implemented: results-summary header, optimization/SCF status, pressure
// correction, reference-geometry distortion block, FINAL HEAT OF FORMATION,
// and DISP total-energy block. Remaining ~2500 lines (geometry/gradient/
// eigenvalue/charge/COSMO/MOZYME orbital tables) are documented stubs pending
// the COSMO and MOZYME output globals.

#include "writmo.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "big_swap.h"          // l_control
#include "bonds.h"
#include "denrot.h"
#include "dimens.h"
#include "dipole.h"
#include "dot.h"
#include "elemts_C.h"
#include "chrge.h"
#include "enpart.h"
#include "meci_C.h"
#include "mecip.h"
#include "mullik.h"
#include "parameters_C.h"
#include "post_scf_corrections.h"
#include "funcon_C.h"
#include "maps_C.h"
#include "molkst_C.h"
#include "timout.h"
#include "to_screen.h"
#include "vecprt.h"
#include "volume.h"
#include "wrttxt.h"
#include "symmetry_C.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace cosmo_C;
using namespace maps_C;
using namespace molkst_C;

// Not-yet-ported helpers called by writmo.
void geout(int mode);
void PM7_TS();
void fdate(std::string& s);

// flepo(1..17): geometry-optimization outcome messages.
static const char* flepo[18] = {
    "",
    " 1SCF WAS SPECIFIED, SO BFGS WAS NOT USED",
    " GRADIENTS WERE INITIALLY ACCEPTABLY SMALL",
    " HERBERTS TEST WAS SATISFIED IN BFGS",
    " THE LINE MINIMIZATION FAILED TWICE IN A ROW.   TAKE CARE!",
    " BFGS FAILED DUE TO COUNTS EXCEEDED. TAKE CARE!",
    " PETERS TEST WAS SATISFIED IN BFGS OPTIMIZATION",
    " THIS MESSAGE SHOULD NEVER APPEAR, THERE IS A BUG IN MOPAC",
    " GRADIENT TEST NOT PASSED, BUT FURTHER WORK NOT JUSTIFIED",
    " A FAILURE HAS OCCURRED, TREAT RESULTS WITH CAUTION!!",
    " GEOMETRY OPTIMIZED USING NLLSQ. GRADIENT NORM MINIMIZED",
    " GEOMETRY OPTIMIZED USING POWSQ. GRADIENT NORM MINIMIZED",
    " CYCLES EXCEEDED, GRADIENT NOT FULLY MINIMIZED IN NLLSQ",
    " 1SCF RUN AFTER RESTART.  GEOMETRY MIGHT NOT BE OPTIMIZED",
    " HEAT OF FORMATION MINIMIZED IN ONE LINE SEARCH",
    " GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (EF).",
    " NO PARAMETERS MARKED FOR OPTIMIZATION, SO 1SCF WAS USED",
    " GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (TS)."};
static const char* iter_msg[3] = {"", " SCF FIELD WAS ACHIEVED",
                                  " ++++----**** FAILED TO ACHIEVE SCF. ****----++++"};

void writmo() {
    bool lprtgra = keywrd.find(" GRAD") != std::string::npos && nvar > 0;
    (void)lprtgra;

    std::string caltyp = "   ";
    if (method_am1) caltyp = "       AM1";
    else if (method_pm3) caltyp = "       PM3";
    else if (method_mndo) caltyp = "      MNDO";
    else if (method_mndod) caltyp = "     MNDO-D";
    else if (method_rm1) caltyp = "       RM1";
    else if (method_pm6) caltyp = "       PM6";
    else if (method_pm7) caltyp = "       PM7";
    // left-trim caltyp
    size_t sp = caltyp.find_first_not_of(' ');
    if (sp != std::string::npos) caltyp = caltyp.substr(sp);

    std::printf("\n ----");
    for (int i = 0; i < 15; ++i) std::printf("-----");
    std::printf("\n");
    wrttxt(iw);

    double gnorm = 0.0;
    if (nvar != 0) gnorm = std::sqrt(dot(grad, grad, nvar));
    (void)gnorm;

    if (iflepo == 15 && keywrd.find(" TS") != std::string::npos) iflepo = 17;
    if (iflepo == 0) iflepo = 7;
    if (iflepo == 1 && nvar == 0 && keywrd.find("1SCF") == std::string::npos)
        iflepo = 16;
    std::printf("\n\n%58s\n", flepo[iflepo]);
    if (iscf < 1) iscf = 1;
    std::printf("%58s\n", iter_msg[iscf]);
    std::printf("\n\n%30s %s CALCULATION\n", "", caltyp.c_str());
    std::printf("%55sMOPAC2016 (Version: %s%s)\n", "", verson.c_str(), ")");

    std::string idate;
    fdate(idate);
    std::printf("%-55s\n", idate.c_str());
    if (ijulian < 500) std::printf("%55s No. of days remaining = %d\n", "", ijulian);

    if (iscf == 2) {
        std::printf(" \n \n FOR SOME REASON THE SCF CALCULATION FAILED.\n \n");
        std::printf(" THE RESULTS WOULD BE MEANINGLESS, SO WILL NOT BE PRINTED.\n");
        std::printf(" TRY TO FIND THE REASON FOR THE FAILURE BY USING 'PL'.\n \n");
        std::printf(" CHECK YOUR GEOMETRY AND ALSO TRY USING SHIFT OR PULAY. \n");
        geout(1);
        return;
    }

    // External pressure contribution.
    if (std::abs(pressure) > 1e-4) {
        if (id == 1) {
            double ln = 0.0;
            for (int d = 0; d < 3; ++d) ln += tvec[d][0] * tvec[d][0];
            escf = escf + pressure * std::sqrt(ln);
        } else if (id == 3) {
            double tv[10];
            for (int d = 0; d < 3; ++d)
                for (int j = 0; j < 3; ++j) tv[d * 3 + j] = tvec[d][j];
            escf = escf + pressure * volume(tv, 3);
        }
    }

    if (use_ref_geo) {
        std::printf("\n\n%10sFINAL H.O.F PLUS STRESS =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", escf, escf * 4.184);
        stress = stress * density;
        std::printf("%10sFINAL STRESS            =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", stress, stress * 4.184);
        std::printf("%10sFINAL HEAT OF FORMATION =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", escf - stress, (escf - stress) * 4.184);
        double sum = 0.0, rms = 0.0;
        for (int i = 1; i <= numat; ++i) {
            double dx = geo[0][i] - geoa[0][i];
            double dy = geo[1][i] - geoa[1][i];
            double dz = geo[2][i] - geoa[2][i];
            sum += std::sqrt(dx * dx + dy * dy + dz * dz);
            rms += dx * dx + dy * dy + dz * dz;
        }
        double distortion = sum / numat;
        std::printf("%10sTOTAL DISTORTION        =%17.5f Angstroms\n", "", sum);
        std::printf("%10sAVERAGE DISTORTION      =%17.5f Angstroms per atom (all atoms)\n",
                    "", distortion);
        std::printf("%10sRMS DISTORTION          =%17.5f Angstroms per atom (all atoms)\n",
                    "", std::sqrt(rms / numat));
    } else {
        if (keywrd.find(" PM7-TS") != std::string::npos) {
            PM7_TS();
            return;
        }
        stress = -1.0;
        std::printf("\n\n%10sFINAL HEAT OF FORMATION =%17.5f KCAL/MOL  =%14.5f KJ/MOL\n",
                    "", escf, escf * 4.184);
    }

    if (keywrd.find(" DISP") != std::string::npos) {
        // Full DISP energy decomposition (F90 writmo L205-230).
        std::printf("\n%10sTOTAL ENERGY            =%17.5f KCAL/MOL\n", "",
                    (elect + enuclr + solv_energy) * funcon_C::fpc_9);
        std::printf("%10sENERGY OF ATOMS         =%17.5f KCAL/MOL\n", "", atheat);
        std::printf("%10s                    SUM =%17.5f KCAL/MOL\n", "",
                    (elect + enuclr) * funcon_C::fpc_9 + atheat +
                        solv_energy * funcon_C::fpc_9);
        if (std::abs(hpress) > 1.0e-5)
            std::printf("%10sENERGY DUE TO PRESSURE  =%17.5f KCAL/MOL\n", "", hpress);
        std::printf("%10sDISPERSION ENERGY       =%17.5f KCAL/MOL\n", "", E_disp);
        if (E_hb < -1.0e-5)
            std::printf("%10sH-BOND ENERGY           =%17.5f KCAL/MOL\n", "", E_hb);
        if (std::abs(nsp2_corr) > 1.0e-5)
            std::printf("%10sMM CORRECTION FOR >N-   =%17.5f KCAL/MOL\n", "", nsp2_corr);
        if (std::abs(Si_O_H_corr) > 1.0e-5)
            std::printf("%10sMM CORR. FOR Si-O-H     =%17.5f KCAL/MOL\n", "", Si_O_H_corr);
        if (std::abs(sum_dihed) > 1.0e-5)
            std::printf("%10sMM CORR. FOR -CO-NH-    =%17.5f KCAL/MOL\n", "", sum_dihed);
        double dsum = (elect + enuclr) * funcon_C::fpc_9 + atheat + hpress +
                      solv_energy * funcon_C::fpc_9 + nsp2_corr +
                      Si_O_H_corr + sum_dihed + E_disp + E_hb;
        std::printf("%30sSUM =%17.5f KCAL/MOL\n\n", "", dsum);
        if (std::abs(dsum - escf) > 1.0e-2)
            std::printf(" WARNING - An energy term is missing!\n");
        if (N_Hbonds > 0)
            std::printf("%10sNo. OF HYDROGEN BONDS   =%11d%7s\n", "", N_Hbonds,
                        "(H-bond Energy < -1.0 Kcal/mol)");
        if (E_hh > 1.0e-5)
            std::printf("%10sH - H CORRECTION ENERGY =%17.5f KCAL/MOL\n", "", E_hh);
        if (keywrd.find(" DISP(") != std::string::npos) {
            l_control("PRT", 3, 1);
            double corr = 0.0;
            post_scf_corrections(corr, false);
        }
    }
    // F90 L231-239: large-HOF warning and screen messages.
    if (numat > 1 && iscf == 1 && escf > 1.e4 &&
        keywrd.find(" CHECK") == std::string::npos)
        std::printf("\n\n%10sCalculated Heat of Formation is very large, re-run using keyword 'CHECK'\n\n", "");
    to_screen(" Job: " + jobnam);
    {
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "Final heat of formation = %16.5f kcal/mol", escf);
        to_screen(buf);
        std::snprintf(buf, sizeof(buf),
                      "GRADIENT NORM           =%17.5f", gnorm);
        to_screen(buf);
    }

    if (keywrd.find(" EPS") != std::string::npos)
        std::printf("%10sVAN DER WAALS AREA      =%14.2f SQUARE ANGSTROMS\n", "", area);

    if (symmetry_C::state_Irred_Rep != "    ")
        std::printf("%10sTOTAL ENERGY            =%17.5f EV   STATE:  %2d %s %s\n",
                    "", elect + enuclr + solv_energy, symmetry_C::state_QN,
                    symmetry_C::state_spin.c_str(), symmetry_C::state_Irred_Rep.c_str());
    else
        std::printf("%10sTOTAL ENERGY            =%17.5f EV\n", "",
                    elect + enuclr + solv_energy);
    if (id == 0 && !mozyme)
        std::printf("%10sELECTRONIC ENERGY       =%17.5f EV  POINT GROUP:  %-4s\n",
                    "", elect, symmetry_C::name.c_str());
    else
        std::printf("%10sELECTRONIC ENERGY       =%17.5f EV\n", "", elect);
    std::printf("%10sCORE-CORE REPULSION     =%17.5f EV\n", "", enuclr);
    if (iseps)
        std::printf("%10sDIELECTRIC ENERGY       =%17.5f EV\n", "", ediel);
    if (fepsi > 1e-10 && area > 1e-3) {
        std::printf("%10sCOSMO AREA              =%14.2f SQUARE ANGSTROMS\n", "", area);
        std::printf("%10sCOSMO VOLUME            =%14.2f CUBIC ANGSTROMS\n", "", cosvol);
    }
    std::printf("%10sGRADIENT NORM           =%17.5f\n", "", gnorm);
    // F90 L345-369: stationary-point check + NOANCI hint.
    bool still = true;
    if (latom == 0) {
        if (keywrd.find(" AIDER") == std::string::npos) {
            if (keywrd.find(" 1SCF") == std::string::npos ||
                keywrd.find(" GRAD") != std::string::npos) {
                double gsum = std::sqrt(dot(dxyz, dxyz, 3 * numat));
                if (gsum > std::max(5.0, 5 * gnorm * gnorm) &&
                    gnorm < 2.0 && nclose == nopen && id == 0 &&
                    keywrd.find("NOANCI") == std::string::npos) {
                    if (nvar != 1 ||
                        keywrd.find(" GRAD") != std::string::npos ||
                        keywrd.find("DERIV") != std::string::npos) {
                        std::printf("%9s WARNING -- GEOMETRY IS NOT AT A STATIONARY POINT\n", "");
                        still = false;
                    }
                }
            }
        }
    }
    // Ionization potential and HOMO/LUMO energies.
    double eionis = 0.0;
    if (nalpha > 0 && nbeta > 0)
        eionis = -std::max(eigs[nalpha], eigb[nbeta]);
    else if (nelecs == 1)
        eionis = -eigs[1];
    else if (nelecs > 1) {
        if (nclose > 0) eionis = -eigs[nclose];
        if (nopen > 0) eionis = std::min(eionis, -eigs[nopen]);
    }
    {
        int iop = nclose;
        if (fract > 1.99) iop = nopen;
        int nopn = nopen - iop;
        if (!mozyme) {
            if (nopn == 1 && meci_C::rjkab.size() > 1)
                eionis += 0.5 * meci_C::rjkab[1][1];  // doublet IP correction
            if (std::abs(eionis) > 1e-5 && nopn < 2)
                std::printf("%10sIONIZATION POTENTIAL    =  %16.6f EV\n", "", eionis);
            if (uhf) {
                if (nalpha >= 1)
                    std::printf("%10sALPHA SOMO LUMO (EV)    =  %13.3f %7.3f\n", "",
                                eigs[nalpha],
                                (nalpha + 1 <= norbs) ? eigs[nalpha + 1] : 0.0);
                if (nbeta >= 1)
                    std::printf("%10sBETA  SOMO LUMO (EV)    =  %13.3f %7.3f\n", "",
                                eigb[nbeta],
                                (nbeta + 1 <= norbs) ? eigb[nbeta + 1] : 0.0);
            } else if (nopen == nclose) {
                if (nopen >= 1)
                    std::printf("%10sHOMO LUMO ENERGIES (EV) =  %13.3f %7.3f\n", "",
                                eigs[nopen],
                                (nopen + 1 <= norbs) ? eigs[nopen + 1] : 0.0);
            } else if (nopn == 1) {
                if (nclose >= 1)
                    std::printf("%10sHOMO (SOMO) LUMO (EV)   =%13.3f (%7.3f) %7.3f\n",
                                "", eigs[nclose],
                                (nclose + 1 <= norbs) ? eigs[nclose + 1] : 0.0,
                                (nclose + 2 <= norbs) ? eigs[nclose + 2] : 0.0);
                else
                    std::printf("%10s     (SOMO) LUMO (EV)   =      (%6.3f)%7.3f\n",
                                "",
                                (nclose + 1 <= norbs) ? eigs[nclose + 1] : 0.0,
                                (nclose + 2 <= norbs) ? eigs[nclose + 2] : 0.0);
            }
        }
        if (uhf) {
            std::printf("%10sNO. OF ALPHA ELECTRONS  =%11d\n", "",
                        nalpha + (int)((nalpha_open - nalpha) * fract + 0.5));
            std::printf("%10sNO. OF BETA  ELECTRONS  =%11d\n", "",
                        nbeta + (int)((nbeta_open - nbeta) * fract + 0.5));
        } else {
            std::printf("%10sNO. OF FILLED LEVELS    =%11d\n", "", nopen - nopn);
            if (nopn != 0)
                std::printf("%10sAND NO. OF OPEN LEVELS  =%11d\n", "", nopn);
        }
    }
    if (mol_weight > 0.1)
        std::printf("%10sMOLECULAR WEIGHT        =%16.4f\n", "", mol_weight);
    dimens(coord, iw);   // F90: call dimens(coord, iw) after MOLECULAR WEIGHT
    std::printf("%10sSCF CALCULATIONS        =   %8d\n", "", nscf);
    timout(iw);   // F90: call timout(iw) after SCF CALCULATIONS
    // Sync optimized geometry back into xparam.
    for (int i = 1; i <= nvar; ++i)
        xparam[i] = geo[loc[2][i]][loc[1][i]];

    // FINAL POINT AND DERIVATIVES table.
    if (prt_gradients && (lprtgra || gnorm > 2.0)) {
        std::printf("\n\n\n%7sFINAL  POINT  AND  DERIVATIVES\n\n", "");
        if (mozyme) {
            // prtgra() not ported.
        } else {
            std::printf("   PARAMETER     ATOM    TYPE            VALUE       GRADIENT\n");
            const double degree = 57.29577951308232;
            for (int i = 1; i <= nvar; ++i) {
                int j = loc[2][i], k = loc[1][i];
                int l = labels[k];
                double xi = xparam[i];
                if (j != 1 && na[k] > 0) xi *= degree;
                const char* gtype = (j == 1 || na[k] == 0) ? "KCAL/ANGSTROM"
                                                          : "KCAL/RADIAN  ";
                const char* type;
                if (na[k] == 0)
                    type = (j == 1) ? "CARTESIAN X" : (j == 2) ? "CARTESIAN Y"
                                                                  : "CARTESIAN Z";
                else
                    type = (j == 1) ? "BOND       " : (j == 2) ? "ANGLE      "
                                                                : "DIHEDRAL   ";
                std::printf("%7d%11d  %s   %-11s%13.6f%13.6f  %s\n",
                            i, k, elemts_C::elemnt[l].c_str(), type, xi, grad[i], gtype);
            }
        }
    }

    // CARTESIAN COORDINATES.
    if (prt_cart) {
        std::printf("\n\n\n%28s CARTESIAN COORDINATES\n", "");
        for (int i = 1; i <= numat; ++i)
            std::printf("%4d   %-2s   %16.9f%16.9f%16.9f\n", i,
                        elemts_C::elemnt[nat[i]].c_str(),
                        coord[0][i], coord[1][i], coord[2][i]);
    }
    // PRTINT: interatomic distances (F90 L527-541).
    if (!mozyme && keywrd.find("PRTINT") != std::string::npos) {
        int npr = numat * (numat + 1) / 2;
        std::vector<double> rxyz(npr, 0.0);
        int l = -1;
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= i; ++j) {
                ++l;
                double dx = coord[0][i] - coord[0][j];
                double dy = coord[1][i] - coord[1][j];
                double dz = coord[2][i] - coord[2][j];
                rxyz[l] = std::sqrt(dx * dx + dy * dy + dz * dz);
            }
        std::printf("\n\n%10s  INTERATOMIC DISTANCES\n", "");
        vecprt(&rxyz[0], numat);
    }

    // Clip wild MO eigenvalues.
    if (!mozyme) {
        for (int i = 1; i <= norbs; ++i) {
            if (!(eigs[i] >= -999.0 && eigs[i] <= 1000.0)) eigs[i] = 1e-12;
            if (eigb[i] < -999.0 || eigb[i] > 1000.0) eigb[i] = 1e-12;
        }
    }
    if (!formula.empty()) std::printf("\n%s\n\n", formula.c_str());

    // Hessian matrix (upper triangle).
    if (keywrd.find(" HESSIAN") != std::string::npos) {
        std::printf("    Hessian matrix from geometry optimization\n");
        for (int i = 1; i <= nvar; ++i) {
            for (int j = 1; j <= i; ++j)
                std::printf("%10.4f", hesinv[(i - 1) * nvar + j - 1]);
            std::printf("\n");
        }
    }
    // MOLECULAR POINT GROUP.
    if (id == 0 && !mozyme)
        std::printf("\n\n      MOLECULAR POINT GROUP   :   %-4s\n",
                    symmetry_C::name.c_str());

    // EIGENVALUES list (unless eigenvectors requested via VEC/ALLVEC).
    if (norbs > 0) {
        if (!(keywrd.find(" VEC") != std::string::npos ||
              keywrd.find(" ALLVEC") != std::string::npos)) {
            if (uhf) std::printf("\n\n%10sALPHA EIGENVALUES\n", "");
            else     std::printf("\n\n%18sEIGENVALUES\n", "");
            for (int i = 1; i <= norbs; ++i) {
                std::printf("%10.5f", eigs[i]);
                if (i % 8 == 0) std::printf("\n");
            }
            if (norbs % 8 != 0) std::printf("\n");
            if (uhf) {
                std::printf("\n\n%10s BETA EIGENVALUES\n", "");
                for (int i = 1; i <= norbs; ++i) {
                    std::printf("%10.5f", eigb[i]);
                    if (i % 8 == 0) std::printf("\n");
                }
                if (norbs % 8 != 0) std::printf("\n");
            }
        }
    }

    // NET ATOMIC CHARGES.
    if (nelecs != 0) {
        // Correct density matrix, if necessary (F90 writmo.F90 L605-610).
        if ((nclose != nopen && std::fabs(fract - 2.0) > 1e-20 && fract > 1e-20) ||
            keywrd.find(" C.I.") != std::string::npos)
            mecip();
        if (prt_charges) {
            std::printf("\n\n%13sNET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS\n\n", "");
            // F90 L614-625: table header with s-/p-/d-Pop labels.
            int nb = 0;
            for (int j = 1; j <= numat; ++j)
                nb = std::max(nb, parameters_C::natorb[nat[j]]);
            nb = (int)std::lround(std::sqrt((double)nb)) * 12;
            const char* pops = "s-Pop       p-Pop       d-Pop";
            if (maxtxt == 26)
                std::printf(" ATOM NO.                 TYPE                        CHARGE      No. of ELECS.   %s\n",
                            std::string(pops).substr(0, std::min(nb, 36)).c_str());
            else
                std::printf(" ATOM NO.   TYPE          CHARGE      No. of ELECS.   %s\n",
                            std::string(pops).substr(0, std::min(nb, 36)).c_str());
        }
        std::vector<double> q2(numat + 1, 0.0);
        chrge(p, q2);
        double sumq = 0.0;
        for (int i = 1; i <= numat; ++i) {
            int l = nat[i];
            q[i] = parameters_C::tore[l] - q2[i];
            sumq += q[i];
            if (!l_atom[i]) continue;
            if (prt_charges)
                std::printf("%5d       %-2s  %15.6f%14.4f\n", i,
                            elemts_C::elemnt[l].c_str(), q[i], q2[i]);
        }
        (void)sumq;  // total charge (kchrge) not printed here
    }
    // DIPOLE section (F90: dipole(p, coord, dumy, 1) after NET ATOMIC CHARGES).
    if (id == 0) {
        std::vector<double> dipvec(4, 0.0);
        fprintf(stderr, "[DIPIN] coord[1..3][1]=%.6f %.6f %.6f  [1..3][2]=%.6f %.6f %.6f\n",
                coord[1][1], coord[2][1], coord[3][1], coord[1][2], coord[2][2], coord[3][2]); fflush(stderr);
        dipole(p, coord, dipvec, 1);
    }
    // FOCK / DENS matrices (F90 L683-698).
    if (norbs > 0) {
        if (keywrd.find(" FOCK") != std::string::npos) {
            std::printf(" FOCK MATRIX \n");
            vecprt(&f[0], norbs);
        }
        if (nelecs != 0 && keywrd.find(" DENS") != std::string::npos) {
            std::printf("\n\n%20s DENSITY MATRIX IS \n", "");
            vecprt(&p[0], norbs);
        }
    }
    // ATOMIC ORBITAL ELECTRON POPULATIONS.
    if (nelecs != 0 && prt_pops) {
        std::printf("\n\n%10sATOMIC ORBITAL ELECTRON POPULATIONS\n\n", "");
        int maxorb = 0;
        for (int i = 1; i <= numat; ++i)
            maxorb = std::max(maxorb, nlast[i] - nfirst[i]);
        if (maxorb == 8)
            std::printf("   Atom     s        px        py        pz      x^2-y^2     xz        z^2       yz        xy\n");
        else if (maxorb == 3)
            std::printf("   Atom     s        px        py        pz\n");
        else
            std::printf("   Atom     s\n");
        for (int i = 1; i <= numat; ++i) {
            std::printf("%5d %-2s", i, elemts_C::elemnt[nat[i]].c_str());
            for (int j = nfirst[i]; j <= nlast[i]; ++j)
                std::printf("%10.5f", p[(j * (j + 1)) / 2]);
            std::printf("\n");
        }
    }
    // PI bond orders (F90 L758-761).
    if (keywrd.find(" PI") != std::string::npos) {
        std::printf("\n\n%10sSIGMA-PI BOND-ORDER MATRIX\n", "");
        denrot();
    }
    // UHF spin analysis: SZ and <S**2>.
    if (nelecs != 0 && uhf) {
        double na_u = nalpha + (nalpha_open - nalpha) * fract;
        double nb_u = nbeta + (nbeta_open - nbeta) * fract;
        sz = (na_u - nb_u) * 0.5;
        ss2 = sz * sz;
        int l = 0;
        for (int i = 1; i <= norbs; ++i) {
            for (int j = 1; j <= i; ++j) {
                ++l;
                pa[l] = pa[l] - pb[l];
                ss2 = ss2 + pa[l] * pa[l];
            }
            ss2 = ss2 - 0.5 * pa[l] * pa[l];
        }
        std::printf("\n\n%20s(SZ)    =%12.6f\n", "", sz);
        if (fract < 1e-5 || fract > 0.9999)
            std::printf("%20s(S**2)  =%12.6f\n", "", ss2);
        else
            std::printf("%10sAverage over configurations used, so (S**2) is not meaningful\n", "");
        // SPIN density matrix (F90 L781-784).
        if (keywrd.find(" SPIN") != std::string::npos) {
            std::printf("\n\n%10sSPIN DENSITY MATRIX\n", "");
            vecprt(&pa[0], norbs);
        }

        // AO spin populations.
        std::printf("\n\n%10sATOMIC ORBITAL SPIN POPULATIONS\n\n", "");
        int maxorb = 0;
        for (int i = 1; i <= numat; ++i)
            maxorb = std::max(maxorb, nlast[i] - nfirst[i]);
        if (maxorb == 8)
            std::printf("    Atom   Total    s     px     py     pz   x^2-y^2  xz     z^2    yz     xy\n");
        else if (maxorb == 3)
            std::printf("    Atom   Total    s     px     py     pz\n");
        else
            std::printf("    Atom   Total    s\n");
        for (int i = 1; i <= numat; ++i) {
            double sum = 0.0;
            for (int j = nfirst[i]; j <= nlast[i]; ++j)
                sum += pa[(j * (j + 1)) / 2];
            std::printf("%5d %-2s%8.4f", i, elemts_C::elemnt[nat[i]].c_str(), sum);
            for (int j = nfirst[i]; j <= nlast[i]; ++j)
                std::printf("%7.3f", pa[(j * (j + 1)) / 2]);
            std::printf("\n");
        }
        for (size_t vi = 0; vi < p.size(); ++vi) pa[vi] = p[vi] - pb[vi];
    }
    // BONDS / ALLBO (F90 L819-832; molval M.O.-contribution part not ported).
    if (keywrd.find(" BONDS") != std::string::npos ||
        keywrd.find(" ALLBO") != std::string::npos) {
        bonds();
    }
    // 1ELE and ENPART (F90 L842-846).
    if (keywrd.find(" 1ELE") != std::string::npos) {
        std::printf(" FINAL ONE-ELECTRON MATRIX \n");
        vecprt(&h[0], norbs);
    }
    if (keywrd.find(" ENPART") != std::string::npos) enpart();
    // MULLIK / GRAPH (F90 L857-887; MOZYME lmo_to_eigenvectors part not ported).
    if ((keywrd.find(" MULLIK") != std::string::npos ||
         keywrd.find(" GRAPH") != std::string::npos) && !gui) {
        if (keywrd.find(" MULLIK") != std::string::npos)
            std::printf("\n%10s MULLIKEN POPULATION ANALYSIS\n", "");
        mullik();
        if (keywrd.find(" GRAPH") != std::string::npos)
            std::printf("\n%10s DATA FOR GRAPH WRITTEN TO DISK\n", "");
    }
    // TODO(translate): ~2500 lines remain: geometry/gradient tables,
    // eigenvalues & MO coefficients, Mulliken charges, COSMO report, MOZYME
    // occupied/virtual orbital blocks, and reference printing. These need
    // pa/pb density matrices, cosmo area/fepsi/cosvol, and MOZYME cocc/cvir.
    ++numcal;
}
