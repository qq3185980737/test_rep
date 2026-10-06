// switch.cpp — C++ translation of MOPAC 2016 "switch.F90".
// Loads the active Hamiltonian's reference parameter set (MNDO, AM1, PM3,
// PM6, PM7, PM7-TS, MNDO/d, RM1) from the parameters_for_*_C modules into
// the runtime arrays of parameters_C, applies SPARKLES overrides, fills the
// alpb/xfac tables and applies fractional-metal-ion fixes.
#include "switch.h"

#include "alpb_and_xfac.h"
#include "chanel_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "parameters_C.h"
#include "parameters_for_AM1_C.h"
#include "parameters_for_AM1_Sparkles_C.h"
#include "parameters_for_mndo_C.h"
#include "parameters_for_mndod_C.h"
#include "parameters_for_PM3_C.h"
#include "parameters_for_PM3_Sparkles_C.h"
#include "parameters_for_PM6_C.h"
#include "parameters_for_PM6_Sparkles_C.h"
#include "parameters_for_PM7_C.h"
#include "parameters_for_PM7_Sparkles_C.h"
#include "parameters_for_PM7_TS_C.h"
#include "parameters_for_RM1_C.h"
#include "parameters_for_RM1_Sparkles_C.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace parameters_C;
using namespace molkst_C;

// ---------------------------------------------------------------- helpers
namespace {

// Fortran whole-array copies become element-wise loops (dst is a C array,
// src is a 1-based std::vector).
inline void copy1D(double* dst, const std::vector<double>& src) {
    for (int i = 1; i <= 107; ++i) dst[i] = src[i];
}
inline void copy2D(double dst[][5], const std::vector<std::vector<double>>& src) {
    for (int i = 1; i <= 107; ++i)
        for (int j = 1; j <= 4; ++j) dst[i][j] = src[i][j];
}
inline void zero1D(double* dst) {
    for (int i = 1; i <= 107; ++i) dst[i] = 0.0;
}
inline void zero2D(double dst[][5]) {
    for (int i = 1; i <= 107; ++i)
        for (int j = 1; j <= 4; ++j) dst[i][j] = 0.0;
}

// Fortran: alpb(:100,57:71) = 0 ; alpb(57:71,:100) = 0 (and same for xfac)
void zero_land_alpb_xfac() {
    for (int i = 1; i <= 100; ++i)
        for (int j = 57; j <= 71; ++j) { alpb[i][j] = 0.0; xfac[i][j] = 0.0; }
    for (int i = 57; i <= 71; ++i)
        for (int j = 1; j <= 100; ++j) { alpb[i][j] = 0.0; xfac[i][j] = 0.0; }
}

// SPARKLES override (alp > 0.1 marks a sparkle ion); clears guess*(:,3:4).
void sparkle_overwrite(int limit, const std::vector<double>& alpSp,
                       const std::vector<double>& gssSp,
                       const std::vector<std::vector<double>>& guesSp1,
                       const std::vector<std::vector<double>>& guesSp2,
                       const std::vector<std::vector<double>>& guesSp3) {
    for (int i = 1; i <= limit; ++i) {
        if (alpSp[i] <= 0.1) continue;
        zd[i] = 0.0; zp[i] = 0.0; zs[i] = 0.0;
        zsn[i] = 0.0; zpn[i] = 0.0; zdn[i] = 0.0;
        uss[i] = 0.0; upp[i] = 0.0; udd[i] = 0.0;
        betas[i] = 0.0; betap[i] = 0.0; betad[i] = 0.0;
        alp[i] = alpSp[i];
        gss[i] = gssSp[i];
        gpp[i] = 0.0; gp2[i] = 0.0; hsp[i] = 0.0; gsp[i] = 0.0;
        f0sd[i] = 0.0; g2sd[i] = 0.0;
        pocord[i] = 0.0;
        guess1[i][1] = guesSp1[i][1];
        guess2[i][1] = guesSp2[i][1];
        guess3[i][1] = guesSp3[i][1];
        guess1[i][2] = guesSp1[i][2];
        guess2[i][2] = guesSp2[i][2];
        guess3[i][2] = guesSp3[i][2];
        guess1[i][3] = guess1[i][4] = 0.0;
        guess2[i][3] = guess2[i][4] = 0.0;
        guess3[i][3] = guess3[i][4] = 0.0;
    }
}

void fractional_metal_ion() {
    tore[85] = -0.5;   // pseudo-halide metal ion (core charge -1/2)
    tore[87] = 0.5;    // pseudo-alkali metal ion (core charge +1/2)
    for (int i = 85; i <= 87; i += 2) {
        upp[i] = 0.0;
        alp[i] = 3.0;
        zs[i] = 0.0;
        zp[i] = 0.0;
        betas[i] = 0.0;
        betap[i] = 0.0;
        gss[i] = 10.0;
        gsp[i] = 0.0;
        gpp[i] = 0.0;
        gp2[i] = 0.0;
        hsp[i] = 0.0;
    }
}

}  // namespace

// ---------------------------------------------------------------- switch
void switch_method() {
    using namespace parameters_for_mndod_C;
    using namespace parameters_for_mndo_C;
    using namespace parameters_for_PM6_C;
    using namespace parameters_for_PM7_C;
    using namespace parameters_for_PM7_TS_C;
    using namespace parameters_for_PM3_C;
    using namespace parameters_for_AM1_C;
    using namespace parameters_for_RM1_C;
    // Sparkles sets are referenced by fully-qualified names below because
    // parameters_C declares C-array versions of the same gss*sp symbols.

    zero1D(pocord);

    if (method_mndo) {
        copy2D(guess1, guesm1); copy2D(guess2, guesm2); copy2D(guess3, guesm3);
        copy1D(polvol, polvolm);
        copy1D(zs, zsm); copy1D(zp, zpm); copy1D(zd, zdm);
        copy1D(zsn, zsnm); copy1D(zpn, zpnm); copy1D(zdn, zdnm);
        copy1D(uss, ussm); copy1D(upp, uppm); copy1D(udd, uddm);
        copy1D(betas, betasm); copy1D(betap, betapm); copy1D(betad, betadm);
        copy1D(alp, alpm);
        copy1D(gss, gssm); copy1D(gpp, gppm); copy1D(gsp, gspm);
        copy1D(gp2, gp2m); copy1D(hsp, hspm);
        copy1D(f0sd, f0sdm); copy1D(g2sd, g2sdm);
        copy1D(pocord, pocm);
        alpb_and_xfac_mndo();
    } else if (method_pm3) {
        copy2D(guess1, guesp1); copy2D(guess2, guesp2); copy2D(guess3, guesp3);
        copy1D(polvol, polvolpm3);
        copy1D(zs, zspm3); copy1D(zp, zppm3); zero1D(zd);
        copy1D(uss, usspm3); copy1D(upp, upppm3);
        copy1D(betas, betasp); copy1D(betap, betapp);
        copy1D(alp, alppm3);
        copy1D(gss, gsspm3); copy1D(gpp, gpppm3); copy1D(gsp, gsppm3);
        copy1D(gp2, gp2pm3); copy1D(hsp, hsppm3);
        alpb_and_xfac_pm3();
        if (keywrd.find(" SPARK") != std::string::npos) {
            sparkle_overwrite(102,
                              parameters_for_PM3_Sparkles_C::alpPM3sp,
                              parameters_for_PM3_Sparkles_C::gssPM3sp,
                              parameters_for_PM3_Sparkles_C::guesPM3sp1,
                              parameters_for_PM3_Sparkles_C::guesPM3sp2,
                              parameters_for_PM3_Sparkles_C::guesPM3sp3);
            zero_land_alpb_xfac();
        }
    } else if (method_pm6) {
        copy2D(guess1, gues61); copy2D(guess2, gues62); copy2D(guess3, gues63);
        copy1D(zs, zs6); copy1D(zp, zp6); copy1D(zd, zd6);
        copy1D(zsn, zsn6); copy1D(zpn, zpn6); copy1D(zdn, zdn6);
        copy1D(uss, uss6); copy1D(upp, upp6); copy1D(udd, udd6);
        copy1D(betas, betas6); copy1D(betap, betap6); copy1D(betad, betad6);
        copy1D(gss, gss6); copy1D(gpp, gpp6); copy1D(gsp, gsp6);
        copy1D(gp2, gp26); copy1D(hsp, hsp6);
        copy1D(f0sd, f0sd6); copy1D(g2sd, g2sd6);
        copy1D(alp, alp6);
        copy1D(pocord, poc_6);
        copy1D(polvol, polvo6);
        alpb_and_xfac_pm6();
        if (keywrd.find(" SPARK") != std::string::npos) {
            sparkle_overwrite(107,
                              parameters_for_PM6_Sparkles_C::alp6sp,
                              parameters_for_PM6_Sparkles_C::gss6sp,
                              parameters_for_PM6_Sparkles_C::gues6sp1,
                              parameters_for_PM6_Sparkles_C::gues6sp2,
                              parameters_for_PM6_Sparkles_C::gues6sp3);
            zero_land_alpb_xfac();
        }
    } else if (method_pm7_ts) {
        copy2D(guess1, gues7_TS1); copy2D(guess2, gues7_TS2); copy2D(guess3, gues7_TS3);
        copy1D(zs, zs7_TS); copy1D(zp, zp7_TS); copy1D(zd, zd7_TS);
        copy1D(zsn, zsn7_TS); copy1D(zpn, zpn7_TS); copy1D(zdn, zdn7_TS);
        copy1D(uss, uss7_TS); copy1D(upp, upp7_TS); copy1D(udd, udd7_TS);
        copy1D(betas, betas7_TS); copy1D(betap, betap7_TS); copy1D(betad, betad7_TS);
        copy1D(gss, gss7_TS); copy1D(gpp, gpp7_TS); copy1D(gsp, gsp7_TS);
        copy1D(gp2, gp27_TS); copy1D(hsp, hsp7_TS);
        copy1D(f0sd, f0sd7_TS); copy1D(g2sd, g2sd7_TS);
        copy1D(alp, alp7_TS);
        copy1D(pocord, poc_7_TS);
        copy1D(polvol, polvo7_TS);
        alpb_and_xfac_pm7_TS();
    } else if (method_pm7) {
        copy2D(guess1, gues71); copy2D(guess2, gues72); copy2D(guess3, gues73);
        copy1D(zs, zs7); copy1D(zp, zp7); copy1D(zd, zd7);
        copy1D(zsn, zsn7); copy1D(zpn, zpn7); copy1D(zdn, zdn7);
        copy1D(uss, uss7); copy1D(upp, upp7); copy1D(udd, udd7);
        copy1D(betas, betas7); copy1D(betap, betap7); copy1D(betad, betad7);
        copy1D(gss, gss7); copy1D(gpp, gpp7); copy1D(gsp, gsp7);
        copy1D(gp2, gp27); copy1D(hsp, hsp7);
        copy1D(f0sd, f0sd7); copy1D(g2sd, g2sd7);
        copy1D(alp, alp7);
        copy1D(pocord, poc_7);
        copy1D(polvol, polvo7);
        alpb_and_xfac_pm7();
        if (keywrd.find(" SPARK") != std::string::npos) {
            sparkle_overwrite(102,
                              parameters_for_PM7_Sparkles_C::alp7sp,
                              parameters_for_PM7_Sparkles_C::gss7sp,
                              parameters_for_PM7_Sparkles_C::gues7sp1,
                              parameters_for_PM7_Sparkles_C::gues7sp2,
                              parameters_for_PM7_Sparkles_C::gues7sp3);
            zero_land_alpb_xfac();
        }
    } else if (method_mndod) {
        copy1D(uss, ussd); copy1D(upp, uppd); copy1D(udd, uddd);
        copy1D(zs, zsd); copy1D(zp, zpd); copy1D(zd, zdd);
        copy1D(zsn, zsnd); copy1D(zpn, zpnd); copy1D(zdn, zdnd);
        copy1D(betas, betasd); copy1D(betap, betapd); copy1D(betad, betadd);
        copy1D(gss, gssd); copy1D(gsp, gspd); copy1D(gpp, gppd);
        copy1D(gp2, gp2d); copy1D(hsp, hspd);
        copy1D(alp, alpd);
        copy1D(pocord, poc_d);
        zero2D(guess1); zero2D(guess2); zero2D(guess3);
        zero1D(f0sd); zero1D(g2sd);
        alpb_and_xfac_mndod();
    } else if (method_rm1) {
        copy2D(guess1, guess1RM1); copy2D(guess2, guess2RM1); copy2D(guess3, guess3RM1);
        copy1D(uss, ussRM1); copy1D(upp, uppRM1);
        copy1D(zs, zsRM1); copy1D(zp, zpRM1);
        copy1D(betas, betasRM1); copy1D(betap, betapRM1);
        copy1D(gss, gssRM1); copy1D(gsp, gspRM1); copy1D(gpp, gppRM1);
        copy1D(gp2, gp2RM1); copy1D(hsp, hspRM1);
        copy1D(alp, alpRM1);
        if (keywrd.find(" SPARK") != std::string::npos) {
            sparkle_overwrite(102,
                              parameters_for_RM1_Sparkles_C::alprm1sp,
                              parameters_for_RM1_Sparkles_C::gssrm1sp,
                              parameters_for_RM1_Sparkles_C::guesrm1sp1,
                              parameters_for_RM1_Sparkles_C::guesrm1sp2,
                              parameters_for_RM1_Sparkles_C::guesrm1sp3);
        }
    } else {
        // AM1 (default)
        copy2D(guess1, guesa1); copy2D(guess2, guesa2); copy2D(guess3, guesa3);
        copy1D(polvol, polvolam1);
        copy1D(zs, zsam1); copy1D(zp, zpam1); copy1D(zd, zdam1);
        copy1D(zsn, zsnam1); copy1D(zpn, zpnam1); copy1D(zdn, zdnam1);
        copy1D(f0sd, f0sdam1); copy1D(g2sd, g2sdam1);
        copy1D(uss, ussam1); copy1D(upp, uppam1); copy1D(udd, uddam1);
        copy1D(betas, betasa); copy1D(betap, betapa); copy1D(betad, betada);
        copy1D(alp, alpam1);
        copy1D(gss, gssam1); copy1D(gpp, gppam1); copy1D(gsp, gspam1);
        copy1D(gp2, gp2am1); copy1D(hsp, hspam1);
        alpb_and_xfac_am1();
        // Unique parameter for Voityuk's Molybdenum
        pocord[42] = 1.334;
        if (keywrd.find(" SPARK") != std::string::npos) {
            sparkle_overwrite(102,
                              parameters_for_AM1_Sparkles_C::alpam1sp,
                              parameters_for_AM1_Sparkles_C::gssam1sp,
                              parameters_for_AM1_Sparkles_C::guesam1sp1,
                              parameters_for_AM1_Sparkles_C::guesam1sp2,
                              parameters_for_AM1_Sparkles_C::guesam1sp3);
            zero_land_alpb_xfac();
        }
    }

    // Back up f0sd/g2sd (used by PARAM, not by MOPAC)
    for (int i = 0; i <= 107; ++i) { f0sd_store[i] = f0sd[i]; g2sd_store[i] = g2sd[i]; }

    // Symmetrize alpb and xfac (lower triangle mirrors upper triangle)
    for (int i = 1; i <= 100; ++i)
        for (int j = 1; j <= i; ++j) {
            alpb[j][i] = alpb[i][j];
            xfac[j][i] = xfac[i][j];
        }

    fractional_metal_ion();

    if (keywrd.find(" EXTERNAL") != std::string::npos) return;
    if (uss[1] > -1.0) {
        std::printf("  THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM\n");
        mopend("THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM");
        return;
    }
}
