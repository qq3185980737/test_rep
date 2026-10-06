// refer.cpp — C++ translation of "refer.F90".
// Prints journal references for elements used, based on method.

#include "refer.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "journal_references_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

void mopend(const std::string& s);

static std::string& allref(int i, int mode) {
    using namespace journal_references_C;
    switch (mode) {
        case 1: return refmn[i];
        case 2: return refpm6[i];
        case 3: return refam[i];
        case 4: return refpm3[i];
        case 5: return refmd[i];
        case 6: return refrm1[i];
        case 7: return refpm7[i];
        default: return refmn[i];
    }
}

void refer() {
    using chanel_C::iw;
    std::ostream& out = std::cout;  // Fortran iw unit; redirect via chanel when I/O wired.

    bool elemns[108] = {};
    for (int i = 1; i <= molkst_C::numat; ++i)
        elemns[common_arrays_C::nat[i]] = true;

    static bool mix = false;

    if (molkst_C::method_pm6) {
        bool allok = true;
        if (molkst_C::sparkle) {
            out << "\n References for PM6:\n";
            for (int i = 57; i <= 71; ++i)
                if (elemns[i]) out << allref(i, 2) << "\n";
        } else {
            for (int i = 58; i <= 70; ++i)
                if (elemns[i]) {
                    out << " DATA ARE NOT AVAILABLE FOR ELEMENT NO." << i << "\n";
                    if (std::abs(parameters_C::gss6sp[i]) > 0.1)
                        out << " (Parameters are available if PM6 SPARKLE is used)\n";
                    allok = false;
                }
        }
        out << "\n General Reference for PM6:\n";
        out << " \"Optimization of Parameters for Semiempirical Methods V: Modification of NDDO Approximations\",\n";
        out << " and Application to 70 Elements\", J. J. P. Stewart, J. Mol. Mod., 13, 1173-1213 (2007)\n";
        out << " URL: http://www.springerlink.com/content/ar33482301010477/fulltext.pdf\n";
        if (molkst_C::keywrd.find("PM6-DH+") != std::string::npos) {
            out << "\n Reference for PM6-DH+:\n";
            out << " \"Third-Generation Hydrogen-Bonding Corrections for Semiempirical QM Methods and Force Fields\"\n";
            out << " Martin Korth, J. Chem. Theory Comput., 6 (12), pp 3808-3816 (2010)\n";
            out << " URL: http://pubs.acs.org/doi/abs/10.1021/ct100408b\n";
        }
        if (allok || molkst_C::is_PARAM) return;
        out << "\n\n\n\n\n          SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH\n";
        out << "          NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.\n";
        mopend("Parameters for some elements are missing");
        return;
    }

    if (molkst_C::method_pm7) {
        bool allok = true;
        if (molkst_C::sparkle) {
            out << "\n References for PM7 Sparkles:\n";
            for (int i = 57; i <= 71; ++i)
                if (elemns[i]) {
                    out << " \"Sparkle/PM7 Lanthanide Parameters for the Modeling of Complexes and Materials\",\n";
                    out << " J. D. L. Dutra, M. A. M. Filho, R. O. Freire, G. B. Rocha, A. M. Simas, and J. J. P. Stewart,\n";
                    out << " Chemistry of Materials (submitted)\n";
                    break;
                }
        } else {
            for (int i = 58; i <= 70; ++i)
                if (elemns[i]) {
                    out << " DATA ARE NOT AVAILABLE FOR ELEMENT NO." << i << "\n";
                    if (std::abs(parameters_C::gss6sp[i]) > 0.1)
                        out << " (Parameters are available if SPARKLE is used)\n";
                    allok = false;
                }
        }
        out << "\n General Reference for PM7:\n";
        out << " \"Optimization of Parameters for Semiempirical Methods VI: More Modifications to the \n";
        out << " NDDO Approximations and Re-optimization of Parameters\", J. J. P. Stewart, J. Mol. Mod., 1:32, 19 (2013)\n";
        out << " http://www.springerlink.com/openurl.asp?genre=article&id=doi:10.1007/s00894-012-1667-x\n";
        if (allok || molkst_C::is_PARAM) return;
        out << "\n\n\n\n\n          SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH\n";
        out << "          NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.\n";
        mopend("Parameters for some elements are missing");
        return;
    }

    bool mixok = molkst_C::keywrd.find("PARASOK") != std::string::npos;
    bool exter = molkst_C::keywrd.find("EXTERNAL") != std::string::npos;

    int mode = 0;
    if (molkst_C::method_pm7) mode = 7;
    if (molkst_C::method_rm1) mode = 6;
    if (molkst_C::method_mndod) mode = 5;
    if (molkst_C::method_pm3) mode = 4;
    if (molkst_C::method_am1) mode = 3;
    if (molkst_C::method_pm6) mode = 2;
    if (molkst_C::method_mndo) mode = 1;

    allref(99, mode) = " DUMMY ATOMS ARE USED; THESE DO NOT AFFECT THE CALCULATION";
    allref(100, mode) = " ";
    bool allok = true;
    out << "\n";
    for (int i = 1; i <= 102; ++i) {
        if (!elemns[i]) continue;
        if (i < 99 && !mix && mode == 3)
            mix = allref(i, 3).find("MNDO") != std::string::npos;
        std::string& r5 = allref(i, 5);
        if (r5.empty() || r5.substr(0, 4) == "    ") allref(i, 5) = allref(i, 1);
        std::string& rm = allref(i, mode);
        if (rm.empty() || rm.substr(0, 4) == "    " || std::abs(parameters_C::gss[i]) < 0.1) {
            if (!exter) {
                out << " DATA ARE NOT AVAILABLE FOR ELEMENT NO." << i << "\n";
                if (std::abs(parameters_C::gssam1sp[i]) > 0.1)
                    out << " (Parameters are available if AM1 SPARKLE is used)\n";
                if (std::abs(parameters_C::gssPM3sp[i]) > 0.1)
                    out << " (Parameters are available if PM3 SPARKLE is used)\n";
                if (std::abs(parameters_C::gss6sp[i]) > 0.1)
                    out << " (Parameters are available if PM6 SPARKLE is used)\n";
                allok = false;
            }
        } else {
            out << rm << "\n";
            if (mode == 7) break;
        }
    }
    if (mix && !mixok) {
        out << "\n\n\n\n\n          SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH ONLY MNDO\n";
        out << "          PARAMETERS ARE AVAILABLE.  SUCH MIXTURES OF METHODS ARE\n";
        out << "          VERY RISKY AND HAVE NOT BEEN FULLY TESTED.  IF YOU FEEL\n";
        out << "          THE RISK IS WORTH WHILE - CHECK THE MANUAL FIRST - THEN\n";
        out << "          SPECIFY \"PARASOK\" IN THE KEYWORDS\n";
        mopend("MIXED PARAMETER SETS.  USE \"PARASOK\" TO CONTINUE");
        return;
    }
    if (allok || molkst_C::is_PARAM ||
        molkst_C::keywrd.find("0SCF") != std::string::npos) return;
    out << "\n\n\n\n\n          SOME ELEMENTS HAVE BEEN SPECIFIED FOR WHICH\n";
    out << "          NO PARAMETERS ARE AVAILABLE.  CALCULATION STOPPED.\n";
    mopend("Parameters for some elements are missing");
}
