// wrtkey.cpp — C++ translation of "wrtkey.F90".
// Entry point: print recognized keyword groups and run conflict checks.
// wrtcon / wrtwor / wrtout are the three printer groups (large; stubbed).
// wrtchk: keyword-conflict validator (implemented from Fortran logic).

#include "wrtkey.h"

#include <algorithm>
#include <string>

#include "chanel_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "reada.h"

using namespace chanel_C;
using namespace molkst_C;

// Printer groups (large keyword printers) — declared elsewhere / stubbed.
static void wrtcon(std::string&) {}
static void wrtwor(std::string&) {}
static void wrtout(std::string&) {}

static bool has(const std::string& s, const std::string& w) {
    return s.find(w) != std::string::npos;
}

static void wrtchk(std::string& allkey) {
    (void)allkey;
    bool birad = has(keywrd, " BIRAD");
    bool exci = has(keywrd, " EXCI");
    bool ci = has(keywrd, " C.I.");
    bool trip = has(keywrd, " TRIP");
    uhf = has(keywrd, " UHF");
    rhf = has(keywrd, " RHF") || has(keywrd, " MECI") || ci ||
          (!uhf && has(keywrd, " OPEN"));

    // F90 wrtkey.F90:652-680: parse/validate the C.I. keyword.  A bare "C.I."
    // (no '=n' or '=(n1,n2)') is an error; otherwise set nmos exactly as the
    // Fortran does (nmos is consumed later by meci/densit).  Without this, an
    // "C.I.=n" job would run with nmos=0 and give a wrong CI energy.
    if (ci) {
        if (keywrd.find(" C.I.=(") != std::string::npos) {
            // C.I.=(n1,n2): nmos = Nint(reada(keywrd, Index("C.I.=(")+5))
            meci_C::nmos = (int)(reada(keywrd, (int)keywrd.find("C.I.=(") + 5) + 0.5);
        } else if (keywrd.find(" C.I.=") != std::string::npos) {
            // C.I.=n: nmos = Nint(reada(keywrd, Index("C.I.=")+5))
            meci_C::nmos = (int)(reada(keywrd, (int)keywrd.find("C.I.=") + 5) + 0.5);
        } else {
            mopend("C.I. keyword must be of form 'C.I.=n' or 'C.I.=(n1,n2)'");
            return;
        }
        if (meci_C::nmos > 29) {
            mopend("Maximum size of open space = 29 M.O.s");
            return;
        }
    }

    if (!(mozyme || has(keywrd, " PDBOUT") || has(keywrd, " RESID") ||
          has(keywrd, " ADD-H")) &&
        !has(keywrd, " 0SCF")) {
        if (has(keywrd, " PDBOUT")) {
            mopend("Keyword PDBOUT only works when MOZYME or 0SCF is also present");
            return;
        }
        if (has(keywrd, " RESID")) {
            mopend("Keyword RESIDUES only works when MOZYME or 0SCF is also present");
            return;
        }
        if (has(keywrd, " CVB")) {
            mopend("Keyword CVB only works with MOZYME");
            return;
        }
        if (has(keywrd, " SETPI")) {
            mopend("Keyword SETPI only works with MOZYME");
            return;
        }
    }
    if (mozyme) {
        if (uhf) mopend("Keyword UHF cannot be used with MOZYME");
        if (has(keywrd, " ENPART")) mopend("Keyword ENPART is not available with MOZYME");
        if (has(keywrd, " LOCAL")) mopend("Keyword LOCAL is not available with MOZYME");
        if (id != 0 && has(keywrd, " CUTOF")) {
            mopend("CUTOFx=n.nn -type keywords do not work with MOZYME for infinite systems");
            return;
        }
        if (has(keywrd, " 1SCF") && has(keywrd, " RAPID")) {
            mopend("RAPID cannot be used with 1SCF");
            return;
        }
    }
    if (has(keywrd, " MULLIK") && uhf) {
        mopend("MULLIKEN POPULATION NOT AVAILABLE WITH UHF");
        return;
    }
    if (uhf) {
        if (rhf) mopend("UHF and RHF cannot both be used");
        if (birad || exci || ci) {
            std::printf("UHF USED WITH EITHER BIRAD, EXCITED OR C.I.\n");
            mopend("IMPOSSIBLE OPTION REQUESTED");
            return;
        }
        if (has(keywrd, " POLAR")) {
            mopend("POLAR does not work with UHF");
            return;
        }
    } else if (exci && trip) {
        std::printf("EXCITED USED WITH TRIPLET\n");
        mopend("IMPOSSIBLE OPTION REQUESTED");
        return;
    }
    if (has(keywrd, " PMEP") && !method_am1) {
        mopend("PMEP only works with AM1");
        return;
    }
    if (has(keywrd, " INT ") && has(keywrd, " XYZ")) {
        mopend("INT cannot be used with XYZ");
        return;
    }
    if (id > 0 && has(keywrd, "EPS=")) {
        if (id == 1) std::printf("COSMO cannot be used with polymers\n");
        if (id == 2) std::printf("COSMO cannot be used with layer systems\n");
        if (id == 3) std::printf("COSMO cannot be used with solids\n");
        mopend("COSMO cannot be used with systems with Tv");
        return;
    }
    if (has(keywrd, " T-PRIO") && !has(keywrd, " DRC")) {
        std::printf("T-PRIO AND NO DRC\n");
        mopend("IMPOSSIBLE OPTION REQUESTED");
        return;
    }
    // Only one method allowed.
    int m = 0;
    if (method_am1) ++m;
    if (method_pm3) ++m;
    if (method_pm6) ++m;
    if (method_pm7) ++m;
    if (method_mndo) ++m;
    if (method_mndod) ++m;
    if (method_rm1) ++m;
    if (m > 1) {
        mopend("ONLY ONE OF MNDO, MNDOD, AM1, PM3, RM1, AND PM6 ALLOWED");
        return;
    }
    // Only one geometry option: strip quoted strings first.
    line = keywrd;
    for (size_t pos = 0; pos < line.size(); ++pos) {
        if (line[pos] == '"') {
            size_t e = line.find('"', pos + 1);
            if (e == std::string::npos) break;
            for (size_t k = pos; k <= e; ++k) line[k] = ' ';
        }
    }
    int g = 0;
    if (has(line, " BFGS")) ++g;
    if (has(line, " LBFGS")) ++g;
    if (has(line, " EF")) ++g;
    if (has(line, " TS")) ++g;
    if (has(line, " SIGMA")) ++g;
    if (has(line, " NLLSQ")) ++g;
    if (has(line, " FORCE") || has(line, " IRC") || has(line, " DRC")) ++g;
    if (g > 1) {
        mopend("MORE THAN ONE GEOMETRY OPTION HAS BEEN SPECIFIED. CONFLICT MUST BE RESOLVED BEFORE JOB WILL RUN.");
        return;
    }
    if (has(keywrd, " HESSIAN") && !has(keywrd, " EF")) {
        mopend("Keyword EF must be present if HESSIAN is used");
        return;
    }
}

void wrtkey() {
    wrtcon(allkey);
    if (moperr) return;
    wrtwor(allkey);
    wrtout(allkey);
    // Fortran wrtkey L921-958: OUTPUT keyword controls the prt_* flags.
    // Without OUTPUT all print flags default to TRUE.
    if (has(keywrd, " OUTPUT")) {
        size_t i = keywrd.find(" OUTPUT");
        size_t j = keywrd.find(") ", i);
        if (j != std::string::npos && keywrd.find("OUTPUT(") != std::string::npos) {
            size_t a = keywrd.find("OUTPUT(") + 6;
            std::string line_ = keywrd.substr(a, j - a);
            prt_coords    = line_.find("C") != std::string::npos;
            prt_gradients = line_.find("G") != std::string::npos;
            prt_cart      = line_.find("X") != std::string::npos;
            prt_charges   = line_.find("Q") != std::string::npos;
            prt_pops      = line_.find("P") != std::string::npos;
            prt_topo      = line_.find("T") != std::string::npos;
        } else {
            prt_coords = false;
            prt_gradients = false;
            prt_cart = false;
            prt_charges = false;
            prt_pops = false;
            prt_topo = false;
        }
    } else {
        prt_coords = true;
        prt_gradients = true;
        prt_cart = true;
        prt_charges = true;
        prt_pops = true;
        prt_topo = true;
    }
    wrtchk(allkey);
}
