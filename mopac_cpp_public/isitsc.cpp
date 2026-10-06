// isitsc.cpp
#include "isitsc.h"
#include "MOZYME_C.h"
#include <cmath>
#include <algorithm>
using namespace MOZYME_C;
namespace molkst_C { extern int iscf; }
using namespace molkst_C;
void isitsc(double escf, double selcon, double emin, int& iemin, int& iemax, bool& okscf, int niter, int itrmax) {
    static bool scf1 = false;
    static double escf0[6] = {0,0,0,0,0,0};
    double energy_test = selcon;
    double fmo_test = selcon * 5.0;
    bool converged = (ovmax < fmo_test && std::abs(energy_diff) < energy_test && scf1) || (niter > itrmax);
    if (converged) {
        okscf = true;
        iscf = scf1 ? 1 : 2;
        return;
    }
    if (ovmax < fmo_test && std::abs(energy_diff) < energy_test) scf1 = true; else scf1 = false;
    if (emin == 0.0) { okscf = false; return; }
    if (escf < emin) {
        iemax = 0;
        int iemin1 = iemin;
        iemin = std::min(5, iemin + 1);
        if (iemin1 == 5) for (int i=2;i<=5;++i) escf0[i-1]=escf0[i];
        escf0[iemin] = escf;
        if (iemin > 3) {
            for (int i=2;i<=iemin;++i)
                if (std::abs(escf0[i]-escf0[i-1]) > 0.1*(emin-escf)) { okscf=false; return; }
            okscf = true; iscf = 1; return;
        }
    } else {
        iemin = 0;
        int iemax1 = iemax;
        iemax = std::min(5, iemax+1);
        if (iemax1 == 5) for (int i=2;i<=5;++i) escf0[i-1]=escf0[i];
        escf0[iemax] = escf;
        if (iemax > 3) {
            for (int i=2;i<=iemax;++i)
                if (std::abs(escf0[i]-escf0[i-1]) > 0.1*(escf-emin)) { okscf=false; return; }
            okscf = true; iscf = 1; return;
        }
    }
    okscf = false;
}
