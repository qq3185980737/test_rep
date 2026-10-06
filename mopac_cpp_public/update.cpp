// update.cpp
#include "update.h"
#include <cmath>
#include "parameters_C.h"
using namespace parameters_C;
extern "C" int iw_;
void update(int iparam, int ielmnt, double param, double c1) {
    int jparam = iparam;
    if (jparam > 21 && jparam < 34) {
        int kfn = (jparam-22)/3;
        jparam = jparam - kfn*3;
        kfn = kfn + 1;
        switch (jparam) {
            case 22: guess1[ielmnt][kfn] = guess1[ielmnt][kfn]*c1 + param; return;
            case 23: guess2[ielmnt][kfn] = guess2[ielmnt][kfn]*c1 + param; return;
            case 24: guess3[ielmnt][kfn] = guess3[ielmnt][kfn]*c1 + param; return;
        }
    }
    switch (jparam) {
        case 2: upp[ielmnt] = upp[ielmnt]*c1 + param; break;
        case 3: udd[ielmnt] = udd[ielmnt]*c1 + param; break;
        case 4: zs[ielmnt] = zs[ielmnt]*c1 + param; break;
        case 5: zp[ielmnt] = zp[ielmnt]*c1 + param; break;
        case 6: zd[ielmnt] = zd[ielmnt]*c1 + param; break;
        case 7: betas[ielmnt] = betas[ielmnt]*c1 + param; break;
        case 8: betap[ielmnt] = betap[ielmnt]*c1 + param; break;
        case 9: betad[ielmnt] = betad[ielmnt]*c1 + param; break;
        case 10: gss[ielmnt] = gss[ielmnt]*c1 + param; break;
        case 11: gsp[ielmnt] = gsp[ielmnt]*c1 + param; break;
        case 12: gpp[ielmnt] = gpp[ielmnt]*c1 + param; break;
        case 13: gp2[ielmnt] = gp2[ielmnt]*c1 + param; break;
        case 14: hsp[ielmnt] = hsp[ielmnt]*c1 + param; break;
        case 15:
            f0sd_store[ielmnt] = f0sd_store[ielmnt]*c1 + param;
            if (c1 < 1e-20) f0sd[ielmnt] = param;
            break;
        case 16:
            g2sd_store[ielmnt] = g2sd_store[ielmnt]*c1 + param;
            if (c1 < 1e-20) g2sd[ielmnt] = param;
            break;
        case 17: pocord[ielmnt] = pocord[ielmnt]*c1 + param; break;
        case 18: alp[ielmnt] = alp[ielmnt]*c1 + param; break;
        case 19: zsn[ielmnt] = zsn[ielmnt]*c1 + param; break;
        case 20: zpn[ielmnt] = zpn[ielmnt]*c1 + param; break;
        case 21: zdn[ielmnt] = zdn[ielmnt]*c1 + param; break;
        case 37:
            natorb[ielmnt] = (int)std::lround(param);
            dorbs[ielmnt] = (natorb[ielmnt] == 9);
            break;
        case 34: {
            int nj = ielmnt/200, ni = ielmnt - nj*200;
            alpb[ni][nj] = alpb[ni][nj]*c1 + param;
            alpb[nj][ni] = alpb[ni][nj];
            break;
        }
        case 35: {
            int nj = ielmnt/200, ni = ielmnt - nj*200;
            xfac[ni][nj] = xfac[ni][nj]*c1 + param;
            xfac[nj][ni] = xfac[ni][nj];
            break;
        }
        default: uss[ielmnt] = uss[ielmnt]*c1 + param; break;
    }
}
