// sp_two_electron.cpp
#include "sp_two_electron.h"
#include "parameters_C.h"
#include "mndod_C.h"
using namespace parameters_C;
using mndod_C::iii;
extern double rsc(int, int, double, int, double, int, double, int, double);
void sp_two_electron() {
    for (int ni=1; ni<=80; ++ni) {
        int ns = iii[ni];
        double es = zsn[ni], ep = zpn[ni];
        if (es<1e-4 || ep<1e-4 || main_group[ni]) continue;
        gss[ni] = rsc(0, ns, es, ns, es, ns, es, ns, es);
        gsp[ni] = rsc(0, ns, es, ns, es, ns, ep, ns, ep);
        hsp[ni] = rsc(1, ns, es, ns, ep, ns, es, ns, ep)/3.0;
        double r033 = rsc(0, ns, ep, ns, ep, ns, ep, ns, ep);
        double r233 = rsc(2, ns, ep, ns, ep, ns, ep, ns, ep);
        gpp[ni] = r033 + 0.16*r233;
        gp2[ni] = r033 - 0.08*r233;
    }
}
