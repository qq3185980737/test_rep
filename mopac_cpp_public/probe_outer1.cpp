// probe_outer1.cpp — debug outer1 mode0/mode1 values
#include <cstdio>
#include <cmath>
#include "outer1.h"
#include "molkst_C.h"
#include "funcon_C.h"
#include "parameters_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace funcon_C;
using namespace parameters_C;
using namespace common_arrays_C;
int main() {
    am[1] = 100.0; tore[1] = 1.0; natorb[1] = 1; method_pm7 = false;
    double w[2026] = {0.0}, e1b[45] = {0.0}, e2a[45] = {0.0};
    double c1[3] = {0,0,0}, c2[3] = {0,0,a0};
    int kr = 0; double enuc = 0;
    outer1(1,1,c1,c2,w,kr,e1b,e2a,enuc,0,false);
    double ww = ev/std::sqrt(1.0+1e-4);
    std::printf("mode0: w[0]=%.12f exp=%.12f e1b[0]=%.12f e2a[0]=%.12f enuc=%.12f kr=%d\n", w[0], ww, e1b[0], e2a[0], enuc, kr);
    int idx9[9] = {0,2,5,9,14,20,27,35,44};
    for (int q = 0; q < 9; ++q)
        std::printf("  e1b[%d]=%.12f e2a[%d]=%.12f\n", idx9[q], e1b[idx9[q]], idx9[q], e2a[idx9[q]]);
    // mode=1
    l1u=l2u=l3u=0; cutofp=10.0; numcal=1; clower=cupper=0;
    for (int i=0;i<4;++i) tvec[i].resize(4,0.0);
    tvec[1][1]=5.0; tvec[2][2]=5.0; tvec[3][3]=5.0;
    double w2[2026]={0.0}, e1b2[45]={0.0}, e2a2[45]={0.0};
    double c2b[3]={1.0,0,0}; kr=0; enuc=0;
    outer1(1,1,c1,c2b,w2,kr,e1b2,e2a2,enuc,1,false);
    std::printf("mode1: w[0]=%.12f exp=%.12f (a0*ev=%.12f) e1b[0]=%.12f e1b[25]=%.12f enuc=%.12f kr=%d clower=%.12f\n", w2[0], a0*ev, a0*ev, e1b2[0], e1b2[25], enuc, kr, clower);
    double bound1=clower/cutofp, range=1.0-bound1;
    double cT=-0.5*bound1*bound1*cutofp/range, crT=1.0+bound1/range, cr2T=-1.0/(cutofp*2.0*range);
    double clim=cT+crT*cupper+cr2T*cupper*cupper;
    std::printf("clim=%.12f bound1=%.12f range=%.12f\n", clim, bound1, range);
    double w3[2026]={0.0}, e1b3[45]={0.0}, e2a3[45]={0.0};
    double c2c[3]={20.0,0,0}; kr=0; enuc=0;
    outer1(1,1,c1,c2c,w3,kr,e1b3,e2a3,enuc,1,false);
    std::printf("mode1-far: w3[0]=%.12f exp=%.12f kr=%d\n", w3[0], a0*ev/clim, kr);
    return 0;
}
