// test_outer1.cpp
#include <cstdio>
#include "outer1.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
namespace molkst_C { bool method_pm7=false; int l1u,l2u,l3u; }
namespace funcon_C { double a0=0.529, ev=27.211; }
extern "C" double trunk_(double* r) { return *r; }
extern "C" void to_point_(double, double& p, double& c) { p=1.0; c=0.0; }
int main() {
    using namespace parameters_C;
    tore[1]=1.0; tore[2]=1.0; am[1]=1.0; am[2]=1.0;
    natorb[1]=1; natorb[2]=1;
    double c1[3]={0,0,0}, c2[3]={1.0,0,0}, w[5], e1[45], e2[45], en;
    int kr=0;
    outer1(1,2,c1,c2,w,kr,e1,e2,en,0,false);
    std::printf("outer1 en=%g kr=%d PASS\n", en, kr);
    return 0;
}
