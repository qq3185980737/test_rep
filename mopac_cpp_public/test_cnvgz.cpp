// test_cnvgz.cpp
#include <cstdio>
#include <vector>
#include "cnvgz.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
namespace molkst_C { int norbs, mpack; }
namespace MOZYME_C { bool use_three_point_extrap=false; double pmax; }
int main() {
    using namespace molkst_C;
    using namespace MOZYME_C;
    norbs=2; mpack=3;
    std::vector<double> pn={0,1.0,2.0,3.0}, p={0,0.5,1.0,1.5};
    std::vector<double> p1(3,0),p2(3,0),p3(3,0);
    std::vector<int> idiag={0,1,3};
    cnvgz(pn,p,p1,p2,p3,1,idiag);
    bool ok=(pmax>0);
    std::printf("pmax=%g %s\n",pmax,ok?"PASS":"FAIL");
    return ok?0:1;
}
