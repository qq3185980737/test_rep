// test_mecid.cpp
#include <cstdio>
#include <vector>
#include "mecid.h"
#include "meci_C.h"
namespace meci_C {
    int nmos, lab;
    std::vector<double> occa;
    std::vector<std::vector<int>> microa, microb;
}
extern "C" double diagi_(int*, int*, double*, double*, int*) { return 0.0; }
int main() {
    using namespace meci_C;
    nmos = 2; lab = 0;
    occa = {0, 1.0, 0.0};
    double eigs[3] = {0, -1.0, -2.0};
    double eiga[3], diag[3], gse;
    double xy[16] = {}; // identity: xy(i,i,j,j)=0 for minimal
    mecid(eigs, gse, eiga, diag, xy);
    // eiga(1) = -1 - 0 = -1; gse = 2*(-1)*1 = -2
    bool ok = (gse == -2.0);
    std::printf("gse=%g %s\n", gse, ok?"PASS":"FAIL");
    return ok?0:1;
}
