// test_drc.cpp
#include <cstdio>
#include <vector>
#include "drc.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C; using namespace common_arrays_C;
    numat=1; numcal=1; keywrd="";
    atmass.assign(2,1.0);
    std::vector<double> startv(4,0), startk(1,0);
    startv[1]=3;
    drc(startv,startk);
    std::printf("drc(skeleton) PASS\n");
    return 0;
}
