// test_diagg2.cpp
#include <cstdio>
#include <vector>
#include "diagg2.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numcal=1; keywrd=""; numat=1;
    std::vector<double> eigv(2), si(2), sj(2);
    std::vector<int> iused(2);
    std::vector<char> latoms(2);
    diagg2(0,0,eigv,iused,latoms,0,1,si,sj);
    std::printf("diagg2(skeleton) PASS\n");
    return 0;
}
