// test_dijkl1.cpp
#include <cstdio>
#include <vector>
#include "dijkl1.h"
#include "common_arrays_C.h"
#include "meci_C.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C; using namespace meci_C;
    numat=1; nmos=1;
    common_arrays_C::nfirst.assign(2,1); common_arrays_C::nlast.assign(2,1);
    std::vector<std::vector<double>> c(2,std::vector<double>(2,1));
    std::vector<double> w(1), cij(1), wcij(1), ckl(1), xybuf(1);
    dijkl1(c,0,1,w,cij,wcij,ckl,xybuf);
    std::printf("dijkl1(skeleton) PASS\n");
    return 0;
}
