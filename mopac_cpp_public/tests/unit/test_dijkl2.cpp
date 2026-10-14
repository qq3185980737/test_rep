// test_dijkl2.cpp
#include <cstdio>
#include <vector>
#include "dijkl2.h"
#include "meci_C.h"
#include "molkst_C.h"
namespace meci_C {
    int nmos;
    std::vector<double> dijkl, xy;
}
namespace molkst_C { int norbs; }
int main() {
    using namespace meci_C; using namespace molkst_C;
    nmos=1; norbs=1;
    dijkl.assign(2,0.0); xy.assign(16,0.0);
    std::vector<std::vector<double>> dc(2,std::vector<double>(2,0));
    dijkl2(dc);
    std::printf("dijkl2 PASS\n");
    return 0;
}
