// test_ijbo.cpp
#include <cstdio>
#include <vector>
#include "ijbo.h"
#include "overlaps_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
namespace overlaps_C {
    double cutof1 = 100.0, cutof2 = 10.0;
}
namespace common_arrays_C {
    std::vector<std::vector<double>> coord;
}
namespace MOZYME_C {
    bool lijbo = false;
    std::vector<std::vector<int>> nijbo;
    std::vector<int> iij, numij, ijall, iijj;
}
int main() {
    using namespace overlaps_C;
    using namespace common_arrays_C;
    using namespace MOZYME_C;
    coord.resize(3, std::vector<double>(10, 0.0));
    coord[0][1] = 0; coord[0][2] = 1.0;
    // set up lookup: atom 2 has partner 1 at index 5 with offset 0
    iij.assign(10, 0); numij.assign(10, 0);
    ijall.assign(10, 0); iijj.assign(10, 0);
    iij[2] = 1; numij[2] = 5;
    ijall[3] = 1; iijj[3] = 0;
    int v = ijbo(1,2);
    std::printf("ijbo(1,2)=%d %s\n", v, v==0?"PASS":"FAIL");
    return v==0?0:1;
}
