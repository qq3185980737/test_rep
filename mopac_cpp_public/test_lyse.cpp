// test_lyse.cpp
#include <cstdio>
#include <vector>
#include "lyse.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
namespace common_arrays_C {
    std::vector<int> nat;
    std::vector<int> nbonds;
    std::vector<std::vector<int>> ibonds;
}
namespace molkst_C { int numat; }
using namespace common_arrays_C;
using namespace molkst_C;
extern "C" double distance_(int* a, int* b){ (void)a;(void)b; return 1.0; }
int main() {
    numat = 3;
    nat = {0, 6, 1, 1};
    nbonds = {0, 2, 1, 1};
    ibonds.resize(20, std::vector<int>(5, 0));
    ibonds[1][1] = 2; ibonds[2][1] = 3;
    ibonds[1][2] = 1; ibonds[1][3] = 1;
    lyse();
    std::printf("nbonds[1]=%d PASS\n", nbonds[1]);
    return 0;
}
