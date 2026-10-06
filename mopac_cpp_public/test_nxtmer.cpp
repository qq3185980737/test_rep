// test_nxtmer.cpp
#include <cstdio>
#include <vector>
#include "nxtmer.h"
#include "common_arrays_C.h"
namespace common_arrays_C {
    std::vector<int> nat;
    std::vector<int> nbonds;
    std::vector<std::vector<int>> ibonds;
}
using namespace common_arrays_C;
int main() {
    nat.assign(5,0); nbonds.assign(5,0);
    ibonds.assign(20, std::vector<int>(5,0));
    nat[1]=7; nbonds[1]=1; ibonds[1][1]=2;
    nat[2]=6; nbonds[2]=1; ibonds[1][2]=1;
    int b[4]={0,0,0,0};
    nxtmer(1,b);
    std::printf("nxtmer b=(%d,%d,%d,%d) PASS\n",b[0],b[1],b[2],b[3]);
    return 0;
}
