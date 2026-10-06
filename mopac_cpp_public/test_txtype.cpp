// test_txtype.cpp
#include <cstdio>
#include <string>
#include <vector>
#include "txtype.h"
#include "common_arrays_C.h"
namespace common_arrays_C {
    std::vector<int> nat;
    std::vector<int> nbonds;
    std::vector<std::vector<int>> ibonds;
    std::vector<std::string> txtatm;
}
using namespace common_arrays_C;
int main() {
    nat = {0, 6, 6, 1};
    nbonds = {0, 2, 2, 0};
    txtatm.resize(4, "............XXX...UNK");
    ibonds.resize(20, std::vector<int>(4, 0));
    int jt[3]={1,2,3}, jj=3;
    txtype(jj, jt, 'A');
    std::printf("jj=%d atom1[14]=%c atom2[14]=%c PASS\n", jj, txtatm[1][14], txtatm[2][14]);
    return 0;
}
