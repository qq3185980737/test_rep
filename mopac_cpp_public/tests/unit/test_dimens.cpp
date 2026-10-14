// test_dimens.cpp
#include <cstdio>
#include <vector>
#include "dimens.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numat=2;
    std::vector<std::vector<double>> coord(3, std::vector<double>(3,0.0));
    coord[0][1]=0; coord[0][2]=3;
    dimens(coord,0);
    std::printf("dimens(skeleton) PASS\n");
    return 0;
}
