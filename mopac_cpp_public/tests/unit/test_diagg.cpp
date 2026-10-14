// test_diagg.cpp
#include <cstdio>
#include <vector>
#include "diagg.h"
int main() {
    std::vector<double> fao(3), partp(3);
    diagg(fao, 0, 0, 0, partp, 0);
    std::printf("diagg(skeleton) PASS\n");
    return 0;
}
