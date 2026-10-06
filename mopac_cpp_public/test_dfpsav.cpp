// test_dfpsav.cpp
#include <cstdio>
#include <vector>
#include "dfpsav.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numcal=1;
    std::vector<double> xp(2), gd(2), xl(2), xdfp(10);
    std::vector<int> mdfp(10,0);
    dfpsav(0, xp, gd, xl, 0, mdfp, xdfp);
    std::printf("dfpsav(skeleton) PASS\n");
    return 0;
}
