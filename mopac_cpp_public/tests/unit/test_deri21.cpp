// test_deri21.cpp
#include <cstdio>
#include <vector>
#include "deri21.h"
int main() {
    std::vector<std::vector<double>> a(2, std::vector<double>(2, 0.0));
    std::vector<double> vnert(3,0), pnert(3,0);
    std::vector<std::vector<double>> b(2, std::vector<double>(2,0));
    int ncut = 0;
    deri21(a, 1, 1, 1, vnert, pnert, b, ncut);
    bool ok = (ncut == 1 && vnert[1] == 1.0);
    std::printf("ncut=%d vnert1=%g %s\n", ncut, vnert[1], ok?"PASS":"FAIL");
    return ok?0:1;
}
