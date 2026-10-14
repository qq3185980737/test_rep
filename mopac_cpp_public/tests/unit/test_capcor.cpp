// test_capcor.cpp
#include <cstdio>
#include <vector>

#include "capcor.h"
#include "molkst_C.h"

int main() {
    molkst_C::numat = 2;
    std::vector<int> nat = {0, 6, 102};
    std::vector<int> nfirst = {0, 1, 2};
    std::vector<int> nlast = {0, 1, 2};
    std::vector<double> p = {0, 1.0, 3.0, 5.0};
    std::vector<double> h = {0, 2.0, 4.0, 6.0};

    double r = capcor(nat, nfirst, nlast, p, h);
    // i=2 cap: j=3, k=1 -> j=2 -> p(2)*h(2)=3*4=12; result=-24.
    bool ok = (r == -24.0);
    std::printf("capcor=%.1f (expect -24.0) %s\n", r, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
