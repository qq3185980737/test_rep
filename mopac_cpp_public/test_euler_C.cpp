// test_euler_C.cpp — euler_C module: l1l/l2l/l3l default-zero and writable;
// also cross-check the split globals (molkst_C::l1u..l3u, common_arrays_C::tvec).
#include "euler_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

#include <cstdio>

using namespace euler_C;
using namespace molkst_C;
using namespace common_arrays_C;

int main() {
    bool ok = true;
    if (l1l != 0 || l2l != 0 || l3l != 0) { std::printf("FAIL default zero\n"); ok = false; }
    l1l = 2; l2l = -3; l3l = 5;
    if (l1l != 2 || l2l != -3 || l3l != 5) { std::printf("FAIL write\n"); ok = false; }
    // split globals are reachable through their home namespaces
    l1u = 1; l2u = 1; l3u = 1; l123 = 1;
    if (l1u != 1 || l123 != 1) { std::printf("FAIL molkst_C cells\n"); ok = false; }
    tvec.assign(4, std::vector<double>(4, 0.0));
    tvec[1][1] = 10.0; tvec[2][2] = 10.0; tvec[3][3] = 10.0;
    if (tvec[3][3] != 10.0) { std::printf("FAIL tvec\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
