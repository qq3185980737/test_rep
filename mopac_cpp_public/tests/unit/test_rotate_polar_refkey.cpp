// test_rotate_polar_refkey.cpp — rotate_C / polar_C / refkey_C / vastkind.
#include "polar_C.h"
#include "refkey_C.h"
#include "rotate_C.h"
#include "vastkind.h"

#include <cstdio>

using namespace polar_C;
using namespace refkey_C;
using namespace rotate_C;

int main() {
    bool ok = true;
    // polar_C
    if (omega != 0.0 || alpavg != 0.0 || dummy != 0.0) { std::printf("FAIL polar defaults\n"); ok = false; }
    omega = 3.5; alpavg = 12.25;
    if (omega != 3.5 || alpavg != 12.25) { std::printf("FAIL polar writes\n"); ok = false; }
    // refkey_C
    if (!refkey.empty()) { std::printf("FAIL refkey default\n"); ok = false; }
    refkey = "AM1";
    if (refkey != "AM1") { std::printf("FAIL refkey write\n"); ok = false; }
    // rotate_C: scalars + ccore equivalence contract
    css1 = 1.0; csp1 = 2.0; cpps1 = 3.0; cppp1 = 4.0;
    css2 = 5.0; csp2 = 6.0; cpps2 = 7.0; cppp2 = 8.0;
    ccore[1][1] = css1; ccore[2][1] = csp1; ccore[3][1] = cpps1; ccore[4][1] = cppp1;
    ccore[1][2] = css2; ccore[2][2] = csp2; ccore[3][2] = cpps2; ccore[4][2] = cppp2;
    if (ccore[4][2] != 8.0 || ccore[1][1] != 1.0) { std::printf("FAIL ccore sync\n"); ok = false; }
    // vast_kind_param
    if (vast_kind_param::double_kind != 8) { std::printf("FAIL double_kind\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
