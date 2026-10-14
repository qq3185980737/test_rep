// test_disp_DnX.cpp
#include <cstdio>
#include <vector>
#include "disp_DnX.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C; using namespace common_arrays_C;
    numat=2; nat.assign(3,0); nat[1]=17; nat[2]=7;
    double e = disp_DnX(false);
    bool ok = (e > 0.0);
    std::printf("disp_DnX=%g %s\n", e, ok?"PASS":"FAIL");
    print_post_scf_corrections();
    return ok?0:1;
}
