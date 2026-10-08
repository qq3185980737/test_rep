// probe_hatom.cpp — 打印单 H 的 uss/uspd/h/eigs 中间值
#include <cstdio>
#include <vector>
#include "molkst_C.h"
#include "parameters_C.h"

// externs in molkst_C.cpp / parameters_C.cpp
namespace molkst_C {
    extern int numat, norbs, mpack, numcal;
    extern std::vector<double> h, pa, pb, p, f;
}
namespace parameters_C {
    extern double uss[];
}
extern double uspd[];   // moldat.cpp 定义? 查

int main() {
    // 仅打印 extern 符号不可行——直接引全局
    using namespace molkst_C;
    std::printf("norbs=%d numat=%d\n", norbs, numat);
    std::printf("uss[1]=%.8f\n", parameters_C::uss[1]);
    std::printf("h[1]=%.8f\n", h[1]);
    std::printf("pa[1]=%.8f pb[1]=%.8f p[1]=%.8f\n", pa[1], pb[1], p[1]);
    return 0;
}
