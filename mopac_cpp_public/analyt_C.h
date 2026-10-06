// analyt_C.h — C++ mapping of Fortran module "analyt_C".
#pragma once
namespace analyt_C {
constexpr int NZTYPE_MAX = 107;
extern int nztype[NZTYPE_MAX + 1];   // 1-based
extern double ds[17];                // dimension(16)
extern double dg[23];                // dimension(22)
extern double dr[101];               // dimension(100)
extern double g[23];                 // dimension(22)
extern double tx[4], ty[4], tz[4];         // dimension(3)
extern double tdx[4], tdy[4], tdz[4];      // dimension(3)
}
