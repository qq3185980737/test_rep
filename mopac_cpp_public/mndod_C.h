// mndod_C.h — C++ translation of MOPAC 2016 "mndod_C.F90" (MNDO/d tables).
// Ported as std::vector containers (1-based indexing, matching mndod.cpp).
#pragma once

#include <vector>

namespace mndod_C {
// Integer data tables (Fortran data statements)
extern std::vector<int> iii;    // (107)
extern std::vector<int> iiid;   // (107)
extern std::vector<int> intij;  // (243)
extern std::vector<int> intkl;  // (243)
extern std::vector<int> intrep; // (243)
extern std::vector<int> isym;   // (491), no data statement; runtime filled

// Index arrays (filled at runtime)
extern std::vector<std::vector<int>> indexd; // (9,9)
extern std::vector<std::vector<int>> indx;   // (9,9)
extern std::vector<std::vector<int>> indpp;  // (3,3)
extern std::vector<std::vector<int>> inddp;  // (5,3)
extern std::vector<std::vector<int>> inddd;  // (5,5)
extern std::vector<int> iaf;                 // (500)
extern std::vector<int> ial;                 // (500)
extern std::vector<std::vector<int>> ind2;   // (45,45)
extern int nalp;

// Double arrays (runtime filled)
extern std::vector<std::vector<std::vector<double>>> ch; // (45,0:2,-2:2)
extern std::vector<double> alpb;   // (500)
extern std::vector<double> xfac;   // (500)
extern std::vector<std::vector<double>> aij;   // (6,107)
extern std::vector<std::vector<double>> repd;  // (52,107)
// thread_local: rewritten on every rotatd/rotmat/spcore pair call; the dcart
// OpenMP atom-pair loop must not share them across threads.
// thread_local POD arrays: rewritten on every rotatd/rotmat/spcore pair call.
extern thread_local double cored[11][3]; // (10,2)
extern thread_local double sp[4][4];     // (3,3)
extern thread_local double sd[6][6];     // (5,5)
extern thread_local double pp[7][4][4];  // (6,3,3)
extern thread_local double dp[16][6][4]; // (15,5,3)
extern thread_local double d_d[16][6][6];// (15,5,5)
extern std::vector<double> fx;               // (30) factorials fx(i)=i!
extern std::vector<std::vector<double>> b;   // (30,30) Pascal triangle
}  // namespace mndod_C
