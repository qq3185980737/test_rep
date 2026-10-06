// test_copyc6.cpp — verify C6 data table integrity + limit/copy logic.
#include <cstdio>
#include <vector>
#include "copyc6.h"
int main() {
  const double* p = copyc6_data::pars;
  // spot checks from source data lines
  struct { int k; double v; } chk[] = {
    {0, 1.0}, {2, 3.0267}, {3, 0.9118}, {4, 0.9118},           // pars(1:5)
    {161920, 482.0}, {161922, 455.2854}, {161924, 3.9098},    // pars(161921:161925)
  };
  int ok = 1;
  for (auto& c : chk) {
    double d = p[c.k];
    if (d != c.v) { std::printf("FAIL pars[%d]=%g want %g\n", c.k, d, c.v); ok = 0; }
  }
  // run copyc6 with max_elem=500, maxc=6 to cover ids up to 482
  int maxc=6, max_elem=500;
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> c6ab(
    max_elem+1, std::vector<std::vector<std::vector<std::vector<double>>>>(max_elem+1,
    std::vector<std::vector<std::vector<double>>>(maxc+1,
    std::vector<std::vector<double>>(maxc+1, std::vector<double>(4, 0.0)))));
  std::vector<int> maxci(max_elem+1, 0);
  copyc6(maxc, max_elem, c6ab, maxci);
  // c6ab(1,1,1,1,*) = 3.0267, 0.9118, 0.9118
  if (c6ab[1][1][1][1][0] != 3.0267 || c6ab[1][1][1][1][1] != 0.9118) { std::printf("FAIL c6ab(1,1)\n"); ok = 0; }
  // c6ab(482,482,5,5,*) = 455.2854, 3.9098, 3.9098 ; maxci(482)=5
  if (c6ab[482][482][5][5][0] != 455.2854 || c6ab[482][482][5][5][1] != 3.9098) { std::printf("FAIL c6ab(482,482)\n"); ok = 0; }
  if (maxci[482] != 5 || maxci[1] != 1) { std::printf("FAIL maxci\n"); ok = 0; }
  // symmetry: c6ab(jat,iat) swaps comp 2/3
  if (c6ab[2][1][1][1][0] != c6ab[1][2][1][1][0]) { std::printf("FAIL sym c6\n"); ok = 0; }
  // untouched cells stay -1
  if (c6ab[10][10][2][2][0] != -1.0) { std::printf("FAIL -1 init\n"); ok = 0; }
  std::printf(ok ? "ALL PASS\n" : "FAILED\n");
  return ok ? 0 : 1;
}
