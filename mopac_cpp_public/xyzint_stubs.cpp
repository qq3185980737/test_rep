// xyzint_stubs.cpp — extern "C" forwarding stubs for the Fortran symbols
// bangle_ / dihed_ that xyzint.cpp references.  They forward to the C++
// translations (bangle.cpp / dihed.cpp), converting the column-major
// double* layout used by the Fortran callers into 1-based vector views.
#include "bangle.h"
#include "dihed.h"

#include <vector>

namespace {

std::vector<std::vector<double>> make_view(const double* xyz, int n) {
  std::vector<std::vector<double>> v(4, std::vector<double>(n + 1, 0.0));
  for (int d = 1; d <= 3; ++d)
    for (int i = 1; i <= n; ++i) v[d][i] = xyz[(i - 1) * 3 + (d - 1)];
  return v;
}

}  // namespace

extern "C" void bangle_(double* xyz, int* a, int* b, int* c, double* sum) {
  int n = *a;
  if (*b > n) n = *b;
  if (*c > n) n = *c;
  auto v = make_view(xyz, n);
  double ang = 0.0;
  bangle(v, *a, *b, *c, ang);
  *sum = ang;
}

extern "C" void dihed_(double* xyz, int* a, int* b, int* c, int* d,
                       double* dih) {
  int n = *a;
  int aa[4] = {*a, *b, *c, *d};
  for (int q = 0; q < 4; ++q)
    if (aa[q] > n) n = aa[q];
  auto v = make_view(xyz, n);
  double ang = 0.0;
  dihed(v, *a, *b, *c, *d, ang);
  *dih = ang;
}
