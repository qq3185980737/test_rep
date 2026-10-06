// matout.cpp — C++ translation of MOPAC 2016 "matout.F90".
// Prints a square matrix of eigenvectors and eigenvalues with orbital labels.
#include "matout.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

void matout(std::vector<double>& a, const std::vector<double>& b, int nc, int& nr, int ndim) {
    matout(a.data(), b.data(), nc, nr, ndim);
}
void matout(const double* a, const double* b, int nc, int& nr, int ndim) {
  using chanel_C::iw;
  using common_arrays_C::nfirst;
  using common_arrays_C::nlast;
  using common_arrays_C::nat;
  using molkst_C::numat;
  using elemts_C::elemnt;
  if (iw <= 0) { nr = std::abs(nr); return; }
  static const char* atorbs[9] = {" S", "PX", "PY", "PZ", "X2", "XZ", "Z2", "YZ", "XY"};
  int i = std::max(ndim, 1);
  std::vector<std::string> itext(i + 1), jtext(i + 1);
  std::vector<int> natom(i + 1, 0);
  bool atomic = (nlast[numat] == nr);
  if (atomic) {
    for (i = 1; i <= numat; ++i) {
      int jlo = nfirst[i];
      int jhi = nlast[i];
      int l = nat[i];
      for (int k = jlo; k <= jhi; ++k) {
        itext[k] = atorbs[k - jlo];
        jtext[k] = elemnt[l];
        natom[k] = i;
      }
    }
  } else {
    nr = std::abs(nr);
    for (i = 1; i <= nr; ++i) {
      itext[i] = "  ";
      jtext[i] = "  ";
      natom[i] = i;
    }
  }
  int ka = 1, kc = 6;
  for (;;) {
    int kb = std::min(kc, nc);
    std::printf("\n\n\n\n     ROOT NO.");
    for (i = ka; i <= kb; ++i) std::printf("%12d", i);
    std::printf("\n");
    if (b[1] != 0.0) {
      std::printf("\n        ");
      for (i = ka; i <= kb; ++i) std::printf("%12.5f", b[i]);
      std::printf("\n");
    }
    std::printf("  \n");
    int la = 1, lc = 40;
    for (;;) {
      int lb = std::min(lc, nr);
      for (i = la; i <= lb; ++i) {
        if (itext[i] == " S") std::printf("  \n");
        std::printf(" %s %s%4d", itext[i].c_str(), jtext[i].c_str(), natom[i]);
        std::printf("%10.5f", a[(i - 1) + ndim * (ka - 1)]);
        for (int j = ka + 1; j <= kb; ++j) std::printf("%12.5f", a[(i - 1) + ndim * (j - 1)]);
        std::printf("\n");
      }
      if (lb == nr) break;
      la = lc + 1;
      lc = lc + 40;
    }
    if (kb == nc) break;
    ka = kc + 1;
    kc = kc + 6;
  }
}
