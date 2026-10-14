// test_batchM04f.cpp  M04 capcor (capped-bond electronic correction) batch.
// Three hand-derived scenes: a lone cap atom (full-row path), a normal atom
// row picking up cap contributions, and the reverse order. The factor -2
// doubles the lower-half sum (capcor appears in 1/2*P(H+F)).
#include <cmath>
#include <cstdio>
#include <vector>

#include "capcor.h"
#include "molkst_C.h"

using namespace molkst_C;

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg)                                     \
  do {                                                    \
    if (c) {                                              \
      ++n_pass;                                           \
      std::printf("  [PASS] %s\n", msg);                  \
    } else {                                              \
      ++n_fail;                                           \
      std::printf("  [FAIL] %s\n", msg);                  \
    }                                                     \
  } while (0)
#define CHK_D(a, b, tol, msg)                                            \
  do {                                                                   \
    double va = (a), vb = (b);                                           \
    if (std::fabs(va - vb) <= (tol)) {                                   \
      ++n_pass;                                                          \
      std::printf("  [PASS] %s (%.6g)\n", msg, va);                      \
    } else {                                                             \
      ++n_fail;                                                          \
      std::printf("  [FAIL] %s: got %.6g want %.6g\n", msg, va, vb);     \
    }                                                                    \
  } while (0)

int main() {
  std::printf("M04f batch - capcor\n");

  // Test 1: single cap atom (nat=102), nlast=2 -> j=3, one step j=2.
  {
    numat = 1;
    std::vector<int> nat = {0, 102};
    std::vector<int> nfirst = {0, 1};
    std::vector<int> nlast = {0, 2};
    std::vector<double> p(4, 0.0), h(4, 0.0);
    p[2] = 1.5; h[2] = 2.0;
    double r = capcor(nat, nfirst, nlast, p, h);
    CHK_D(r, -2.0 * 1.5 * 2.0, 1e-12, "capcor lone cap atom");
  }

  // Test 2: atom 2 is cap with nlast=3 -> j=6, two steps j=5,j=4.
  {
    numat = 2;
    std::vector<int> nat = {0, 1, 102};
    std::vector<int> nfirst = {0, 1, 2};
    std::vector<int> nlast = {0, 1, 3};
    std::vector<double> p(7, 0.0), h(7, 0.0);
    p[4] = 1.0; h[4] = 3.0;
    p[5] = 2.0; h[5] = 4.0;
    double r = capcor(nat, nfirst, nlast, p, h);
    CHK_D(r, -2.0 * (1.0 * 3.0 + 2.0 * 4.0), 1e-12, "capcor cap row two terms");
  }

  // Test 3: cap atom first, normal atom second -> off-diagonal kk=2.
  {
    numat = 2;
    std::vector<int> nat = {0, 102, 1};
    std::vector<int> nfirst = {0, 1, 2};
    std::vector<int> nlast = {0, 1, 2};
    std::vector<double> p(4, 0.0), h(4, 0.0);
    p[2] = 0.5; h[2] = 6.0;
    double r = capcor(nat, nfirst, nlast, p, h);
    // i=1 (cap, iu=1): ii=0 -> nothing. i=2: j=1 cap -> kk=1+1=2.
    CHK_D(r, -2.0 * 0.5 * 6.0, 1e-12, "capcor off-diagonal kk=2");
  }

  // Test 4: no cap atoms -> zero.
  {
    numat = 2;
    std::vector<int> nat = {0, 1, 6};
    std::vector<int> nfirst = {0, 1, 2};
    std::vector<int> nlast = {0, 1, 2};
    std::vector<double> p(4, 1.0), h(4, 1.0);
    double r = capcor(nat, nfirst, nlast, p, h);
    CHK_D(r, 0.0, 1e-15, "capcor no cap atoms -> 0");
  }

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
