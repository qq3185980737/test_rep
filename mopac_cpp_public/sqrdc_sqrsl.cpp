// sqrdc_sqrsl.cpp — C++ port of LINPACK SQRDC / SQRSL (single precision),
// LAPACK-style Householder storage:
//   after sqrdc, R(i,j) (i<=j) lives at x[(j-1)*ldx + (i-1)];
//   the Householder vector v (v[1]=1 implicit) lives at x[(j-1)*ldx + (i-1)]
//   for i>j, and tau (qraux) is stored in qraux(j).
// Column-major 1-based indexing mirrors how new_esp.F90 feeds its flat matrix.
#include "sqrdc_sqrsl.h"
#include <algorithm>
#include <cmath>

namespace {

float col_nrm2(const std::vector<float>& x, int ldx, int n, int jp, int j) {
  double s = 0.0;
  for (int i = j; i <= n; ++i) {
    float v = x[(jp - 1) * ldx + i];
    s += (double)v * v;
  }
  return (float)std::sqrt(s);
}

void col_swap(std::vector<float>& x, int ldx, int n, int j, int jp) {
  for (int i = 1; i <= n; ++i)
    std::swap(x[(j - 1) * ldx + i], x[(jp - 1) * ldx + i]);
}

}  // namespace

void sqrdc(std::vector<float>& x, int ldx, int n, int p,
           std::vector<float>& qraux, std::vector<int>& jpvt,
           std::vector<float>& work, int job) {
  int j, jp, lup;
  bool pivoting = (job != 0);
  if (pivoting) {
    // LINPACK: jpvt(k)=0 marks a free column; the routine first records the
    // original column indices so pivot interchanges leave a true permutation.
    for (j = 1; j <= p; ++j) jpvt[j] = j;
  }
  lup = std::min(n, p);
  for (j = 1; j <= lup; ++j) {
    if (pivoting) {
      for (jp = j; jp <= p; ++jp) work[jp] = col_nrm2(x, ldx, n, jp, j);
      int maxj = j;
      for (jp = j + 1; jp <= p; ++jp)
        if (work[jp] > work[maxj]) maxj = jp;
      if (maxj != j) {
        col_swap(x, ldx, n, j, maxj);
        std::swap(jpvt[j], jpvt[maxj]);
        work[maxj] = work[j];
      }
    }
    float nrm2 = col_nrm2(x, ldx, n, j, j);
    float x1 = x[(j - 1) * ldx + j];
    if (nrm2 != 0.0f && x1 != 0.0f) {
      float alpha = (x1 < 0.0f) ? nrm2 : -nrm2;   // -sign(x1)*nrm2
      float tau = (alpha - x1) / alpha;
      float v1 = x1 - alpha;
      // Build v components (i>j) and apply H = I - tau v v^T to cols j+1..p.
      for (jp = j + 1; jp <= p; ++jp) {
        float t = x[(jp - 1) * ldx + j];   // v[1]=1 * x(j,jp)
        for (int i = j + 1; i <= n; ++i)
          t += x[(j - 1) * ldx + i] / v1 * x[(jp - 1) * ldx + i];
        t *= tau;
        x[(jp - 1) * ldx + j] -= t;
        for (int i = j + 1; i <= n; ++i)
          x[(jp - 1) * ldx + i] -= t * x[(j - 1) * ldx + i] / v1;
      }
      // Store R(j,j) and the normalized Householder vector.
      x[(j - 1) * ldx + j] = alpha;
      for (int i = j + 1; i <= n; ++i)
        x[(j - 1) * ldx + i] /= v1;
      qraux[j] = tau;
    } else {
      x[(j - 1) * ldx + j] = x1;  // keep as R(j,j)
      qraux[j] = 0.0f;
    }
  }
  if (lup < p) {
    for (j = lup + 1; j <= p; ++j) qraux[j] = 0.0f;
  }
}

void sqrsl(std::vector<float>& x, int ldx, int n, int k,
           std::vector<float>& qraux, std::vector<float>& y,
           std::vector<float>& qy, std::vector<float>& qty,
           std::vector<float>& b, std::vector<float>& rsd,
           std::vector<float>& xb, int job, int& info) {
  info = 0;
  // LINPACK job coding: cqy = job/10000, cqty = mod(job,10000),
  // cb = mod(job,1000)/100, cr = mod(job,100)/10, cxb = mod(job,10).
  bool cqy = (job / 10000 != 0);
  bool cqty = (job % 10000 != 0);
  bool cb = ((job % 1000) / 100 != 0);
  bool crsd = ((job % 100) / 10 != 0);
  bool cxb = (job % 10 != 0);
  if (!cqy && !cqty && !cb && !crsd && !cxb) return;
  if (k < 1 || k > n) return;
  // Apply the Householder transforms to y -> qy, qty.
  for (int j = 1; j <= k; ++j) {
    float tau = qraux[j];
    if (tau != 0.0f) {
      float t = y[j];
      for (int i = j + 1; i <= n; ++i)
        t += x[(j - 1) * ldx + i] * y[i];
      t *= tau;
      y[j] -= t;
      for (int i = j + 1; i <= n; ++i)
        y[i] -= t * x[(j - 1) * ldx + i];
    }
  }
  if (cqy || cqty) {
    for (int i = 1; i <= n; ++i) {
      qy[i] = y[i];
      qty[i] = y[i];
    }
  }
  // Solve R xb = qty(1..k); R upper triangular at x(i,j), i<=j.
  if (cxb || crsd || cb) {
    for (int j = k; j >= 1; --j) {
      xb[j] = qty[j] / x[(j - 1) * ldx + j];
      for (int i = 1; i <= j - 1; ++i)
        qty[i] -= x[(j - 1) * ldx + i] * xb[j];
    }
  }
  if (cb) {
    for (int j = 1; j <= k; ++j) b[j] = xb[j];
  }

  if (crsd) {
    for (int i = 1; i <= n; ++i) {
      float s = 0.0f;
      for (int j = 1; j <= k; ++j) s += x[(j - 1) * ldx + i] * xb[j];
      rsd[i] = qy[i] - s;
    }
  }
}
