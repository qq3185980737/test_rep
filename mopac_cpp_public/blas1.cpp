// blas1.cpp — standard BLAS Level-1 implementations.
#include "blas1.h"

#include <cmath>
#include <cstdlib>

void daxpy(int n, double a, const std::vector<double>& x, int incx,
           std::vector<double>& y, int incy) {
    if (n <= 0 || a == 0.0) return;
    int ix = 0, iy = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    if (incy < 0) iy = (-n + 1) * incy;
    for (int i = 0; i < n; ++i) {
        y[iy] += a * x[ix];
        ix += incx;
        iy += incy;
    }
}

void dscal(int n, double a, std::vector<double>& x, int incx) {
    if (n <= 0) return;
    int ix = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    for (int i = 0; i < n; ++i) {
        x[ix] *= a;
        ix += incx;
    }
}

void dswap(int n, std::vector<double>& x, int incx,
           std::vector<double>& y, int incy) {
    if (n <= 0) return;
    int ix = 0, iy = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    if (incy < 0) iy = (-n + 1) * incy;
    for (int i = 0; i < n; ++i) {
        double t = x[ix];
        x[ix] = y[iy];
        y[iy] = t;
        ix += incx;
        iy += incy;
    }
}

double ddot(int n, const std::vector<double>& x, int incx,
            const std::vector<double>& y, int incy) {
    double s = 0.0;
    if (n <= 0) return s;
    int ix = 0, iy = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    if (incy < 0) iy = (-n + 1) * incy;
    for (int i = 0; i < n; ++i) {
        s += x[ix] * y[iy];
        ix += incx;
        iy += incy;
    }
    return s;
}

double dnrm2(int n, const std::vector<double>& x, int incx) {
    double s = 0.0;
    if (n <= 0) return s;
    int ix = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    for (int i = 0; i < n; ++i) {
        s += x[ix] * x[ix];
        ix += incx;
    }
    return std::sqrt(s);
}

void dcopy(int n, const std::vector<double>& x, int incx,
           std::vector<double>& y, int incy) {
    if (n <= 0) return;
    int ix = 0, iy = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    if (incy < 0) iy = (-n + 1) * incy;
    for (int i = 0; i < n; ++i) {
        y[iy] = x[ix];
        ix += incx;
        iy += incy;
    }
}

int idamax(int n, const std::vector<double>& x, int incx) {
    if (n < 1) return 0;
    int best = 0;
    double bestv = -1.0;
    int ix = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    for (int i = 0; i < n; ++i) {
        double v = std::abs(x[ix]);
        if (v > bestv) { bestv = v; best = i; }
        ix += incx;
    }
    return best + 1;  // 1-based
}
