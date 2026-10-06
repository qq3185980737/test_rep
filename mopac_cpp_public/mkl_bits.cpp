// mkl_bits.cpp — C++ translation of MOPAC 2016 "mkl_bits.F90".
// BLAS level-1/2/3 and LAPACK routines, column-major Fortran layout.
// Semantics follow the standard Netlib BLAS/LAPACK definitions that
// mkl_bits.F90 re-implements (ddot/dcopy/dscal/daxpy/dswap/dnrm2,
// dgemm/dgemv/dger/dsyrk/dtrmm/dtrmv/dtrsm, dgesv/dgetf2/dgetrf/dgetrs/
// dlaswp, dpotf2/dpotrf/dpotri/dlauu2/dlauum/dtrti2/dtrtri, xerbla).
#include "mkl_bits.h"
#include <cmath>
#include <cstdio>

// Case-insensitive character compare (mirrors LAPACK LSAME).
static bool lsame(char a, char b) {
    if (a >= 'a' && a <= 'z') a = a - 'a' + 'A';
    if (b >= 'a' && b <= 'z') b = b - 'a' + 'A';
    return a == b;
}
static bool upper(char u) { return lsame(u, 'U'); }

// ------------------------------------------------------------------ level 1
double ddot(int n, const double* dx, int incx, const double* dy, int incy) {
    double s = 0.0;
    if (n <= 0) return s;
    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) s += dx[i] * dy[i];
    } else {
        int ix = 0, iy = 0;
        if (incx < 0) ix = (-n + 1) * incx;
        if (incy < 0) iy = (-n + 1) * incy;
        for (int i = 0; i < n; ++i) { s += dx[ix] * dy[iy]; ix += incx; iy += incy; }
    }
    return s;
}

void dcopy(int n, const double* dx, int incx, double* dy, int incy) {
    if (n <= 0) return;
    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) dy[i] = dx[i];
    } else {
        int ix = 0, iy = 0;
        if (incx < 0) ix = (-n + 1) * incx;
        if (incy < 0) iy = (-n + 1) * incy;
        for (int i = 0; i < n; ++i) { dy[iy] = dx[ix]; ix += incx; iy += incy; }
    }
}

void dscal(int n, double da, double* dx, int incx) {
    if (n <= 0 || incx <= 0) return;
    for (int i = 0; i < n; ++i) dx[i * incx] *= da;
}

void daxpy(int n, double da, const double* dx, int incx, double* dy, int incy) {
    if (n <= 0 || da == 0.0) return;
    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) dy[i] += da * dx[i];
    } else {
        int ix = 0, iy = 0;
        if (incx < 0) ix = (-n + 1) * incx;
        if (incy < 0) iy = (-n + 1) * incy;
        for (int i = 0; i < n; ++i) { dy[iy] += da * dx[ix]; ix += incx; iy += incy; }
    }
}

void dswap(int n, double* dx, int incx, double* dy, int incy) {
    if (n <= 0) return;
    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) { double t = dx[i]; dx[i] = dy[i]; dy[i] = t; }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0) ix = (-n + 1) * incx;
        if (incy < 0) iy = (-n + 1) * incy;
        for (int i = 0; i < n; ++i) { double t = dx[ix]; dx[ix] = dy[iy]; dy[iy] = t; ix += incx; iy += incy; }
    }
}

double dnrm2(int n, const double* x, int incx) {
    double s = 0.0;
    if (n <= 0) return s;
    int ix = 0;
    if (incx < 0) ix = (-n + 1) * incx;
    for (int i = 0; i < n; ++i) { s += x[ix] * x[ix]; ix += incx; }
    return std::sqrt(s);
}

// ------------------------------------------------------------------ level 2
void dgemv(char trans, int m, int n, double alpha, const double* a, int lda,
           const double* x, int incx, double beta, double* y, int incy) {
    if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0)) return;
    bool tr = lsame(trans, 'T') || lsame(trans, 'C');
    int len = tr ? m : n;
    for (int j = 0; j < len; ++j) {
        if (beta == 0.0) y[j * incy] = 0.0;
        else y[j * incy] *= beta;
    }
    if (alpha == 0.0) return;
    if (!tr) {
        // y = alpha*A*x + beta*y ; A is m x n column-major
        for (int j = 0; j < n; ++j) {
            double xj = x[j * incx];
            if (xj == 0.0) continue;
            for (int i = 0; i < m; ++i) y[i * incy] += alpha * xj * a[i + j * lda];
        }
    } else {
        // y = alpha*A^T*x + beta*y
        for (int i = 0; i < m; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) sum += a[i + j * lda] * x[j * incx];
            y[i * incy] += alpha * sum;
        }
    }
}

void dger(int m, int n, double alpha, const double* x, int incx,
          const double* y, int incy, double* a, int lda) {
    if (m == 0 || n == 0 || alpha == 0.0) return;
    for (int j = 0; j < n; ++j) {
        double yj = alpha * y[j * incy];
        if (yj == 0.0) continue;
        for (int i = 0; i < m; ++i) a[i + j * lda] += x[i * incx] * yj;
    }
}

void dtrmv(char uplo, char trans, char diag, int n, const double* a, int lda,
           double* x, int incx) {
    if (n == 0) return;
    bool unit = lsame(diag, 'U');
    bool tr = lsame(trans, 'T') || lsame(trans, 'C');
    if (!tr) {
        if (upper(uplo)) {
            for (int j = 0; j < n; ++j) {
                double temp = x[j * incx];
                if (!unit && temp != 0.0) temp *= a[j + j * lda];
                for (int i = j + 1; i < n; ++i) x[i * incx] += temp * a[i + j * lda];
            }
        } else {
            for (int j = n - 1; j >= 0; --j) {
                double temp = x[j * incx];
                if (!unit && temp != 0.0) temp *= a[j + j * lda];
                for (int i = 0; i < j; ++i) x[i * incx] += temp * a[i + j * lda];
            }
        }
    } else {
        if (upper(uplo)) {
            for (int j = n - 1; j >= 0; --j) {
                double temp = x[j * incx];
                if (!unit && temp != 0.0) temp *= a[j + j * lda];
                for (int i = 0; i < j; ++i) x[i * incx] += temp * a[j + i * lda];
            }
        } else {
            for (int j = 0; j < n; ++j) {
                double temp = x[j * incx];
                if (!unit && temp != 0.0) temp *= a[j + j * lda];
                for (int i = j + 1; i < n; ++i) x[i * incx] += temp * a[j + i * lda];
            }
        }
    }
}

// ------------------------------------------------------------------ level 3
void dgemm(char transa, char transb, int m, int n, int k, double alpha,
           const double* a, int lda, const double* b, int ldb, double beta,
           double* c, int ldc) {
    bool ta = lsame(transa, 'T') || lsame(transa, 'C');
    bool tb = lsame(transb, 'T') || lsame(transb, 'C');
    if (m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0)) return;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i) {
            if (beta == 0.0) c[i + j * ldc] = 0.0;
            else c[i + j * ldc] *= beta;
        }
    if (alpha == 0.0) return;
    if (!ta && !tb) {
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l) {
                double blj = b[l + j * ldb];
                if (blj == 0.0) continue;
                for (int i = 0; i < m; ++i) c[i + j * ldc] += alpha * a[i + l * lda] * blj;
            }
    } else if (!ta && tb) {
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                double sum = 0.0;
                for (int l = 0; l < k; ++l) sum += a[i + l * lda] * b[j + l * ldb];
                c[i + j * ldc] += alpha * sum;
            }
    } else if (ta && !tb) {
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l) {
                double blj = b[l + j * ldb];
                if (blj == 0.0) continue;
                for (int i = 0; i < m; ++i) c[i + j * ldc] += alpha * a[l + i * lda] * blj;
            }
    } else {
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) {
                double sum = 0.0;
                for (int l = 0; l < k; ++l) sum += a[l + i * lda] * b[j + l * ldb];
                c[i + j * ldc] += alpha * sum;
            }
    }
}

void dsyrk(char uplo, char trans, int n, int k, double alpha, const double* a,
           int lda, double beta, double* c, int ldc) {
    if (n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0)) return;
    bool upper_ = upper(uplo);
    bool tr = lsame(trans, 'T') || lsame(trans, 'C');
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            if (!upper_ && i < j) continue;
            if (upper_ && j < i) continue;
            if (beta == 0.0) c[i + j * ldc] = 0.0;
            else c[i + j * ldc] *= beta;
        }
    if (alpha == 0.0) return;
    if (!tr) {
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l) {
                double al = alpha * a[j + l * lda];
                if (al == 0.0) continue;
                int i_lo = upper_ ? 0 : j, i_hi = upper_ ? j : n - 1;
                for (int i = i_lo; i <= i_hi; ++i) c[i + j * ldc] += a[i + l * lda] * al;
            }
    } else {
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < k; ++l) {
                double al = alpha * a[l + j * lda];
                if (al == 0.0) continue;
                int i_lo = upper_ ? 0 : j, i_hi = upper_ ? j : n - 1;
                for (int i = i_lo; i <= i_hi; ++i) c[i + j * ldc] += a[l + i * lda] * al;
            }
    }
}

void dtrmm(char side, char uplo, char transa, char diag, int m, int n,
           double alpha, const double* a, int lda, double* b, int ldb) {
    if (m == 0 || n == 0) return;
    bool left = lsame(side, 'L');
    bool upper_ = upper(uplo);
    bool tr = lsame(transa, 'T') || lsame(transa, 'C');
    bool unit = lsame(diag, 'U');
    if (left) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) b[i + j * ldb] *= alpha;
            if (upper_) {
                for (int kk = 0; kk < m; ++kk) {
                    double bkj = b[kk + j * ldb];
                    if (!unit && bkj != 0.0) bkj *= a[kk + kk * lda];
                    if (bkj != 0.0)
                        for (int i = 0; i < kk; ++i) b[i + j * ldb] += bkj * a[i + kk * lda];
                    b[kk + j * ldb] = bkj;
                }
            } else {
                for (int kk = m - 1; kk >= 0; --kk) {
                    double bkj = b[kk + j * ldb];
                    if (!unit && bkj != 0.0) bkj *= a[kk + kk * lda];
                    if (bkj != 0.0)
                        for (int i = kk + 1; i < m; ++i) b[i + j * ldb] += bkj * a[i + kk * lda];
                    b[kk + j * ldb] = bkj;
                }
            }
        }
    } else {
        for (int j = 0; j < n; ++j) {
            if (upper_) {
                for (int kk = 0; kk < n; ++kk) {
                    double bjk = b[j + kk * ldb];
                    if (!unit && bjk != 0.0) bjk *= a[kk + kk * lda];
                    if (bjk != 0.0)
                        for (int i = kk + 1; i < n; ++i) b[j + i * ldb] += bjk * a[i + kk * lda];
                    b[j + kk * ldb] = bjk;
                }
            } else {
                for (int kk = n - 1; kk >= 0; --kk) {
                    double bjk = b[j + kk * ldb];
                    if (!unit && bjk != 0.0) bjk *= a[kk + kk * lda];
                    if (bjk != 0.0)
                        for (int i = 0; i < kk; ++i) b[j + i * ldb] += bjk * a[i + kk * lda];
                    b[j + kk * ldb] = bjk;
                }
            }
        }
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i) b[i + j * ldb] *= alpha;
    }
}

void dtrsm(char side, char uplo, char transa, char diag, int m, int n,
           double alpha, const double* a, int lda, double* b, int ldb) {
    if (m == 0 || n == 0) return;
    bool left = lsame(side, 'L');
    bool upper_ = upper(uplo);
    bool tr = lsame(transa, 'T') || lsame(transa, 'C');
    bool unit = lsame(diag, 'U');
    if (left) {
        for (int j = 0; j < n; ++j) {
            if (upper_) {
                for (int kk = m - 1; kk >= 0; --kk) {
                    double bkj = b[kk + j * ldb] / (unit ? 1.0 : a[kk + kk * lda]);
                    b[kk + j * ldb] = bkj;
                    if (bkj != 0.0 && kk > 0)
                        for (int i = 0; i < kk; ++i) b[i + j * ldb] -= bkj * a[i + kk * lda];
                }
            } else {
                for (int kk = 0; kk < m; ++kk) {
                    double bkj = b[kk + j * ldb] / (unit ? 1.0 : a[kk + kk * lda]);
                    b[kk + j * ldb] = bkj;
                    if (bkj != 0.0 && kk + 1 < m)
                        for (int i = kk + 1; i < m; ++i) b[i + j * ldb] -= bkj * a[i + kk * lda];
                }
            }
            for (int i = 0; i < m; ++i) b[i + j * ldb] *= alpha;
        }
    } else {
        for (int j = 0; j < n; ++j) {
            if (upper_) {
                for (int kk = 0; kk < n; ++kk) {
                    double bjk = b[j + kk * ldb] / (unit ? 1.0 : a[kk + kk * lda]);
                    b[j + kk * ldb] = bjk;
                    if (bjk != 0.0 && kk + 1 < n)
                        for (int i = kk + 1; i < n; ++i) b[j + i * ldb] -= bjk * a[i + kk * lda];
                }
            } else {
                for (int kk = n - 1; kk >= 0; --kk) {
                    double bjk = b[j + kk * ldb] / (unit ? 1.0 : a[kk + kk * lda]);
                    b[j + kk * ldb] = bjk;
                    if (bjk != 0.0 && kk > 0)
                        for (int i = 0; i < kk; ++i) b[j + i * ldb] -= bjk * a[i + kk * lda];
                }
            }
            for (int i = 0; i < m; ++i) b[i + j * ldb] *= alpha;
        }
    }
}

// ------------------------------------------------------------------ LAPACK
void xerbla(const char* srname, int info) {
    std::fprintf(stderr, " ** On entry to %s parameter number %2d had an illegal value\n",
                 srname, info);
}

void dgetf2(int m, int n, double* a, int lda, int* ipiv, int& info) {
    info = 0;
    if (m < 0 || n < 0) { xerbla("DGETF2", 1); return; }
    for (int j = 0; j < n && j < m; ++j) {
        int jp = j;
        double amax = std::fabs(a[j + j * lda]);
        for (int i = j + 1; i < m; ++i) {
            double t = std::fabs(a[i + j * lda]);
            if (t > amax) { amax = t; jp = i; }
        }
        ipiv[j] = jp + 1;
        if (a[jp + j * lda] != 0.0) {
            if (jp != j) for (int c = 0; c < n; ++c) { double t = a[j + c * lda]; a[j + c * lda] = a[jp + c * lda]; a[jp + c * lda] = t; }
            double ajj = 1.0 / a[j + j * lda];
            for (int i = j + 1; i < m; ++i) a[i + j * lda] *= ajj;
            for (int c = j + 1; c < n; ++c) {
                double t = a[j + c * lda];
                for (int i = j + 1; i < m; ++i) a[i + c * lda] -= a[i + j * lda] * t;
            }
        } else if (info == 0) {
            info = j + 1;
        }
    }
}

void dgetrf(int m, int n, double* a, int lda, int* ipiv, int& info) {
    info = 0;
    if (m < 0 || n < 0) { xerbla("DGETRF", 1); return; }
    dgetf2(m, n, a, lda, ipiv, info);
}

void dlaswp(int n, double* a, int lda, int k1, int k2, const int* ipiv,
            int incx) {
    if (incx > 0) {
        for (int k = k1 - 1; k < k2; ++k) {
            int ip = ipiv[k] - 1;
            if (ip != k) for (int c = 0; c < n; ++c) { double t = a[k + c * lda]; a[k + c * lda] = a[ip + c * lda]; a[ip + c * lda] = t; }
        }
    } else {
        for (int k = k2 - 1; k >= k1 - 1; --k) {
            int ip = ipiv[k] - 1;
            if (ip != k) for (int c = 0; c < n; ++c) { double t = a[k + c * lda]; a[k + c * lda] = a[ip + c * lda]; a[ip + c * lda] = t; }
        }
    }
}

void dgetrs(char trans, int n, int nrhs, const double* a, int lda,
            const int* ipiv, double* b, int ldb, int& info) {
    info = 0;
    if (n < 0 || nrhs < 0) { xerbla("DGETRS", 1); return; }
    bool tr = lsame(trans, 'T') || lsame(trans, 'C');
    for (int j = 0; j < nrhs; ++j) {
        if (!tr) {
            // P*b
            for (int k = 0; k < n; ++k) {
                int ip = ipiv[k] - 1;
                if (ip != k) { double t = b[k + j * ldb]; b[k + j * ldb] = b[ip + j * ldb]; b[ip + j * ldb] = t; }
            }
            // L forward (unit lower)
            for (int k = 0; k < n - 1; ++k)
                for (int i = k + 1; i < n; ++i)
                    b[i + j * ldb] -= a[i + k * lda] * b[k + j * ldb];
            // U back
            for (int k = n - 1; k >= 0; --k) {
                for (int i = k + 1; i < n; ++i) b[k + j * ldb] -= a[k + i * lda] * b[i + j * ldb];
                b[k + j * ldb] /= a[k + k * lda];
            }
        } else {
            // U^T forward
            for (int k = 0; k < n; ++k) {
                for (int i = 0; i < k; ++i) b[k + j * ldb] -= a[i + k * lda] * b[i + j * ldb];
                b[k + j * ldb] /= a[k + k * lda];
            }
            // L^T back
            for (int k = n - 2; k >= 0; --k)
                for (int i = k + 1; i < n; ++i)
                    b[i + j * ldb] -= a[i + k * lda] * b[k + j * ldb];
            // P^T b
            for (int k = n - 1; k >= 0; --k) {
                int ip = ipiv[k] - 1;
                if (ip != k) { double t = b[k + j * ldb]; b[k + j * ldb] = b[ip + j * ldb]; b[ip + j * ldb] = t; }
            }
        }
    }
}

void dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb,
           int& info) {
    info = 0;
    if (n < 0 || nrhs < 0) { xerbla("DGESV", 1); return; }
    if (n == 0 || nrhs == 0) return;
    dgetrf(n, n, a, lda, ipiv, info);
    if (info == 0) dgetrs('N', n, nrhs, a, lda, ipiv, b, ldb, info);
}

void dpotf2(char uplo, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DPOTF2", 1); return; }
    bool upper_ = upper(uplo);
    if (upper_) {
        for (int j = 0; j < n; ++j) {
            double ajj = a[j + j * lda] - ddot(j, a + j * lda, 1, a + j * lda, 1);
            if (ajj <= 0.0) { info = j + 1; return; }
            ajj = std::sqrt(ajj);
            a[j + j * lda] = ajj;
            if (j + 1 < n) {
                for (int i = j + 1; i < n; ++i)
                    a[j + i * lda] -= ddot(j, a + i * lda, 1, a + j * lda, 1);
                double inv = 1.0 / ajj;
                for (int i = j + 1; i < n; ++i) a[j + i * lda] *= inv;
            }
        }
    } else {
        for (int j = 0; j < n; ++j) {
            double ajj = a[j + j * lda] - ddot(j, a + j * lda, lda, a + j * lda, lda);
            if (ajj <= 0.0) { info = j + 1; return; }
            ajj = std::sqrt(ajj);
            a[j + j * lda] = ajj;
            if (j + 1 < n) {
                for (int i = j + 1; i < n; ++i) {
                    double sum = a[i + j * lda];
                    for (int k = j; k < i; ++k) sum -= a[k + j * lda] * a[k + i * lda];
                    a[i + j * lda] = sum;
                }
                double inv = 1.0 / ajj;
                for (int i = j + 1; i < n; ++i) a[i + j * lda] *= inv;
            }
        }
    }
}

void dpotrf(char uplo, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DPOTRF", 1); return; }
    dpotf2(uplo, n, a, lda, info);
}

void dlauu2(char uplo, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DLAUU2", 1); return; }
    if (upper(uplo)) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < j; ++i) {
                double sum = a[i + j * lda] * a[j + j * lda];
                for (int k = 0; k < i; ++k) sum += a[k + i * lda] * a[k + j * lda];
                a[i + j * lda] = sum;
            }
            double sum = a[j + j * lda] * a[j + j * lda];
            for (int k = 0; k < j; ++k) sum += a[k + j * lda] * a[k + j * lda];
            a[j + j * lda] = sum;
        }
    } else {
        for (int j = n - 1; j >= 0; --j) {
            for (int i = j + 1; i < n; ++i) {
                double sum = a[j + i * lda] * a[j + j * lda];
                for (int k = j + 1; k < i; ++k) sum += a[k + i * lda] * a[k + j * lda];
                a[j + i * lda] = sum;
            }
            double sum = a[j + j * lda] * a[j + j * lda];
            for (int k = j + 1; k < n; ++k) sum += a[k + j * lda] * a[k + j * lda];
            a[j + j * lda] = sum;
        }
    }
}

void dlauum(char uplo, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DLAUUM", 1); return; }
    if (upper(uplo)) {
        for (int j = 0; j < n; ++j)
            for (int i = 0; i <= j; ++i) {
                double sum = 0.0;
                for (int k = j; k < n; ++k) sum += a[i + k * lda] * a[j + k * lda];
                a[i + j * lda] = sum;
            }
    } else {
        for (int j = n - 1; j >= 0; --j)
            for (int i = j; i < n; ++i) {
                double sum = 0.0;
                for (int k = 0; k <= j; ++k) sum += a[i + k * lda] * a[j + k * lda];
                a[i + j * lda] = sum;
            }
    }
}

void dtrti2(char uplo, char diag, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DTRTI2", 1); return; }
    bool unit = lsame(diag, 'U');
    if (upper(uplo)) {
        // X = U^{-1}: column j solved by back-substitution; a[i][i] already
        // holds X[i][i] for i < j (processed in reverse column order).
        for (int j = n - 1; j >= 0; --j) {
            if (!unit) {
                if (a[j + j * lda] == 0.0) { info = j + 1; return; }
                a[j + j * lda] = 1.0 / a[j + j * lda];
            }
            for (int i = j - 1; i >= 0; --i) {
                double temp = a[i + j * lda];           // U[i][j]
                for (int k = i + 1; k < j; ++k)         // U[i][k] * X[k][j]
                    temp += a[i + k * lda] * a[k + j * lda];
                a[i + j * lda] = -temp * a[j + j * lda] / a[i + i * lda];
            }
        }
    } else {
        // X = L^{-1}: column j forward-substitution; a[i][i] is still the
        // original diagonal for i > j, so divide by it.
        for (int j = 0; j < n; ++j) {
            if (!unit) {
                if (a[j + j * lda] == 0.0) { info = j + 1; return; }
                a[j + j * lda] = 1.0 / a[j + j * lda];
            }
            for (int i = j + 1; i < n; ++i) {
                double temp = a[i + j * lda];           // L[i][j]
                for (int k = j + 1; k < i; ++k)         // L[i][k] * X[k][j]
                    temp += a[i + k * lda] * a[k + j * lda];
                a[i + j * lda] = -temp * a[j + j * lda] / a[i + i * lda];
            }
        }
    }
}

void dtrtri(char uplo, char diag, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DTRTRI", 1); return; }
    dtrti2(uplo, diag, n, a, lda, info);
}

void dpotri(char uplo, int n, double* a, int lda, int& info) {
    info = 0;
    if (n < 0) { xerbla("DPOTRI", 1); return; }
    bool upper_ = upper(uplo);
    if (upper_) {
        dtrtri('U', 'N', n, a, lda, info);
        if (info != 0) return;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i <= j; ++i) {
                double sum = 0.0;
                for (int k = j; k < n; ++k) sum += a[i + k * lda] * a[j + k * lda];
                a[i + j * lda] = sum;
            }
    } else {
        dtrtri('L', 'N', n, a, lda, info);
        if (info != 0) return;
        for (int j = 0; j < n; ++j)
            for (int i = j; i < n; ++i) {
                double sum = 0.0;
                for (int k = 0; k <= j; ++k) sum += a[i + k * lda] * a[j + k * lda];
                a[i + j * lda] = sum;
            }
    }
}


// Fortran-ABI bridge used by powsq.cpp (static inline ddot forwards to ddot_).
extern "C" double ddot_(int* n, const double* dx, int* incx,
                        const double* dy, int* incy) {
    return ddot(*n, dx, *incx, dy, *incy);
}
