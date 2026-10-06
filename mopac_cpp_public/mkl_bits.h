// mkl_bits.h — C++ translation of MOPAC 2016 "mkl_bits.F90".
// Self-contained BLAS level-1/2/3 and LAPACK routines (column-major,
// Fortran data layout), mirroring the standard Netlib signatures.
#pragma once

// Level 1
double ddot(int n, const double* dx, int incx, const double* dy, int incy);
void dcopy(int n, const double* dx, int incx, double* dy, int incy);
void dscal(int n, double da, double* dx, int incx);
void daxpy(int n, double da, const double* dx, int incx, double* dy, int incy);
void dswap(int n, double* dx, int incx, double* dy, int incy);
double dnrm2(int n, const double* x, int incx);

// Level 2
void dgemv(char trans, int m, int n, double alpha, const double* a, int lda,
           const double* x, int incx, double beta, double* y, int incy);
void dger(int m, int n, double alpha, const double* x, int incx,
          const double* y, int incy, double* a, int lda);
void dtrmv(char uplo, char trans, char diag, int n, const double* a, int lda,
           double* x, int incx);

// Level 3
void dgemm(char transa, char transb, int m, int n, int k, double alpha,
           const double* a, int lda, const double* b, int ldb, double beta,
           double* c, int ldc);
void dsyrk(char uplo, char trans, int n, int k, double alpha, const double* a,
           int lda, double beta, double* c, int ldc);
void dtrmm(char side, char uplo, char transa, char diag, int m, int n,
           double alpha, const double* a, int lda, double* b, int ldb);
void dtrsm(char side, char uplo, char transa, char diag, int m, int n,
           double alpha, const double* a, int lda, double* b, int ldb);

// LAPACK linear equations
void dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb,
           int& info);
void dgetf2(int m, int n, double* a, int lda, int* ipiv, int& info);
void dgetrf(int m, int n, double* a, int lda, int* ipiv, int& info);
void dgetrs(char trans, int n, int nrhs, const double* a, int lda,
            const int* ipiv, double* b, int ldb, int& info);
void dlaswp(int n, double* a, int lda, int k1, int k2, const int* ipiv,
            int incx);

// LAPACK Cholesky / triangular inverse
void dpotf2(char uplo, int n, double* a, int lda, int& info);
void dpotrf(char uplo, int n, double* a, int lda, int& info);
void dpotri(char uplo, int n, double* a, int lda, int& info);
void dlauu2(char uplo, int n, double* a, int lda, int& info);
void dlauum(char uplo, int n, double* a, int lda, int& info);
void dtrti2(char uplo, char diag, int n, double* a, int lda, int& info);
void dtrtri(char uplo, char diag, int n, double* a, int lda, int& info);

// Error handler
void xerbla(const char* srname, int info);
