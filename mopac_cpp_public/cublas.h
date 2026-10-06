// cublas.h — C++ bindings for the external CUBLAS/CUDA C API used by MOPAC.
// Mirrors mod_calls_cublas.F90. Host-side declarations; the implementations
// live in CUDA/CUBLAS runtime libraries and are linked in at build time.
#pragma once

extern "C" {

// BLAS-1 vector operations.
void call_asum_cublas(int n, double* vecx, int incx, double* res);
void call_axpy_cublas(int n, double alpha, double* vecx, int incx,
                      double* vecy, int incy);
void call_copy_cublas(int n, double* vecx, int incx, double* vecy, int incy);
void call_dot_cublas(int n, double* vecx, int incx, double* vecy, int incy,
                     double* res);
void call_scal_cublas(int n, double alpha, double* vecx, int incx);
void call_swap_cublas(int n, double* vecx, int incx, double* vecy, int incy);
void call_rot_cublas(int n, double* x, int incx, double* y, int incy,
                     double c, double s);
void call_iamax_cublas(int n, double* x, int incx, int* res);
void call_iamin_cublas(int n, double* x, int incx, int* res);

// BLAS-2 / BLAS-3 matrix operations.
void call_gemm_cublas(char tra, char trb, int m, int n, int k, double alpha,
                      double* a, int lda, double* b, int ldb, double beta,
                      double* c, int ldc);
void call_gemm_cublas_mgpu(char tra, char trb, int m, int n, int k,
                           double alpha, double* a, int lda, double* b,
                           int ldb, double beta, double* c, int ldc);
void call_gemm_cublas_thrust(char tra, char trb, int m, int n, int k,
                             double alpha, double* a, int lda, double* b,
                             int ldb, double beta, double* c, int ldc);
void phigemm(char tra, char trb, int m, int n, int k, double alpha,
             double* a, int lda, double* b, int ldb, double beta, double* c,
             int ldc);
void call_gemv_cublas(char tra, int m, int n, double alpha, double* a,
                      int lda, double* x, int incx, double beta, double* y,
                      int incy);
void call_ger_cublas(int m, int n, double alpha, double* x, int incx,
                     double* y, int incy, double* a, int lda);
void call_trmm_cublas(char side, char uplo, char transa, char diag, int m,
                      int n, double alpha, double* a, int lda, double* b,
                      int ldb);
void call_trmv_cublas(char uplo, char transa, char diag, int n, double* a,
                      int lda, double* x, int incx);
void call_trsm_cublas(char side, char uplo, char transa, char diag, int m,
                      int n, double alpha, double* a, int lda, double* b,
                      int ldb);
void call_syrk_cublas(char uplo, char trans, int n, int k, double alpha,
                      double* a, int lda, double beta, double* c, int ldc);
void call_syrk_cublas_thrust(char uplo, char trans, int n, int k,
                             double alpha, double* a, int lda, double beta,
                             double* c, int ldc);

// CUDA device context.
void create_handle();
void destroy_handle();
}
