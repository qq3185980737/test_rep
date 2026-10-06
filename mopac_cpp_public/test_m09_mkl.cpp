// test_m09_mkl.cpp — verification for mkl_bits.cpp (MOPAC 2016 mkl_bits.F90).
// BLAS/LAPACK routines tested on small matrices; expected values computed
// independently in python (column-major Fortran layout).
#include <cstdio>
#include <cmath>
#include "mkl_bits.h"

static int g_checks = 0; static int g_fail = 0; static double g_err = 0.0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    double e = std::fabs(got - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1.0);
    if (e > g_err) g_err = e;
    if (e > tol) { ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.12g want %.12g\n", name, got, want); }
}

int main() {
    std::fprintf(stderr, "[t1] level-1 (ddot/daxpy/dcopy/dscal/dswap/dnrm2)\n");
    {
        double x[4] = {1.0, 2.0, 3.0, 4.0}, y[4] = {4.0, 3.0, 2.0, 1.0};
        chk_rel(ddot(4, x, 1, y, 1), 20.0, 1e-12, "ddot=20");
        double z[4] = {0.0, 0.0, 0.0, 0.0};
        daxpy(4, 2.0, x, 1, z, 1);
        chk_rel(z[3], 8.0, 1e-12, "daxpy z3=8");
        double w[4];
        dcopy(4, x, 1, w, 1);
        chk_rel(w[2], 3.0, 1e-12, "dcopy w2=3");
        dscal(4, 0.5, x, 1);
        chk_rel(x[1], 1.0, 1e-12, "dscal x1=1");
        dswap(4, x, 1, y, 1);
        chk_rel(x[0], 4.0, 1e-12, "dswap x0=4");
        chk_rel(y[0], 0.5, 1e-12, "dswap y0=0.5");
        chk_rel(dnrm2(3, y, 1), std::sqrt(0.25 + 1.0 + 2.25), 1e-12, "dnrm2");
    }

    std::fprintf(stderr, "[t2] dgemm N,N (A*B)\n");
    {
        // A=[[1,2],[3,4]], B=[[5,6],[7,8]] column-major
        double a[4] = {1, 3, 2, 4}, b[4] = {5, 7, 6, 8}, c[4] = {0, 0, 0, 0};
        dgemm('N', 'N', 2, 2, 2, 1.0, a, 2, b, 2, 0.0, c, 2);
        // C=[[19,22],[43,50]] colmajor {19,43,22,50}
        chk_rel(c[0], 19.0, 1e-12, "dgemm c0=19");
        chk_rel(c[1], 43.0, 1e-12, "dgemm c1=43");
        chk_rel(c[2], 22.0, 1e-12, "dgemm c2=22");
        chk_rel(c[3], 50.0, 1e-12, "dgemm c3=50");
        // beta path: C = 1.0*A*B + 2.0*C0
        double c0[4] = {1, 1, 1, 1};
        dgemm('N', 'N', 2, 2, 2, 1.0, a, 2, b, 2, 2.0, c0, 2);
        chk_rel(c0[0], 21.0, 1e-12, "dgemm beta c0=21");
    }

    std::fprintf(stderr, "[t3] dgemv / dger / dsyrk\n");
    {
        double a[4] = {1, 3, 2, 4}, x[2] = {1, 1}, y[2] = {0, 0};
        dgemv('N', 2, 2, 1.0, a, 2, x, 1, 0.0, y, 1);
        chk_rel(y[0], 3.0, 1e-12, "dgemv y0=3");
        chk_rel(y[1], 7.0, 1e-12, "dgemv y1=7");
        // dger: A += x*y^T, x={1,2}, y={1,1}
        double xg[2] = {1, 2}, yg[2] = {1, 1}, ag[4] = {1, 3, 2, 4};
        dger(2, 2, 1.0, xg, 1, yg, 1, ag, 2);
        chk_rel(ag[0], 2.0, 1e-12, "dger a0=2");
        chk_rel(ag[3], 6.0, 1e-12, "dger a3=6");
        // dsyrk upper: C = A*A^T -> [[5,11],[11,25]] upper colmajor {5,0,11,25}
        double ak[4] = {1, 3, 2, 4}, ck[4] = {0, 0, 0, 0};
        dsyrk('U', 'N', 2, 2, 1.0, ak, 2, 0.0, ck, 2);
        chk_rel(ck[0], 5.0, 1e-12, "dsyrk c0=5");
        chk_rel(ck[2], 11.0, 1e-12, "dsyrk c2=11");
        chk_rel(ck[3], 25.0, 1e-12, "dsyrk c3=25");
    }

    std::fprintf(stderr, "[t4] dgesv (Ax=b)\n");
    {
        // A=[[2,1],[1,3]] colmajor {2,1,1,3}, b={3,5} -> x={0.8,1.4}
        double a[4] = {2, 1, 1, 3}, b[2] = {3, 5};
        int ipiv[2], info = 0;
        dgesv(2, 1, a, 2, ipiv, b, 2, info);
        chk(info == 0, "dgesv info=0");
        chk_rel(b[0], 0.8, 1e-12, "dgesv x0=0.8");
        chk_rel(b[1], 1.4, 1e-12, "dgesv x1=1.4");
        // dgetrs standalone with a fresh factor
        double a2[4] = {2, 1, 1, 3}, b2[2] = {1, 2};
        int ip2[2], info2 = 0;
        dgetrf(2, 2, a2, 2, ip2, info2);
        dgetrs('N', 2, 1, a2, 2, ip2, b2, 2, info2);
        chk_rel(b2[0], 0.2, 1e-12, "dgetrs x0=0.2");
        chk_rel(b2[1], 0.6, 1e-12, "dgetrs x1=0.6");
    }

    std::fprintf(stderr, "[t5] dpotrf/dpotri (Cholesky inverse)\n");
    {
        // A=[[4,2],[2,3]] upper colmajor {4,0,2,3}
        double a[4] = {4, 0, 2, 3};
        int info = 0;
        dpotrf('U', 2, a, 2, info);
        chk(info == 0, "dpotrf info=0");
        chk_rel(a[0], 2.0, 1e-12, "dpotrf U11=2");
        chk_rel(a[2], 1.0, 1e-12, "dpotrf U12=1");
        chk_rel(a[3], std::sqrt(2.0), 1e-12, "dpotrf U22=sqrt2");
        dpotri('U', 2, a, 2, info);
        chk(info == 0, "dpotri info=0");
        chk_rel(a[0], 0.375, 1e-12, "dpotri inv11=0.375");
        chk_rel(a[2], -0.25, 1e-12, "dpotri inv12=-0.25");
        chk_rel(a[3], 0.5, 1e-12, "dpotri inv22=0.5");
    }

    std::fprintf(stderr, "[t6] dtrsm / dtrmm / dlauum\n");
    {
        // L=[[1,0],[0.5,1]] lower colmajor {1,0.5,0,1}; solve L*X=I -> X=L^-1={1,-0.5,0,1}
        double l[4] = {1, 0.5, 0, 1}, b[4] = {1, 0, 0, 1};
        dtrsm('L', 'L', 'N', 'N', 2, 2, 1.0, l, 2, b, 2);
        chk_rel(b[0], 1.0, 1e-12, "dtrsm b0=1");
        chk_rel(b[1], -0.5, 1e-12, "dtrsm b1=-0.5");
        chk_rel(b[3], 1.0, 1e-12, "dtrsm b3=1");
        // U=[[1,2],[0,3]] upper colmajor {1,0,2,3}; U*U^T = [[5,6],[6,9]]
        double u[4] = {1, 0, 2, 3};
        int uinfo = 0;
        dlauum('U', 2, u, 2, uinfo);
        chk(uinfo == 0, "dlauum info=0");
        chk_rel(u[0], 5.0, 1e-12, "dlauum u0=5");
        chk_rel(u[2], 6.0, 1e-12, "dlauum u2=6");
        chk_rel(u[3], 9.0, 1e-12, "dlauum u3=9");
        // dtrmm: L*B with unit diag; L=[[1,0],[0.5,1]], B=[[1,2],[3,4]] colmajor {1,3,2,4}
        double lm[4] = {1, 0.5, 0, 1}, bm[4] = {1, 3, 2, 4};
        dtrmm('L', 'L', 'N', 'U', 2, 2, 1.0, lm, 2, bm, 2);
        // L*B = [[1,2],[3.5,5]] -> colmajor {1,3.5,2,5}
        chk_rel(bm[0], 1.0, 1e-12, "dtrmm bm0=1");
        chk_rel(bm[1], 3.5, 1e-12, "dtrmm bm1=3.5");
        chk_rel(bm[2], 2.0, 1e-12, "dtrmm bm2=2");
        chk_rel(bm[3], 5.0, 1e-12, "dtrmm bm3=5");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_fail ? 1 : 0;
}
