// test_m12_batch.cpp — M12 线性代数/数值工具 数值对拍测试
#include <cstdio>
#include <cmath>
#include <vector>

#include "minv.h"
#include "mxm.h"
#include "mxv.h"
#include "mxmt.h"
#include "mtxm.h"
#include "dot.h"
#include "mult.h"
#include "mamult.h"
#include "mult33.h"
#include "mat33.h"
#include "supdot.h"
#include "osinv.h"
#include "rsp.h"
#include "linpack.h"
#include "schmit.h"
#include "schmib.h"
#include "interp.h"
#include "swap.h"
#include "molkst_C.h"
#include "funcon_C.h"
#include "symmetry_C.h"

using molkst_C::norbs;
using molkst_C::numcal;
using molkst_C::keywrd;

static int fails = 0;
static int passes = 0;

static void check(const char* name, bool ok) {
    if (ok) { ++passes; std::printf("  [PASS] %s\n", name); }
    else    { ++fails; std::printf("  [FAIL] %s\n", name); }
}
static bool close(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) <= tol * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

int main() {
    std::printf("M12 batch tests\n");

    // 1) minv: A = [[2,1,1],[1,2,1],[1,1,2]] (col-major), det=4, inv=0.25*[[3,-1,-1],[-1,3,-1],[-1,-1,3]]
    {
        double a[9] = {2,1,1, 1,2,1, 1,1,2};
        double d = 0;
        minv(a, 3, d);
        double a0[9] = {2,1,1, 1,2,1, 1,1,2};
        // A * inv = I
        bool ok = close(d, 4.0, 1e-8);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0;
                for (int k = 0; k < 3; ++k) s += a0[k + i * 3] * a[j + k * 3];
                ok = ok && close(s, (i == j) ? 1.0 : 0.0, 1e-8);
            }
        check("minv 3x3 inverse + det", ok);
    }

    // 2) mxm: C = A*B
    {
        const double A[6] = {1,2,3, 4,5,6};   // 2x3 col-major
        const double B[6] = {7,9, 8,10, 11,12}; // 3x2 col-major
        double C[4];
        mxm(A, 2, B, 3, C, 2);
        // C(1,1)=1*7+2*8+3*11=56? A row1=(1,2,3)?? col-major A(i,j)=A[i-1+(j-1)*2]
        // A(1,1)=1 A(2,1)=2 A(1,2)=3 A(2,2)=4 A(1,3)=5 A(2,3)=6
        // C(1,1)=1*7+3*8+5*11=86 ; C(2,1)=2*7+4*8+6*11=112 ; C(1,2)=1*9+3*10+5*12=99 ; C(2,2)=2*9+4*10+6*12=130
        bool ok = close(C[0], 74.0) && close(C[1], 98.0) && close(C[2], 103.0) && close(C[3], 136.0);
        check("mxm 2x3 * 3x2", ok);
    }

    // 3) mxv: y = A*x, A 2x3, x 3
    {
        const double A[6] = {1,2,3, 4,5,6};
        const double x[4] = {0, 1, 2, 3};
        double y[3];
        mxv(A, 2, x, 3, y);
        // y1 = 1*1+3*2+5*3 = 22 ; y2 = 2*1+4*2+6*3 = 28
        check("mxv 2x3 * x", close(y[1], 22.0) && close(y[2], 28.0));
    }

    // 4) mxmt: C = A*B', A 2x3, B 3x2
    {
        const double A[6] = {1,2,3, 4,5,6};
        const double B[6] = {1,2, 3,4, 5,6};  // 3x2 col-major
        double C[4];
        mxmt(A, 2, B, 3, C, 2);
        // C(i,j)=sum_k A(i,k)*B(j,k): C(1,1)=1*1+3*3+5*5=35 ; C(2,1)=2*1+4*3+6*5=44
        // C(1,2)=1*2+3*4+5*6=44 ; C(2,2)=2*2+4*4+6*6=56
        check("mxmt A*B'", close(C[0], 35.0) && close(C[1], 44.0) && close(C[2], 44.0) && close(C[3], 56.0));
    }

    // 5) mtxm: C = A'*B, A 3x2, B 3x2
    {
        const double A[6] = {1,2, 3,4, 5,6};  // 3x2
        const double B[6] = {7,8, 9,10, 11,12}; // 3x2
        double C[4];
        mtxm(A, 2, B, 3, C, 2);
        // C(1,1)=1*7+3*9+5*11=83 ; C(2,1)=2*7+4*9+6*11=116 ; C(1,2)=1*8+3*10+5*12=98 ; C(2,2)=2*8+4*10+6*12=124
        check("mtxm A'*B", close(C[0], 50.0) && close(C[1], 122.0) && close(C[2], 68.0) && close(C[3], 167.0));
    }

    // 6) dot
    {
        std::vector<double> x = {0, 1, 2, 3}, y = {0, 4, 5, 6};
        check("dot = 1*4+2*5+3*6 = 32", close(dot(x, y, 3), 32.0));
    }

    // 7) mult: VECS = S*C (col-major n x n)
    {
        const double S[4] = {2,0, 0,3};
        const double C[4] = {1,0, 0,1};
        double V[4];
        mult(C, S, V, 2);
        check("mult S*C", close(V[0], 2.0) && close(V[1], 0.0) && close(V[2], 0.0) && close(V[3], 3.0));
    }

    // 8) mamult: symmetric packed A*B + one*C, n=2
    {
        // A,B symmetric 2x2 packed (1-based): A={a11,a21,a22}, B={b11,b21,b22}
        const double A[4] = {0, 2.0, 1.0, 3.0};   // [[2,1],[1,3]]
        const double B[4] = {0, 1.0, 0.0, 2.0};   // [[1,0],[0,2]]
        double C[4] = {0, 0, 0, 0};
        mamult(A, B, C, 2, 1.0);
        // C = A*B = [[2,2],[1,6]] packed {2,1,6} at C[1],C[2],C[3]
        check("mamult 2x2", close(C[1], 2.0) && close(C[2], 1.0) && close(C[3], 6.0));
    }

    // 9) mult33: fmat = fmat*fmat with elem = identity
    {
        for (int i = 0; i < 21; ++i)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    symmetry_C::elem[r][c][i] = (r == c) ? 1.0 : 0.0;
        double fmat[9] = {1,0,0, 0,1,0, 0,0,1};
        mult33(fmat, 1);
        bool ok = true;
        for (int i = 0; i < 9; ++i) ok = ok && close(fmat[i], (i % 4 == 0) ? 1.0 : 0.0);
        check("mult33 identity*identity", ok);
    }

    // 10) mat33: C = A'*B*A with B=I => C=A'*A
    {
        const double A[9] = {1,0,0, 2,1,0, 3,4,1};  // row-major upper-triangular
        const double B[9] = {1,0,0, 0,1,0, 0,0,1};
        double C[9];
        mat33(A, B, C);
        // A'A(1,1)=1*1+2*2+3*3=14 ; A'A(2,2)=1+16=17 ; A'A(3,3)=1
        check("mat33 A'*A", close(C[0], 14.0) && close(C[4], 17.0) && close(C[8], 1.0));
    }

    // 11) supdot: s = H*g, H=[[2,1],[1,3]] packed {2,1,3}, g={1,1}
    {
        double h[4] = {0, 2.0, 1.0, 3.0};
        double g[3] = {0, 1.0, 1.0};
        double s[3] = {0, 0, 0};
        supdot(s, h, g, 2);
        // s1 = 2*1+1*1 = 3 ; s2 = 1*1+3*1 = 4
        check("supdot 2x2", close(s[1], 3.0) && close(s[2], 4.0));
    }

    // 12) osinv: symmetric inverse of [[2,1],[1,2]] => 0.5*[[2,-1],[-1,2]]
    {
        double a[4] = {2,1, 1,2};
        double d = 0;
        osinv(a, 2, d);
        bool ok = close(a[0], 2.0 / 3.0, 1e-8) && close(a[1], -1.0 / 3.0, 1e-8) &&
                  close(a[2], -1.0 / 3.0, 1e-8) && close(a[3], 2.0 / 3.0, 1e-8);
        check("osinv 2x2", ok);
    }

    // 13) linpack: dgefa+dgedi inverse vs analytic
    {
        double a[9] = {2,1,1, 1,2,1, 1,1,2};
        int ipvt[4];
        int info = 0;
        dgefa(a, 3, 3, ipvt, info);
        double det[2];
        double work[4];
        dgedi(a, 3, 3, ipvt, det, work, 10);
        bool ok = (info == 0) && close(det[0] * std::pow(10.0, det[1]), 4.0, 1e-6);
        check("linpack dgefa/dgedi det", ok);
    }

    // 14) rsp: [[2,1],[1,2]] eigenvalues {1,3}
    {
        double a[3] = {2, 1, 2};  // packed lower triangle
        double root[2];
        double vect[4];
        rsp(a, 2, root, vect);
        check("rsp eigenvalues 1,3",
              (close(root[0], 1.0, 1e-6) && close(root[1], 3.0, 1e-6)) ||
              (close(root[0], 3.0, 1e-6) && close(root[1], 1.0, 1e-6)));
    }

    // 15) schmit: orthogonalize columns {1,0} and {1,1}
    {
        double u[4] = {1,0, 1,1};  // col1={1,0}, col2={1,1}
        schmit(u, 2, 2);
        bool ok = close(u[0] * u[0] + u[1] * u[1], 1.0, 1e-8) &&
                  close(u[2] * u[2] + u[3] * u[3], 1.0, 1e-8) &&
                  close(u[0] * u[2] + u[1] * u[3], 0.0, 1e-8);
        check("schmit orthonormal", ok);
    }

    // 16) schmib
    {
        double u[4] = {1,0, 1,1};
        schmib(u, 2, 2);
        bool ok = close(u[0] * u[0] + u[1] * u[1], 1.0, 1e-8) &&
                  close(u[2] * u[2] + u[3] * u[3], 1.0, 1e-8) &&
                  close(u[0] * u[2] + u[1] * u[3], 0.0, 1e-8);
        check("schmib orthonormal", ok);
    }

    // 17) interp: mode=1 entry should return with mode=2 (no spline)
    {
        norbs = 2;
        numcal = 1;
        keywrd = "";
        double fp[3] = {1.0, 0.5, 2.0};
        double cp[4] = {1,0, 0,1};
        double theta[4] = {0,0,0,0};
        double vec_interp[4] = {0,0,0,0};
        double fock_interp[4] = {0,0,0,0};
        double p_interp[4] = {0,0,0,0};
        double h_interp[6] = {0,0,0,0,0,0};
        double vecl[4] = {0,0,0,0};
        double eold_ref = 0.0;
        int mode = 1;
        interp(1, 1, mode, -10.0, fp, cp, theta, vec_interp, fock_interp,
               p_interp, h_interp, vecl, eold_ref);
        check("interp mode1 -> mode2", mode == 2 && close(eold_ref, -10.0));
    }

    // 18) spline: parabola f=x^2 at -1,0,1 => xmin ~ 0
    {
        double x[4] = {0, -1.0, 0.0, 1.0};
        double f[4] = {0, 1.0, 0.0, 1.0};
        double df[4] = {0, -2.0, 0.0, 2.0};
        double xmin = 0;
        spline(x, f, df, 1.0, -1.0, xmin, 3);
        check("spline parabola xmin~0", close(xmin, 0.0, 1e-6));
    }

    // 19) swap: initialize psi with ifill<0, then swap virtual into occupied
    {
        double c2[6] = {1,0,0, 0,1,0};  // 3 rows x 2 cols
        int ifill = -1;
        swap(c2, 3, 3, 1, ifill);   // stores column 1 as psi/stdpsi, ifill -> 1
        bool ok = (ifill == 1);
        // second call with a fresh target; must not crash, no ASan
        ifill = 2;
        swap(c2, 3, 3, 1, ifill);
        // occupied column (col 0) must now hold a virtual column's content
        // (either original col 1 or a searched best-match column)
        ok = ok && close(c2[0], 1.0) && close(c2[1], 0.0);
        check("swap init + virtual swap smoke", ok);
    }

    std::printf("M12 batch: %d passed, %d failed\n", passes, fails);
    return fails == 0 ? 0 : 1;
}
