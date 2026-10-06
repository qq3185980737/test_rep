// esp_utilities.cpp — C++ translation of MOPAC 2016 "esp_utilities.F90".

#include "esp_utilities.h"

#include <algorithm>
#include <cmath>

void rys(double x, int nroots, double u[4], double w[4]) {
    // Rys quadrature — direct translation of esp_utilities.F90 "rys"
    // (Chebyshev asymptotic fits by x range).
    const double pie4 = 7.85398163397448e-01;
    const double r12 = 2.75255128608411e-01;
    const double w22 = 9.17517095361369e-02;
    const double r22 = 2.72474487139158e+00;
    double Ex, f1, y;
    for (int i = 0; i < 4; ++i) { u[i] = 0.0; w[i] = 0.0; }
    if (x > 33.0) {
        w[1] = std::sqrt(pie4 / x);
        if (nroots == 1) {
            u[1] = 0.5 / (x - 0.5);
            w[1] = std::sqrt(pie4 / x);
        } else {
            if (x > 40.0) {
                u[1] = r12 / (x - r12);
                u[2] = r22 / (x - r22);
                w[2] = w22 * w[1];
                w[1] = w[1] - w[2];
            } else {
                Ex = std::exp(-x);
                u[1] = (-8.78947307498880e-01 * x + 1.09243702330261e+01) * Ex +
                       r12 / (x - r12);
                u[2] = (-9.28903924275977e+00 * x + 8.10642367843811e+01) * Ex +
                       r22 / (x - r22);
                w[2] = (4.46857389308400e+00 * x - 7.79250653461045e+01) * Ex +
                       w22 * w[1];
                w[1] = w[1] - w[2];
            }
        }
    } else if (x > 15.0) {
        Ex = std::exp(-x);
        w[1] = (((1.9623264149430e-01 / x - 4.9695241464490e-01) / x -
                 6.0156581186481e-05) * Ex + std::sqrt(pie4 / x));
        f1 = (w[1] - Ex) / (2.0 * x);
        if (nroots == 1) {
            u[1] = f1 / (w[1] - f1);
        } else {
            u[1] = ((((-1.14906395546354e-06 * x + 1.76003409708332e-04) * x -
                      1.71984023644904e-02) * x - 1.37292644149838e-01) * x +
                    (-4.75742064274859e+01 / x + 9.21005186542857e+00) / x -
                    2.31080873898939e-02) * Ex + r12 / (x - r12);
            u[2] = (((3.64921633404158e-04 * x - 9.71850973831558e-02) * x -
                     4.02886174850252e+00) * x +
                    (-1.35831002139173e+02 / x - 8.66891724287962e+01) / x +
                    2.98011277766958e+00) * Ex + r22 / (x - r22);
            w[2] = ((f1 - w[1]) * u[1] + f1) * (1.0 + u[2]) / (u[2] - u[1]);
            w[1] = w[1] - w[2];
        }
    } else if (x > 10.0) {
        Ex = std::exp(-x);
        w[1] = ((((-1.8784686463512e-01 / x + 2.2991849164985e-01) / x -
                  4.9893752514047e-01) / x - 2.1916512131607e-05) * Ex +
                std::sqrt(pie4 / x));
        f1 = (w[1] - Ex) / (2.0 * x);
        if (nroots == 1) {
            u[1] = f1 / (w[1] - f1);
        } else {
            u[1] = ((((-1.01041157064226e-05 * x + 1.19483054115173e-03) * x -
                      6.73760231824074e-02) * x + 1.25705571069895e+00) * x +
                    (((-8.57609422987199e+03 / x + 5.91005939591842e+03) / x -
                      1.70807677109425e+03) / x + 2.64536689959503e+02) / x -
                    2.38570496490846e+01) * Ex + r12 / (x - r12);
            u[2] = (((3.39024225137123e-04 * x - 9.34976436343509e-02) * x -
                     4.22216483306320e+00) * x +
                    (((-2.08457050986847e+03 / x - 1.04999071905664e+03) / x +
                      3.39891508992661e+02) / x - 1.56184800325063e+02) / x +
                    8.00839033297501e+00) * Ex + r22 / (x - r22);
            w[2] = ((f1 - w[1]) * u[1] + f1) * (1.0 + u[2]) / (u[2] - u[1]);
            w[1] = w[1] - w[2];
        }
    } else if (x > 5.0) {
        Ex = std::exp(-x);
        w[1] = (((((((4.6897511375022e-01 / x - 6.9955602298985e-01) / x +
                     5.3689283271887e-01) / x - 3.2883030418398e-01) / x -
                    2.4645596956002e-01) / x - 4.9984072848436e-01) / x -
                  3.1501078774085e-06) * Ex + std::sqrt(pie4 / x));
        f1 = (w[1] - Ex) / (2.0 * x);
        if (nroots == 1) {
            u[1] = f1 / (w[1] - f1);
        } else {
            y = x - 7.5;
            u[1] = (((((((((((((-1.43632730148572e-16 * y +
                                2.38198922570405e-16) * y +
                               1.358319618800e-14) * y - 7.064522786879e-14) * y -
                              7.719300212748e-13) * y + 7.802544789997e-12) * y +
                             6.628721099436e-11) * y - 1.775564159743e-09) * y +
                            1.713828823990e-08) * y - 1.497500187053e-07) * y +
                           2.283485114279e-06) * y - 3.76953869614706e-05) * y +
                          4.74791204651451e-04) * y - 4.60448960876139e-03) * y +
                         3.72458587837249e-02;
            u[2] = ((((((((((((2.48791622798900e-14 * y - 1.36113510175724e-13) *
                             y - 2.224334349799e-12) * y + 4.190559455515e-11) * y -
                           2.222722579924e-10) * y - 2.624183464275e-09) * y +
                          6.128153450169e-08) * y - 4.383376014528e-07) * y -
                         2.49952200232910e-06) * y + 1.03236647888320e-04) * y -
                        1.44614664924989e-03) * y + 1.35094294917224e-02) * y -
                       9.53478510453887e-02) * y + 5.44765245686790e-01;
            w[2] = ((f1 - w[1]) * u[1] + f1) * (1.0 + u[2]) / (u[2] - u[1]);
            w[1] = w[1] - w[2];
        }
    } else if (x > 3.0) {
        y = x - 4.0;
        f1 = ((((((((((-2.62453564772299e-11 * y + 3.24031041623823e-10) * y -
                      3.614965656163e-09) * y + 3.760256799971e-08) * y -
                     3.553558319675e-07) * y + 3.022556449731e-06) * y -
                    2.290098979647e-05) * y + 1.526537461148e-04) * y -
                   8.81947375894379e-04) * y + 4.33207949514611e-03) * y -
                  1.75257821619926e-02) * y + 5.28406320615584e-02;
        w[1] = 2.0 * x * f1 + std::exp(-x);
        if (nroots == 1) {
            w[1] = 2.0 * x * f1 + std::exp(-x);
            u[1] = f1 / (w[1] - f1);
        } else {
            u[1] = ((((((((-4.11560117487296e-12 * y + 7.10910223886747e-11) * y -
                          1.73508862390291e-09) * y + 5.93066856324744e-08) * y -
                         9.76085576741771e-07) * y + 1.08484384385679e-05) * y -
                        1.12608004981982e-04) * y + 1.16210907653515e-03) * y -
                       9.89572595720351e-03) * y + 6.12589701086408e-02;
            u[2] = (((((((((-1.80555625241001e-10 * y + 5.44072475994123e-10) * y +
                           1.603498045240e-08) * y - 1.497986283037e-07) * y -
                          7.017002532106e-07) * y + 1.85882653064034e-05) * y -
                         2.04685420150802e-05) * y - 2.49327728643089e-03) * y +
                        3.56550690684281e-02) * y - 2.60417417692375e-01) * y +
                       1.12155283108289e+00;
            w[2] = ((f1 - w[1]) * u[1] + f1) * (1.0 + u[2]) / (u[2] - u[1]);
            w[1] = w[1] - w[2];
        }
    } else if (x > 1.0) {
        y = x - 2.0;
        f1 = ((((((((((-1.61702782425558e-10 * y + 1.96215250865776e-09) * y -
                      2.14234468198419e-08) * y + 2.17216556336318e-07) * y -
                     1.98850171329371e-06) * y + 1.62429321438911e-05) * y -
                    1.16740298039895e-04) * y + 7.24888732052332e-04) * y -
                   3.79490003707156e-03) * y + 1.61723488664661e-02) * y -
                  5.29428148329736e-02) * y + 1.15702180856167e-01;
        w[1] = 2.0 * x * f1 + std::exp(-x);
        if (nroots == 1) {
            u[1] = f1 / (w[1] - f1);
        } else {
            u[1] = (((((((((-6.36859636616415e-12 * y + 8.47417064776270e-11) * y -
                          5.152207846962e-10) * y - 3.846389873308e-10) * y +
                         8.472253388380e-08) * y - 1.85306035634293e-06) * y +
                        2.47191693238413e-05) * y - 2.49018321709815e-04) * y +
                       2.19173220020161e-03) * y - 1.63329339286794e-02) * y +
                      8.68085688285261e-02;
            u[2] = (((((((((1.45331350488343e-10 * y + 2.07111465297976e-09) * y -
                          1.878920917404e-08) * y - 1.725838516261e-07) * y +
                         2.247389642339e-06) * y + 9.76783813082564e-06) * y -
                        1.93160765581969e-04) * y - 1.58064140671893e-03) * y +
                       4.85928174507904e-02) * y - 4.30761584997596e-01) * y +
                      1.80400974537950e+00;
            w[2] = ((f1 - w[1]) * u[1] + f1) * (1.0 + u[2]) / (u[2] - u[1]);
            w[1] = w[1] - w[2];
        }
    } else if (x > 3.0e-7) {
        f1 = ((((((((-8.36313918003957e-08 * x + 1.21222603512827e-06) * x -
                    1.15662609053481e-05) * x + 9.25197374512647e-05) * x -
                   6.40994113129432e-04) * x + 3.78787044215009e-03) * x -
                  1.85185172458485e-02) * x + 7.14285713298222e-02) * x -
                 1.99999999997023e-01) * x + 3.33333333333318e-01;
        w[1] = 2.0 * x * f1 + std::exp(-x);
        if (nroots == 1) {
            u[1] = f1 / (w[1] - f1);
        } else {
            u[1] = (((((((-2.35234358048491e-09 * x + 2.49173650389842e-08) * x -
                         4.558315364581e-08) * x - 2.447252174587e-06) * x +
                        4.743292959463e-05) * x - 5.33184749432408e-04) * x +
                       4.44654947116579e-03) * x - 2.90430236084697e-02) * x +
                      1.30693606237085e-01;
            u[2] = (((((((-2.47404902329170e-08 * x + 2.36809910635906e-07) * x +
                         1.835367736310e-06) * x - 2.066168802076e-05) * x -
                        1.345693393936e-04) * x - 5.88154362858038e-05) * x +
                       5.32735082098139e-02) * x - 6.37623643056745e-01) * x +
                      2.86930639376289e+00;
            w[2] = ((f1 - w[1]) * u[1] + f1) * (1.0 + u[2]) / (u[2] - u[1]);
            w[1] = w[1] - w[2];
        }
    } else {
        // NOT verified in the Fortran source
        if (nroots == 1) {
            u[1] = 0.5 - x / 5.0;
            w[1] = 1.0 - x / 3.0;
        } else {
            u[1] = 1.30693606237085e-01 - 2.90430236082028e-02 * x;
            u[2] = 2.86930639376291e+00 - 6.37623643058102e-01 * x;
            w[1] = 6.52145154862545e-01 - 1.22713621927067e-01 * x;
            w[2] = 3.47854845137453e-01 - 2.10619711404725e-01 * x;
        }
    }
}

// Gauss-Hermite quadrature (21-point table from esp_utilities.F90 "vint").
// root(k) has weight w(k); k = 1..21 (1-based).
static const double vint_roots[22] = {
    0.0,
    0.0, -0.7071067811865, 0.7071067811865,
    -1.224744871391, 0.0, 1.2247448713915,
    -1.6506801238857, -0.52464762327529, 0.52464762327529, 1.6506801238857,
    -2.020182870456, -0.9585724646138, 0.0, 0.9585724646138, 2.020182870456,
    -2.350604973674, -1.335849074014, -0.436077411928, 0.436077411928,
    1.335849074014, 2.350604973674};
static const double vint_wts[22] = {
    0.0,
    1.772453850905,
    0.886226925452, 0.886226925452,
    0.295408975150, 1.18163590060, 0.295408975150,
    8.131283544725e-02, 0.8049140900055, 0.8049140900055, 8.131283544725e-02,
    1.995324205905e-02, 3.936193231522e-01, 9.453087204829e-01, 3.936193231522e-01, 1.995324205905e-02,
    4.530009905509e-03, 1.570673203229e-01, 7.246295952244e-01, 7.246295952244e-01,
    1.570673203229e-01, 4.530009905509e-03};
static const int vint_start[7] = {0, 0, 1, 3, 6, 10, 15};
static const int vint_end[7] = {0, 0, 2, 5, 9, 14, 20};

void vint(double& xint, double& yint, double& zint, int ni, int nj,
          double x0, double y0, double z0,
          double xi, double yi, double zi,
          double xj, double yj, double zj, double t) {
    // Discrete Gauss-Hermite quadrature, direct translation of
    // esp_utilities.F90 "vint": shell products evaluated at the nodes.
    xint = 0.0; yint = 0.0; zint = 0.0;
    int npt = (ni + nj) / 2 + 1;
    if (npt < 1) npt = 1;
    if (npt > 6) npt = 6;
    int lo = vint_start[npt] + 1;
    int hi = vint_end[npt] + 1;
    for (int i = lo; i <= hi; ++i) {
        double px = 1.0, py = 1.0, pz = 1.0;
        double dum = vint_roots[i] * t;
        double ptx = dum + x0, pty = dum + y0, ptz = dum + z0;
        double ax = ptx - xi, ay = pty - yi, az = ptz - zi;
        double bx = ptx - xj, by = pty - yj, bz = ptz - zj;
        for (int j = 2; j <= ni; ++j) { px *= ax; py *= ay; pz *= az; }
        for (int j = 2; j <= nj; ++j) { px *= bx; py *= by; pz *= bz; }
        double dw = vint_wts[i];
        xint += dw * px;
        yint += dw * py;
        zint += dw * pz;
    }
}

void evec(std::vector<float>& aVector, double x, double y, double z,
          const std::vector<std::vector<double>>& coord, int numat) {
    int j = 0;
    for (int i = 1; i <= numat; ++i) {
        float u = (float)(x - coord[0][i]);
        float v = (float)(y - coord[1][i]);
        float w = (float)(z - coord[2][i]);
        float uu = u*u, vv = v*v, ww = w*w;
        float shellR2 = std::max(uu + vv + ww, 1e-2f);
        float shellR2i = 1.0f / (shellR2 + 1e-7f);
        float shellRi = std::sqrt(shellR2i);
        float shellR3i = shellRi * shellR2i;
        aVector[j + 0] = shellRi;
        aVector[j + 1] = u * shellR3i;
        aVector[j + 2] = v * shellR3i;
        aVector[j + 3] = w * shellR3i;
        aVector[j + 4] = shellR2i;
        aVector[j + 5] = shellR3i;
        aVector[j + 6] = shellR2i * shellR2i;
        j += 7;
    }
}

void saxpy(int n, float sa, const std::vector<float>& sx, int incx,
           std::vector<float>& sy, int incy) {
    if (n <= 0 || sa == 0.0f) return;
    int ix = 1, iy = 1;
    if (incx < 0) ix = (1 - n) * incx + 1;
    if (incy < 0) iy = (1 - n) * incy + 1;
    for (int k = 0; k < n; ++k) {
        sy[iy + k * incy] += sa * sx[ix + k * incx];
    }
}

float sdot(int n, const std::vector<float>& sx, int incx,
           const std::vector<float>& sy, int incy) {
    if (n <= 0) return 0.0f;
    int ix = 1, iy = 1;
    if (incx < 0) ix = (1 - n) * incx + 1;
    if (incy < 0) iy = (1 - n) * incy + 1;
    float s = 0.0f;
    for (int k = 0; k < n; ++k) s += sx[ix + k * incx] * sy[iy + k * incy];
    return s;
}

void sscal(int n, float da, std::vector<float>& dx, int incx) {
    if (n <= 0) return;
    int ix = 1;
    if (incx < 0) ix = (1 - n) * incx + 1;
    for (int k = 0; k < n; ++k) dx[ix + k * incx] *= da;
}

float snrm2(int n, const std::vector<float>& sx, int incx) {
    if (n <= 0) return 0.0f;
    int ix = 1;
    if (incx < 0) ix = (1 - n) * incx + 1;
    double sum = 0;
    for (int k = 0; k < n; ++k) { double v = sx[ix + k * incx]; sum += v * v; }
    return (float)std::sqrt(sum);
}

void sswap(int n, std::vector<float>& sx, int incx,
           std::vector<float>& sy, int incy) {
    if (n <= 0) return;
    int ix = 1, iy = 1;
    if (incx < 0) ix = (1 - n) * incx + 1;
    if (incy < 0) iy = (1 - n) * incy + 1;
    for (int k = 0; k < n; ++k) {
        float t = sx[ix + k * incx];
        sx[ix + k * incx] = sy[iy + k * incy];
        sy[iy + k * incy] = t;
    }
}

void scopy(int n, const std::vector<float>& sx, int incx,
           std::vector<float>& sy, int incy) {
    if (n <= 0) return;
    int ix = 1, iy = 1;
    if (incx < 0) ix = (1 - n) * incx + 1;
    if (incy < 0) iy = (1 - n) * incy + 1;
    for (int k = 0; k < n; ++k) sy[iy + k * incy] = sx[ix + k * incx];
}
