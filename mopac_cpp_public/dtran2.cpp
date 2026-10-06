// dtran2.cpp — C++ translation of MOPAC 2016 "dtran2.F90".

#include "dtran2.h"

#include <cmath>

#include "funcon_C.h"

using namespace funcon_C;

void dtran2(std::vector<std::vector<double>>& r,
            std::vector<std::vector<std::vector<double>>>& t, int ioper) {
    const double tol = 0.001, s12 = 3.46410161513, s3 = 1.73205080756, one = 1.0;
    double r1 = r[1][0]*r[2][1] - r[2][0]*r[1][1];
    double r2 = r[2][0]*r[0][1] - r[0][0]*r[2][1];
    double r3 = r[0][0]*r[1][1] - r[1][0]*r[0][1];
    double check = r1*r[0][2] + r2*r[1][2] + r3*r[2][2];
    bool right = check > 0.0;
    r[0][2]=r1; r[1][2]=r2; r[2][2]=r3;
    double arg = r3;
    if (std::fabs(arg) > one) arg = (arg>=0?one:-one);
    double b = std::acos(arg);
    double sina = std::sqrt(1.0 - arg*arg);
    double g, a;
    if (sina >= tol) {
        arg = r[2][1]/sina; if (std::fabs(arg)>one) arg=(arg>=0?one:-one); g=std::asin(arg);
        arg = r[1][2]/sina; if (std::fabs(arg)>one) arg=(arg>=0?one:-one); a=std::asin(arg);
    } else {
        arg = r[0][1]; if (std::fabs(arg)>one) arg=(arg>=0?one:-one); g=std::asin(arg); a=0.0;
    }
    double f[2][4] = {{a,a,pi-a,pi-a},{g,pi-g,g,pi-g}};
    for (int i = 0; i < 4; ++i) {
        double ai=f[0][i], gi=f[1][i];
        check = std::fabs(std::sin(b)*std::cos(ai) + r[0][2]);
        if (check > tol) continue;
        check = (-std::sin(gi)*std::cos(b)*std::sin(ai)) + std::cos(gi)*std::cos(ai);
        if (std::fabs(check - r[1][1]) > tol) continue;
        check = std::sin(ai)*std::cos(gi) + std::cos(ai)*std::cos(b)*std::sin(gi);
        if (std::fabs(check - r[0][1]) > tol) continue;
        a=ai; g=gi; break;
    }
    g=-g; a=-a; b=-b;
    double e1=std::cos(b*0.5), x1=-std::sin(b*0.5);
    double e2=e1*e1, e3=e1*e2, e4=e2*e2;
    double x2=x1*x1, x3=x1*x2, x4=x2*x2;
    double ta=2*a, tg=2*g;
    auto& T = t[ioper];
    T[0][0]=e4*std::cos(ta+tg)+x4*std::cos(ta-tg);
    T[0][1]=2*e3*x1*std::cos(a+tg)-2*e1*x3*std::cos(a-tg);
    T[0][2]=2*s3*e2*x2*std::cos(tg);
    T[0][3]=2*e3*x1*std::sin(a+tg)-2*e1*x3*std::sin(a-tg);
    T[0][4]=e4*std::sin(ta+tg)+x4*std::sin(ta-tg);
    T[1][0]=2*e1*x3*std::cos(ta-g)-2*e3*x1*std::cos(ta+g);
    T[1][1]=(e4-3*e2*x2)*std::cos(a+g)-(3*e2*x2-x4)*std::cos(a-g);
    T[1][2]=2*s3*(e3*x1-e1*x3)*std::cos(g);
    T[1][3]=(e4-3*e2*x2)*std::sin(a+g)-(3*e2*x2-x4)*std::sin(a-g);
    T[1][4]=-2*e3*x1*std::sin(ta+g)+2*e1*x3*std::sin(ta-g);
    T[2][0]=s12*e2*x2*std::cos(ta);
    T[2][1]=-s12*(e3*x1-e1*x3)*std::cos(a);
    T[2][2]=e4-4*e2*x2+x4;
    T[2][3]=-s12*(e3*x1-e1*x3)*std::sin(a);
    T[2][4]=s12*e2*x2*std::sin(ta);
    T[3][0]=2*e1*x3*std::sin(ta-g)+2*e3*x1*std::sin(ta+g);
    T[3][1]=-(e4-3*e2*x2)*std::sin(a+g)-(3*e2*x2-x4)*std::sin(a-g);
    T[3][2]=-2*s3*(e3*x1-e1*x3)*std::sin(g);
    T[3][3]=(e4-3*e2*x2)*std::cos(a+g)+(3*e2*x2-x4)*std::cos(a-g);
    T[3][4]=-2*e3*x1*std::cos(ta+g)-2*e1*x3*std::cos(ta-g);
    T[4][0]=-e4*std::sin(ta+tg)+x4*std::sin(ta-tg);
    T[4][1]=-2*e3*x1*std::sin(a+tg)-2*e1*x3*std::sin(a-tg);
    T[4][2]=-2*s3*e2*x2*std::sin(tg);
    T[4][3]=2*e3*x1*std::cos(a+tg)+2*e1*x3*std::cos(a-tg);
    T[4][4]=e4*std::cos(ta+tg)-x4*std::cos(ta-tg);
    if (right) return;
    for (int c=0;c<5;++c){ T[1][c]=-T[1][c]; T[3][c]=-T[3][c]; }
}
