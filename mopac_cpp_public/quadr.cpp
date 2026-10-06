// quadr.cpp — solve quadratic through three points.
#include "quadr.h"
void quadr(double f0, double f1, double f2, double x1, double x2,
           double& a, double& b, double& c) {
    c = (x2*(f1-f0) - x1*(f2-f0))/(x2*x1*(x1-x2));
    b = (f1 - f0 - c*x1*x1)/x1;
    a = f0;
}
