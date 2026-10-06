// dang.cpp — C++ translation of MOPAC 2016 "dang.F90".
#include "dang.h"

#include <algorithm>
#include <cmath>

void dang(double& a1, double& a2, double& b1, double& b2, double& rcos) {
    const double zero = 1e-6;
    if (std::abs(a1) >= zero || std::abs(a2) >= zero) {
        if (std::abs(b1) >= zero || std::abs(b2) >= zero) {
            double anorm = 1.0 / std::sqrt(a1 * a1 + a2 * a2);
            double bnorm = 1.0 / std::sqrt(b1 * b1 + b2 * b2);
            a1 *= anorm; a2 *= anorm;
            b1 *= bnorm; b2 *= bnorm;
            double sinth = a1 * b2 - a2 * b1;
            double costh = std::min(1.0, a1 * b1 + a2 * b2);
            costh = std::max(-1.0, costh);
            rcos = std::acos(costh);
            if (std::abs(rcos) >= 4e-5) {
                if (sinth > 0.0) rcos = 6.28318530717959 - rcos;
                rcos = -rcos;
                return;
            }
        }
    }
    rcos = 0.0;
}
