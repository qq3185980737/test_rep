// bfn.cpp — C++ translation of MOPAC 2016 "bfn.F90".

#include "bfn.h"

#include <cmath>

#include "overlaps_C.h"

void bfn(double x, double bf[14]) {
    int k = 12, io = 0;
    double absx = std::fabs(x);
    int last;

    if (absx <= 3.0) {
        if (absx > 2.0) {
            last = 15;
        } else if (absx > 1.0) {
            last = 12;
        } else if (absx > 0.5) {
            last = 7;
        } else if (absx <= 1.0e-6) {
            // small-x limit (label 90)
            for (int i = io; i <= k; ++i)
                bf[i + 1] = (2.0 * ((i + 1) % 2)) / (i + 1.0);
            return;
        } else {
            last = 6;
        }
        // label 60: series expansion
        for (int i = io; i <= k; ++i) {
            double y = 0.0;
            for (int m = io; m <= last; ++m) {
                double xf = (m != 0) ? overlaps_C::fact[m] : 1.0;
                y += std::pow(-x, m) * (2 * ((m + i + 1) % 2)) / (xf * (m + i + 1));
            }
            bf[i + 1] = y;
        }
        return;
    }

    // |x| > 3: recurrence from exp
    double expx = std::exp(x);
    double expmx = 1.0 / expx;
    bf[1] = (expx - expmx) / x;
    for (int i = 1; i <= k; ++i) {
        double sign = (i % 2 == 0) ? 1.0 : -1.0;
        bf[i + 1] = (i * bf[i] + sign * expx - expmx) / x;
    }
}
