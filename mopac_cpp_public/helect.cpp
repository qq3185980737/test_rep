#include <vector>
// helect.cpp
#include "helect.h"
// Fortran packed lower-triangle arrays are 1-based; p/h/f[0] is padding.
double helect(int n, const std::vector<double>& p, const std::vector<double>& h, const std::vector<double>& f) { return helect(n, p.data(), h.data(), f.data()); }
double helect(int n, const double* p, const double* h, const double* f) {
    double ed = 0.0, ee = 0.0;
    int k = 0;  // Fortran k starts at 1 after first increment
    int nn = n + 1;
    for (int i = 2; i <= nn; ++i) {
        ++k;  // now Fortran k = 1, 2, ...
        int jj = i - 1;
        ed += p[k] * (h[k] + f[k]);
        if (i == nn) continue;
        if (jj > 0) {
            double s = 0.0;
            for (int t = 1; t <= jj; ++t) s += p[k+t] * (h[k+t] + f[k+t]);
            ee += s;
            k = jj + k;
        }
    }
    ee += 0.5 * ed;
    return ee;
}
