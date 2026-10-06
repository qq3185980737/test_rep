// density_for_MOZYME.cpp — C++ translation of MOPAC 2016.
#include "density_for_MOZYME.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "ijbo.h"
#include "molkst_C.h"

using namespace molkst_C;
using namespace MOZYME_C;

// Bare globals defined in MOZYME_C.cpp outside the namespace.

void density_for_MOZYME(std::vector<double>& p, int mode, int nclose_loc,
                        const std::vector<double>& partpin) {
    static bool first = true, prnt = false;
    if (first) { first = false; prnt = keywrd.find(" DIAG") != std::string::npos; }

    if (mode == 0) {
        for (int i = 0; i < (int)p.size(); ++i) p[i] = 0.0;
    } else if (mode == -1) {
        for (int i = 0; i < (int)p.size(); ++i) p[i] = -0.5 * partpin[i];
    } else {
        for (int i = 0; i < (int)p.size(); ++i) p[i] = 0.5 * partpin[i];
    }

    for (int i = 1; i <= nclose_loc; ++i) {
        int loopvar = ncocc[i];
        int ja = 0;
        if (lijbo) {
            for (int jj = nncf[i] + 1; jj <= nncf[i] + ncf[i]; ++jj) {
                int j = icocc[jj];
                int nj = iorbs[j];
                int ka = loopvar;
                for (int kk = nncf[i] + 1; kk <= nncf[i] + ncf[i]; ++kk) {
                    int k = icocc[kk];
                    if (j == k) {
                        int l = nijbo[j][k];
                        for (int j1 = 1; j1 <= nj; ++j1) {
                            double sum = cocc[ja + j1 + loopvar];
                            for (int k1 = 1; k1 <= j1; ++k1) {
                                int k2 = ka + k1;
                                l++;
                                p[l] += cocc[k2] * sum;
                            }
                        }
                    } else if (j > k && nijbo[j][k] >= 0) {
                        int l = nijbo[j][k];
                        for (int j1 = 1; j1 <= nj; ++j1) {
                            double sum = cocc[ja + j1 + loopvar];
                            for (int k1 = 1; k1 <= iorbs[k]; ++k1) {
                                int k2 = ka + k1;
                                l++;
                                p[l] += cocc[k2] * sum;
                            }
                        }
                    }
                    ka += iorbs[k];
                }
                ja += nj;
            }
        } else {
            for (int jj = nncf[i] + 1; jj <= nncf[i] + ncf[i]; ++jj) {
                int j = icocc[jj];
                int nj = iorbs[j];
                int ka = loopvar;
                for (int kk = nncf[i] + 1; kk <= nncf[i] + ncf[i]; ++kk) {
                    int k = icocc[kk];
                    int l = ijbo(j, k);
                    if (j == k) {
                        for (int j1 = 1; j1 <= nj; ++j1) {
                            double sum = cocc[ja + j1 + loopvar];
                            for (int k1 = 1; k1 <= j1; ++k1) {
                                int k2 = ka + k1;
                                l++;
                                p[l] += cocc[k2] * sum;
                            }
                        }
                    } else if (j > k && l >= 0) {
                        for (int j1 = 1; j1 <= nj; ++j1) {
                            double sum = cocc[ja + j1 + loopvar];
                            for (int k1 = 1; k1 <= iorbs[k]; ++k1) {
                                int k2 = ka + k1;
                                l++;
                                p[l] += cocc[k2] * sum;
                            }
                        }
                    }
                    ka += iorbs[k];
                }
                ja += nj;
            }
        }
    }

    double spinfa;
    if (mode == 0 || mode == 1) spinfa = 2.0;
    else if (mode == -1) spinfa = -2.0;
    else spinfa = 1.0;
    if (std::abs(spinfa - 1.0) > 1e-10)
        for (int i = 0; i < (int)p.size(); ++i) p[i] *= spinfa;
    if (prnt) {
        double sum = 0.0;
        for (int i = 1; i <= numat; ++i) {
            int j = (lijbo ? nijbo[i][i] : ijbo(i, i)) + 1;
            if (iorbs[i] > 0) sum += p[j];
            if (iorbs[i] == 4) sum += p[j + 2] + p[j + 5] + p[j + 9];
        }
        std::printf(" COMPUTED NUMBER OF ELECTRONS: %.0f\n", sum);
    }
}
