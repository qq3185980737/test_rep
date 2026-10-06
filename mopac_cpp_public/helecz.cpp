// helecz.cpp
#include "helecz.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace MOZYME_C;
using namespace common_arrays_C;
extern int ijbo(int, int);
extern "C" int ijbo_(int* i, int* j) { return ijbo(*i, *j); }
double helecz() {
    double ed = 0.0, ee = 0.0;
    for (int i = 1; i <= numat; ++i) {
        for (int j = 1; j <= i - 1; ++j) {
            int k = ijbo_(&i, &j);
            if (k >= 0) {
                int l = k + iorbs[i] * iorbs[j];
                for (int m = k + 1; m <= l; ++m)
                    ee += p[m] * (h[m] + f[m]);
            }
        }
        int k = ijbo_(&i, &i);
        for (int l1 = 1; l1 <= iorbs[i]; ++l1) {
            for (int l2 = 1; l2 <= l1 - 1; ++l2) {
                k++;
                ee += p[k] * (h[k] + f[k]);
            }
            k++;
            ed += p[k] * (h[k] + f[k]);
        }
    }
    ee += 0.5 * ed;
    return ee;
}
