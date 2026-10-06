// eimp.cpp — C++ translation of MOPAC 2016 "eimp.F90".

#include "eimp.h"

#include "ijbo.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;


void eimp() {
    for (int i = 1; i <= numat; ++i)
        for (int j = 1; j < i; ++j) {
            int k = ijbo(i, j);
            if (k >= 0) {
                int l = MOZYME_C::iorbs[i] * MOZYME_C::iorbs[j];
                if (l != 0) {
                    l = k + l;
                    double sum = 0.0;
                    for (int m = k + 1; m <= l; ++m) sum += f[m] * f[m];
                    p[k + 1] = sum;
                }
            }
        }
}
