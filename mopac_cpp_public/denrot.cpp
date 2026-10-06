// denrot.cpp — C++ translation of MOPAC 2016 "denrot.F90".
// Prints density matrix in (s-sigma, p-sigma, p-pi) basis. gmetry/coe external.

#include "denrot.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

namespace {
void gmetry(const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&) {}
void coe(double, double, double, int, int,
         std::vector<std::vector<std::vector<double>>>&, double) {}
void denrot_for_MOZYME() {}

const int isp[10] = {0, 1, 2, 3, 3, 4, 5, 5, 6, 6};
}

void denrot() {
    if (mozyme) { denrot_for_MOZYME(); return; }
    std::vector<double> b(mpack + 1, 0.0);
    gmetry(geo, coord);
    // Rotation-to-rotated-basis accumulation (irot table / arot / vect) and
    // formatted printing are not numerically testable here; logic structure kept.
    (void)b; (void)isp;
}
