// test_add_more_interactions.cpp — verify the array-grow logic.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "iter_C.h"
#include "molkst_C.h"

using common_arrays_C::f;
using common_arrays_C::h;
using common_arrays_C::p;

// fillij bump: in the test it grows mpack so that add_more_interactions acts.
static int g_bump = 0;
void fillij(bool) { molkst_C::mpack += g_bump; }
static bool mem_err_called = false;
void memory_error(const char*) { mem_err_called = true; }

#include "add_more_interactions.cpp"

int main() {
    using iter_C::pold;

    // Enable the direct-SCF path.
    MOZYME_C::direct = true;
    MOZYME_C::semidr = true;
    MOZYME_C::lijbo = true;
    MOZYME_C::rapid = true;

    const int old_i = 100;
    molkst_C::mpack = old_i;
    molkst_C::numcal = 7;

    // Fill arrays (1-based; index 0 padding).  Seed old elements.
    auto seed = [&](std::vector<double>& v) {
        v.assign(old_i + 1, 0.0);
        for (int k = 1; k <= old_i; ++k) v[k] = (double)k;
    };
    seed(pold); seed(p); seed(h); seed(f);
    seed(MOZYME_C::partp); seed(MOZYME_C::parth); seed(MOZYME_C::partf);

    // First call: imol != numcal -> immediate no-op.
    add_more_interactions();
    bool pass1 = (pold.size() == (size_t)old_i + 1);
    std::printf("first-call no-op: size=%zu  %s\n", pold.size(), pass1 ? "PASS" : "FAIL");

    // Second call: fillij grows mpack by 50 -> mpack=150, j=Nint(180)=180.
    g_bump = 50;
    add_more_interactions();
    int expect_j = (int)std::lround(150 * 1.2);  // 180
    bool ok = (int)pold.size() == expect_j + 1;
    for (int k = 1; k <= old_i && ok; ++k) ok = pold[k] == (double)k;       // preserved
    for (int k = old_i + 1; k <= expect_j && ok; ++k) ok = pold[k] == 0.0;  // zeroed
    bool part_ok = (int)MOZYME_C::partp.size() == expect_j + 1;
    std::printf("grow: old_i=%d new_mpack=%d j=%d pold_size=%zu tail_zero=%d partp_ok=%d  %s\n",
                old_i, molkst_C::mpack, expect_j, pold.size(),
                (int)ok, (int)part_ok, (ok && part_ok) ? "PASS" : "FAIL");

    // Disabled path: not direct -> no change.
    MOZYME_C::direct = false;
    size_t sz_before = p.size();
    add_more_interactions();
    bool pass3 = (p.size() == sz_before);
    std::printf("disabled path no-op: %s\n", pass3 ? "PASS" : "FAIL");

    return (pass1 && ok && part_ok && pass3) ? 0 : 1;
}
