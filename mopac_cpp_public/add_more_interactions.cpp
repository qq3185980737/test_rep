// add_more_interactions.cpp — C++ translation of MOPAC 2016
// "add_more_interactions.F90" (127 lines, single subroutine).
//
// 1-based Fortran indexing is kept for the packed arrays (index 0 is
// padding).  The Fortran "copy to temp, deallocate, reallocate, copy back"
// pattern maps to std::vector::resize, which preserves elements 1..i and
// value-initialises the new tail to zero.

#include "add_more_interactions.h"

#include <cmath>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "iter_C.h"
#include "molkst_C.h"

// External routines (defined in not-yet-ported files).
void fillij(bool flag);
void memory_error(const char* where);

void add_more_interactions() {
    int i;
    static int imol = 0;  // save imol

    // Method is limited to direct SCF.
    if (!MOZYME_C::direct || !MOZYME_C::semidr || !MOZYME_C::lijbo) return;
    if (imol != molkst_C::numcal) {
        imol = molkst_C::numcal;
        return;  // first time this job is called: do nothing
    }

    i = molkst_C::mpack;
    fillij(false);

    if (molkst_C::mpack > i) {
        // Some arrays are now too small — grow them by 20%, preserving contents.
        int j = (int)std::lround(molkst_C::mpack * 1.2);

        // grow(v, i, j): preserve elements 1..i, zero the new tail i+1..j.
        auto grow = [j](std::vector<double>& v) {
            v.resize((size_t)j + 1, 0.0);
        };

        try {
            grow(iter_C::pold);
            grow(common_arrays_C::p);
            grow(common_arrays_C::h);
            grow(common_arrays_C::f);
            if (MOZYME_C::rapid) {
                grow(MOZYME_C::partp);
                grow(MOZYME_C::parth);
                grow(MOZYME_C::partf);
            }
        } catch (const std::bad_alloc&) {
            memory_error("add_more_interactions");
            return;
        }
    }
}
