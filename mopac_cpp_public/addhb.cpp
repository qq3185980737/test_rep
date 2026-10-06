// addhb.cpp — C++ translation of MOPAC 2016 "addhb.F90"
//
// Original Fortran:
//
//   subroutine addhb (nocc1, nvir1, idiagg, nij, nhb)
//     use molkst_C, only : numat, norbs
//     use common_arrays_C, only: eigs, f
//     integer, intent (in) :: idiagg, nhb, nocc1, nvir1
//     integer, intent (out) :: nij
//     integer :: alloc_stat
//     logical, dimension (:), allocatable :: latom_loc
//     double precision, dimension(:), allocatable :: storei, storej
//     integer, allocatable, dimension(:) :: iused
//     double precision, dimension (4) :: hblims
//     data hblims / 1.d0, 0.1d0, 0.01d0, 0.001d0 /
//     allocate (latom_loc(numat), iused(Max(norbs, numat)),
//               storei(norbs), storej(norbs), stat=alloc_stat)
//     if (alloc_stat /= 0) then
//       call memory_error ("addhb")
//       return
//     end if
//     call hbonds (f, nocc1, nvir1, iused, nij, hblims(nhb))
//     if (nij /= 0) then
//       call diagg2 (nocc1, nvir1, eigs(nocc1+1:), iused, latom_loc,
//                    nij, idiagg, storei, storej)
//     end if
//     deallocate (latom_loc, iused, storei, storej)
//     if (alloc_stat /= 0) then
//       call memory_error ("addhb")
//     end if
//   end subroutine addhb
//
// Translation notes:
//  * Fortran 1-based arrays are mirrored with an unused slot 0; Fortran
//    array sections are passed as base + (start_index - 1).
//  * Allocation status is emulated by catching std::bad_alloc.
//  * The trailing alloc_stat check is a latent Fortran bug (deallocate has
//    no stat=); it is preserved as a documented dead check, never fires.

#include "addhb.h"

#include <algorithm>
#include <new>
#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"

// ---- Forward declarations: routines not yet ported ----------------------
// memory_error(name): fatal allocation-error handler.
// hbonds(fao, nocc, nvir, iused, nij_loc, cutoff): identify H-bonds to make.
// diagg2(nocc, nvir, eigv, iused, latoms, nij, idiagg, storei, storej):
//   2x2 Jacobi annihilation; eigv points at eigs(nocc1+1).
void memory_error(const char* name);
void hbonds(const double* fao, int nocc, int nvir, int* iused,
            int& nij_loc, double cutoff);
void diagg2(int nocc, int nvir, const double* eigv, int* iused,
            char* latoms, int nij, int idiagg, double* storei,
            double* storej);

void addhb(int nocc1, int nvir1, int idiagg, int& nij, int nhb) {
    // hblims(1..4) = 1, 0.1, 0.01, 0.001  (slot 0 is padding for 1-based use)
    static const double hblims[5] = {0.0, 1.0, 0.1, 0.01, 0.001};

    int alloc_stat = 0;

    // Temporaries, Fortran 1-based (index 0 unused).
    std::vector<char> latom_loc;        // logical(numat)
    std::vector<int> iused;             // (max(norbs, numat))
    std::vector<double> storei, storej; // (norbs)
    try {
        latom_loc.resize(molkst_C::numat + 1);
        iused.resize(std::max(molkst_C::norbs, molkst_C::numat) + 1);
        storei.resize(molkst_C::norbs + 1);
        storej.resize(molkst_C::norbs + 1);
    } catch (const std::bad_alloc&) {
        alloc_stat = 1;
    }

    if (alloc_stat != 0) {
        memory_error("addhb");
        return;
    }

    // call hbonds(f, nocc1, nvir1, iused, nij, hblims(nhb))
    // f.data() points at padding[0]; callees index fao(1) at f[1].
    hbonds(common_arrays_C::f.data(), nocc1, nvir1, iused.data(), nij,
           hblims[nhb]);

    if (nij != 0) {
        // call diagg2(..., eigs(nocc1+1:), ...)
        // The Fortran slice starts at Fortran index nocc1+1. With the 1-based
        // vector, the callee's eigv(1) must alias eigs[nocc1+1], so we pass
        // &eigs[0] + nocc1 (one slot before the slice start).
        diagg2(nocc1, nvir1,
               common_arrays_C::eigs.data() + nocc1,
               iused.data(), latom_loc.data(), nij, idiagg,
               storei.data(), storej.data());
    }

    // deallocate(...) runs automatically at scope exit.
    //
    // The Fortran's trailing
    //     if (alloc_stat /= 0) call memory_error("addhb")
    // is a latent bug: deallocate has no stat=, so alloc_stat still holds the
    // (already handled) allocate status and this branch is dead. It is
    // preserved here as a faithful artifact.
}
