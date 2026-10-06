// set_up_MOZYME_arrays.cpp — C++ translation of MOPAC 2016
// "set_up_MOZYME_arrays.F90" (allocates MOZYME working arrays).
#include "set_up_MOZYME_arrays.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "iter_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
using namespace MOZYME_C;

namespace molkst_C {
extern int n2elec, mpack, numat, nelecs, norbs, l123, id;
extern std::string keywrd;
}
namespace MOZYME_C {
}
namespace common_arrays_C {
extern std::vector<double> h, w, wk, p, f, q, dxyz, eigs;
extern std::vector<int> ifact, nfirst, nlast;
extern std::vector<double> errfn, fb;
}
namespace chanel_C {
extern int iw;
}

// Global (non-namespace) MOZYME linkage.
using MOZYME_C::nvirtual;
using MOZYME_C::noccupied;

extern void fillij(bool);
extern void mopend(const std::string&);

static void delete_MOZYME_arrays_internal() {
    using namespace common_arrays_C;
    using namespace iter_C;
    h.clear(); w.clear(); wk.clear(); f.clear(); fb.clear(); p.clear();
    pold.clear(); pold2.clear(); pold3.clear(); pbold.clear(); pbold2.clear();
    ncocc.clear(); ncvir.clear(); nncf.clear(); nnce.clear(); icocc.clear();
    icvir.clear(); kopt.clear(); ncf.clear(); nce.clear(); ifact.clear();
    cocc.clear(); cvir.clear(); q.clear(); dxyz.clear(); errfn.clear();
    p1.clear(); p2.clear(); p3.clear(); MOZYME_C::ws.clear(); fmo.clear();
    ifmo.clear(); MOZYME_C::partf.clear(); MOZYME_C::partp.clear(); MOZYME_C::parth.clear();
    part_dxyz.clear(); nfmo.clear(); eigs.clear(); idiag.clear();
    MOZYME_C::iorbs.clear(); MOZYME_C::nijbo.clear();
    jopt.clear();
}

// Public interface: delete_MOZYME_arrays (set_up_MOZYME_arrays.F90 lines 120-182).
// Fortran also deallocates ions/iopt when numat==0; iopt has no C++ mapping.
void delete_MOZYME_arrays() {
    if (molkst_C::numat == 0) MOZYME_C::ions.clear();
    delete_MOZYME_arrays_internal();
}

void set_up_MOZYME_arrays() {
    using namespace common_arrays_C;
    using namespace iter_C;
    using molkst_C::n2elec; using molkst_C::mpack; using molkst_C::numat;
    using molkst_C::nelecs; using molkst_C::norbs; using molkst_C::l123;
    using molkst_C::id; using molkst_C::keywrd;
    delete_MOZYME_arrays_internal();
    MOZYME_C::iorbs.resize(numat + 1);
    for (int i = 1; i <= numat; ++i) MOZYME_C::iorbs[i] = nlast[i] - nfirst[i] + 1;
    fillij(true);
    noccupied = nelecs / 2;
    nvirtual = norbs - noccupied;
    int ipad2 = std::min(numat + 10, (int)std::lround(std::pow((double)numat, 0.25) * 50));
    int ipad4 = (int)std::lround((double)ipad2 * norbs / numat);
    icocc_dim = ipad2 * noccupied + 1000;
    icvir_dim = ipad2 * nvirtual + 1000;
    cocc_dim = ipad4 * noccupied + 1000;
    cvir_dim = ipad4 * nvirtual + 1000;
    h.resize(mpack + 1); w.resize(n2elec + 1);
    if (id > 0) wk.resize(n2elec + 1);
    pold.resize(mpack + 1); pold2.resize(2); pold3.resize(2);
    pbold.resize(2); pbold2.resize(2);
    ncocc.resize(noccupied + 1); ncvir.resize(nvirtual + 1);
    nncf.resize(noccupied + 1); nnce.resize(nvirtual + 1);
    icocc.resize(icocc_dim + 1); icvir.resize(icvir_dim + 1);
    jopt.resize(numat + 1); kopt.resize(numat + 1);
    ncf.resize(noccupied + 1); nce.resize(nvirtual + 1);
    p.resize(mpack + 1); f.resize(mpack + 1); q.resize(numat + 1);
    ifact.resize(norbs + 1);
    cocc.resize(cocc_dim + 1); cvir.resize(cvir_dim + 1);
    if (keywrd.find(" RESTART") != std::string::npos ||
        keywrd.find(" 1SCF") == std::string::npos ||
        keywrd.find(" GRAD") != std::string::npos) {
        if (MOZYME_C::rapid) part_dxyz.resize(3 * numat * l123 + 1);
        dxyz.resize(3 * numat * l123 + 1);
        errfn.resize(3 * numat * l123 + 1);
    }
    p1.resize(norbs + 1); p2.resize(norbs + 1); p3.resize(norbs + 1);
    nfmo.resize(norbs + 1); idiag.resize(norbs + 1);
    MOZYME_C::ws.resize(norbs + 1); eigs.resize(norbs + 1);
    if (MOZYME_C::rapid) {
        MOZYME_C::partf.resize(mpack + 1); MOZYME_C::partp.resize(mpack + 1); MOZYME_C::parth.resize(mpack + 1);
    } else {
        MOZYME_C::partf.resize(2); MOZYME_C::partp.resize(2); MOZYME_C::parth.resize(2);
        mode = 0;
    }
    fmo_dim = std::min(300, norbs) * norbs;
    fmo.resize(fmo_dim + 1);
    ifmo.resize(2, std::vector<int>(fmo_dim + 1, 0));
    std::fill(icocc.begin(), icocc.end(), 0);
    std::fill(icvir.begin(), icvir.end(), 0);
    fillij(false);
    ifact[1] = 0;
    for (int i = 1; i <= norbs - 1; ++i) ifact[i + 1] = i + ifact[i];
    numred = numat;
    norred = norbs;
    nelred = nelecs;
}
