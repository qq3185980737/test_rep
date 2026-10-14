// test_batchB1.cpp — batch tests: den_in_out/forsav/parsav/mult_symm_AB/
// set_up_MOZYME_arrays/local_for_MOZYME/diag_for_GPU/thermo.
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include "iter_C.h"

namespace molkst_C {
extern int numat, norbs, nelecs, n2elec, mpack, l123, id, natoms, nvar, ndep;
extern std::string keywrd, line, title, koment;
extern double mol_weight;
extern double cutofp;
extern bool uhf, mozyme;
extern int ilim;
}
namespace chanel_C {
extern std::string restart_fn, density_fn;
}
namespace common_arrays_C {
extern std::vector<double> p, pa, pb, xparam, atmass, q;
extern std::vector<int> nfirst, nlast;
extern std::vector<std::vector<double>> coord;
extern std::vector<std::string> simbol, txtatm;
}
namespace MOZYME_C {
extern std::vector<int> iorbs;
extern bool rapid;
}
namespace iter_C {
extern std::vector<double> pold;
}
extern std::vector<int> nce, ncf, ncocc, ncvir, nnce, nncf, icocc, icvir;
extern std::vector<double> cocc, cvir;
extern int cocc_dim, cvir_dim, icocc_dim, icvir_dim;
extern std::vector<std::vector<int>> nijbo;

// Stubs for routines not under test (real impls for the rest are linked).
void geout(int) {}
void memory_error(const std::string&) {}
void web_message(const char*) {}
void xxx(char, int, int, int, int, std::string&) {}

// Fortran-name wrappers required by fillij / xyzint.
extern double reada(const std::string&, int);
#include "bangle.h"
#include "dihed.h"
extern "C" double reada_(const char* s, int* i, int len) {
    return reada(std::string(s, (size_t)len), *i);
}
extern "C" void bangle_(double* xyz, int* a, int* b, int* c, double* sum) {
    std::vector<std::vector<double>> v(4, std::vector<double>(4, 0.0));
    for (int d = 1; d <= 3; ++d) { v[d][*a] = xyz[(*a - 1) * 3 + d - 1]; v[d][*b] = xyz[(*b - 1) * 3 + d - 1]; v[d][*c] = xyz[(*c - 1) * 3 + d - 1]; }
    bangle(v, *a, *b, *c, *sum);
}
extern "C" void dihed_(double* xyz, int* a, int* b, int* c, int* d, double* dih) {
    std::vector<std::vector<double>> v(4, std::vector<double>(5, 0.0));
    for (int k = 0; k < 3 * 5; ++k) v[k % 3 + 1][k / 3 + 1] = xyz[k];
    dihed(v, *a, *b, *c, *d, *dih);
}

#include "den_in_out.h"
#include "forsav.h"
#include "parsav.h"
#include "mult_symm_AB.h"
#include "set_up_MOZYME_arrays.h"
#include "local_for_MOZYME.h"
#include "diag_for_GPU.h"
#include "thermo.h"

int main() {
    bool ok = true;
    // ---------------- den_in_out (RHF round-trip) ----------------
    molkst_C::mozyme = false; molkst_C::uhf = false;
    molkst_C::norbs = 2; molkst_C::numat = 1;
    molkst_C::keywrd = "";
    chanel_C::density_fn = "tB1.density";
    common_arrays_C::pa.assign(3, 0.0);
    common_arrays_C::pb.assign(3, 0.0);
    common_arrays_C::p.assign(3, 0.0);
    common_arrays_C::pa[1] = 0.25; common_arrays_C::pa[2] = 0.5;
    den_in_out(1);
    common_arrays_C::pa[1] = common_arrays_C::pa[2] = -1.0;
    common_arrays_C::p[1] = common_arrays_C::p[2] = -1.0;
    den_in_out(0);
    if (std::fabs(common_arrays_C::p[1] - 0.5) > 1e-9 || std::fabs(common_arrays_C::p[2] - 1.0) > 1e-9) {
        std::printf("FAIL den_in_out p=%g %g\n", common_arrays_C::p[1], common_arrays_C::p[2]);
        ok = false;
    }

    // ---------------- forsav round-trip ----------------
    chanel_C::restart_fn = "tB1.restart";
    molkst_C::norbs = 4; molkst_C::numat = 2;
    int nvar = 3, ipt = 2, jstart = 5;
    double time = 2e7;          // >1e7, must be folded
    double refh = -12.34;
    std::vector<std::vector<double>> deldip(4, std::vector<double>(3, 0.0));
    std::vector<double> fmatrx(7, 0.0), coordv(4, 0.0), evecs(10, 0.0), fconst(4, 0.0);
    for (int i = 1; i <= 3; ++i) { coordv[i] = i * 0.5; fconst[i] = i; }
    fmatrx[1] = 1.1; fmatrx[2] = 2.2; fmatrx[3] = 3.3;
    deldip[1][1] = 0.1; deldip[2][1] = 0.2; deldip[3][1] = 0.3;
    deldip[1][2] = 0.4; deldip[2][2] = 0.5; deldip[3][2] = 0.6;
    evecs[1] = 9.9;
    double time_saved = time;
    forsav(time, deldip, ipt, fmatrx, coordv, nvar, refh, evecs, jstart, fconst);
    // Re-read with fresh buffers.
    double time2 = 0, refh2 = 0; int ipt2 = 0, jstart2 = 0;
    std::vector<double> fmatrx2(7, -1.0), coordv2(4, -1.0), evecs2(10, -1.0), fconst2(4, -1.0);
    std::vector<std::vector<double>> deldip2(4, std::vector<double>(3, -1.0));
    int nvar2 = nvar;
    double time_in = 0;
    forsav(time_in, deldip2, ipt2, fmatrx2, coordv2, nvar2, refh2, evecs2, jstart2, fconst2);
    if (std::fabs(time2) + std::fabs(ipt2 - 2) + std::fabs(refh2 + 12.34) > 1e-9 ||
        std::fabs(fmatrx2[1] - 1.1) > 1e-9 || std::fabs(coordv2[2] - 1.0) > 1e-9 ||
        std::fabs(deldip2[3][2] - 0.6) > 1e-9 || std::fabs(evecs2[1] - 9.9) > 1e-9 ||
        jstart2 != 5 || std::fabs(fconst2[3] - 3.0) > 1e-9) {
        std::printf("FAIL forsav round-trip\n"); ok = false;
    }
    if (time_in != 0) time2 = time_in;  // keep read value visible

    // ---------------- parsav round-trip ----------------
    molkst_C::norbs = 4; molkst_C::numat = 2;
    chanel_C::restart_fn = "tB1.parsav";
    int n = 3, m = 2;
    double q[4], r[25], efslst[4], xlast[4];
    int iiium[6] = {1, 2, 3, 4, 5, 6};
    for (int i = 0; i < 4; ++i) q[i] = 0.0;
    for (int i = 0; i < 25; ++i) r[i] = 0.0;
    for (int i = 1; i <= 3; ++i) { efslst[i] = i; xlast[i] = 3 * i; }
    q[0] = 1.5; q[2] = 2.5; q[1] = 3.5; q[3] = 4.5;  // m=2: q(1,1),q(2,1),q(1,2),q(2,2)
    r[0] = 7.5; r[6] = 8.5;
    molkst_C::keywrd = "";
    parsav(1, n, m, q, r, efslst, xlast, iiium);
    double q2[4], r2[25], efslst2[4], xlast2[4];
    int iiium2[6] = {0};
    for (int i = 0; i < 4; ++i) { q2[i] = -1; efslst2[i] = -1; xlast2[i] = -1; }
    for (int i = 0; i < 25; ++i) r2[i] = -1;
    int m2 = 0, n2 = 3;  // n must carry the expected size on entry (Fortran inout)
    parsav(0, n2, m2, q2, r2, efslst2, xlast2, iiium2);
    if (std::fabs(q2[0] - 1.5) > 1e-9 || std::fabs(r2[0] - 7.5) > 1e-9 ||
        std::fabs(efslst2[2] - 2.0) > 1e-9 || std::fabs(xlast2[1] - 3.0) > 1e-9 ||
        m2 != 2 || iiium2[5] != 6) {
        std::printf("FAIL parsav round-trip q=%g r=%g m=%d\n", q2[0], r2[0], m2);
        ok = false;
    }

    // ---------------- mult_symm_AB (iopc=3) ----------------
    {
        double a[3] = {2, 1, 3};   // [[2,1],[1,3]] upper-packed
        double b[3] = {4, 2, 5};   // [[4,2],[2,5]]
        double c[3] = {0, 0, 0};
        mult_symm_AB(a, b, 1.0, 2, 3, c, 0.0, 3);
        if (std::fabs(c[0] - 10.0) > 1e-9 || std::fabs(c[1] - 9.0) > 1e-9 || std::fabs(c[2] - 17.0) > 1e-9) {
            std::printf("FAIL mult_symm_AB c=%g %g %g\n", c[0], c[1], c[2]);
            ok = false;
        }
    }

    // ---------------- set_up_MOZYME_arrays ----------------
    molkst_C::numat = 2; molkst_C::norbs = 4; molkst_C::nelecs = 2;
    molkst_C::l123 = 1; molkst_C::id = 0;
    molkst_C::cutofp = 0.0;
    molkst_C::keywrd = "";
    MOZYME_C::rapid = false;
    common_arrays_C::nfirst.assign(3, 0); common_arrays_C::nlast.assign(3, 0);
    common_arrays_C::nfirst[1] = 1; common_arrays_C::nlast[1] = 2;
    common_arrays_C::nfirst[2] = 3; common_arrays_C::nlast[2] = 4;
    common_arrays_C::coord.assign(3, std::vector<double>(3, 0.0));
    common_arrays_C::coord[0][1] = 0; common_arrays_C::coord[1][1] = 0; common_arrays_C::coord[2][1] = 0;
    common_arrays_C::coord[0][2] = 1; common_arrays_C::coord[1][2] = 0; common_arrays_C::coord[2][2] = 0;
    set_up_MOZYME_arrays();
    if (MOZYME_C::iorbs[1] != 2 || MOZYME_C::iorbs[2] != 2) {
        std::printf("FAIL set_up iorbs\n"); ok = false;
    }
    if (molkst_C::n2elec < 2000 || molkst_C::mpack < 10) {
        std::printf("FAIL set_up sizes n2elec=%d mpack=%d\n", molkst_C::n2elec, molkst_C::mpack);
        ok = false;
    }
    if (fmo_dim < 16 || cocc_dim < 1000) {
        std::printf("FAIL set_up dims fmo_dim=%d cocc_dim=%d\n", fmo_dim, cocc_dim);
        ok = false;
    }

    // ---------------- local_for_MOZYME (OCCUPIED; single LMO -> no rotation) --
    noccupied = 1; nvirtual = 2;
    icocc_dim = 100; cocc_dim = 100;
    ncf.assign(2, 0); ncf[1] = 1;            // 1 atom in the LMO
    nncf.assign(2, 0); nncf[1] = 0;
    ncocc.assign(2, 0); ncocc[1] = 0;
    icocc.assign(101, 0); icocc[1] = 1;      // atom 1
    cocc.assign(101, 0.0); cocc[1] = 0.5; cocc[2] = 0.5;
    local_for_MOZYME("OCCUPIED");
    // Coefficients must be unchanged (single LMO).
    if (std::fabs(cocc[1] - 0.5) > 1e-12) { std::printf("FAIL local_for_MOZYME\n"); ok = false; }

    // ---------------- diag_for_GPU (2x2 analytic) ----------------
    {
        double fao[3] = {1.0, 0.5, 2.0};   // packed lower: F11,F21,F22
        double vec[4] = {1, 0, 0, 1};      // identity (column-major 2x2)
        double eig[2] = {1.0, 2.0};
        diag_for_GPU(fao, vec, 1, eig, 2, 3);
        double n1 = vec[0] * vec[0] + vec[1] * vec[1];
        double n2 = vec[2] * vec[2] + vec[3] * vec[3];
        double dot = vec[0] * vec[2] + vec[1] * vec[3];
        if (std::fabs(n1 - 1.0) > 1e-9 || std::fabs(n2 - 1.0) > 1e-9 || std::fabs(dot) > 1e-9) {
            std::printf("FAIL diag_for_GPU n1=%g n2=%g dot=%g\n", n1, n2, dot);
            ok = false;
        }
    }

    // ---------------- thermo (atom: nvibs=0) ----------------
    molkst_C::title = "TEST"; molkst_C::koment = "THERMO TEST";
    molkst_C::mol_weight = 1.0079;
    common_arrays_C::T_range.assign(301, 0.0);
    common_arrays_C::HOF_tot.assign(301, 0.0);
    common_arrays_C::H_tot.assign(301, 0.0);
    common_arrays_C::Cp_tot.assign(301, 0.0);
    common_arrays_C::S_tot.assign(301, 0.0);
    double vibs1[1] = {0.0};
    thermo(1.0, 1.0, 1.0, 0, 1.0, vibs1, 0, 0.0);
    if (molkst_C::ilim <= 1) { std::printf("FAIL thermo ilim=%d\n", molkst_C::ilim); ok = false; }
    if (std::fabs(common_arrays_C::Cp_tot[1] - 4.9681) > 1e-3) {  // 2.5R
        std::printf("FAIL thermo Cp=%g\n", common_arrays_C::Cp_tot[1]);
        ok = false;
    }

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
