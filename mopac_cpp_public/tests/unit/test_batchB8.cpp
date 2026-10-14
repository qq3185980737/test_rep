// test_batchB8.cpp — tests: force() control flow (water, Cartesian input)
// and phase_lock (independent, fully verified).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "force.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"

namespace molkst_C {
extern int natoms, ndep, nvar, numat, id, numcal, last;
extern int itemp_1, itemp_2;
extern bool moperr, mozyme, uhf;
extern double gnorm, escf, zpe;
extern std::string keywrd;
extern int iflepo;
}
namespace common_arrays_C {
extern std::vector<double> xparam, grad, errfn, fmatrx, p, q;
extern std::vector<int> na, nb, nc, labels, nat;
extern std::vector<std::vector<int>> loc;
extern std::vector<std::vector<double>> geo, geoa, coord;
}
namespace chanel_C { extern int iw; }

// ---- heavy-dependency stubs (real implementations not yet ported) ----
void compfg(const std::vector<double>&, bool, double& escf, bool, std::vector<double>& grad, bool lgrad) {
    escf = -57.8;
    if (lgrad)
        for (size_t i = 1; i < grad.size(); ++i) grad[i] = 0.1;
}
double second(int) { static double t = 0.0; t += 0.01; return t; }
void axis(double& a, double& b, double& c, double evec[4][4]) {
    a = b = c = 1.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) evec[i][j] = 0.0;
    evec[1][1] = 1.0;
}
void fmat(std::vector<double>& fmatrx, int& nreal, double, double,
          std::vector<std::vector<double>>&, double, std::vector<double>&, bool) {
    for (int k = 1; k <= nreal; ++k) fmatrx[(k * (k - 1)) / 2 + k] = 1.0;
}
void freqcy(std::vector<double>& fmatrx, std::vector<double>& freq, std::vector<double>&,
            bool, std::vector<std::vector<double>>&, std::vector<double>&,
            std::vector<double>& oldf, bool) {
    for (size_t i = 1; i < freq.size(); ++i) freq[i] = 400.0 + 100.0 * i;
    for (size_t i = 0; i < oldf.size(); ++i) oldf[i] = fmatrx[i];
}
void anavib(std::vector<double>&, std::vector<double>&, int,
            std::vector<std::vector<double>>&, std::vector<double>&, int,
            std::vector<double>&, std::vector<double>&) {}
void symtrz(double*, double*, int, int) {}
void vecprt(double*, int) {}
void write_trajectory(double*, int, double*, double, double, double, double) {}
void frame(std::vector<double>&, int, int) {}
void matou1(double*, double*, int, int&, int, int) {}
void matout(const double*, const double*, int, int&, int) {}
void drc(std::vector<double>&, const std::vector<double>&) {}
double reada(const std::string&, int) { return 0.0; }
void thermo(double, double, double, int, double, double*, int, double) {}
void intfc(double*, double*, double*, int*, int*, int*) {}
void mullik() {}
void chrge(const std::vector<double>&, std::vector<double>&) {}
double dipole(const std::vector<double>&, std::vector<std::vector<double>>&, std::vector<double>&, int) { return 0.0; }
void rsp(double*, int, double*, double*) {}
void write_path_html() {}
void reverse_aux() {}
void to_screen(const std::string&) {}
void mopend(const std::string&) {}
void upcase(std::string&, int) {}

int main() {
    bool ok = true;
    // ---- phase_lock test (packed n*n column-major) ----
    {
        std::vector<double> v = { 0.0, 1.0, 2.0, -3.0, 1.0 };  // col1=[1,2] col2=[-3,1]
        phase_lock(v, 2);
        if (v[1] != 1.0 || v[2] != 2.0) { std::fprintf(stderr, "FAIL P1 col1 %f %f\n", v[1], v[2]); ok = false; }
        if (v[3] != 3.0 || v[4] != -1.0) { std::fprintf(stderr, "FAIL P1 col2 %f %f\n", v[3], v[4]); ok = false; }
        // col with largest coeff negative -> flipped; largest positive -> kept.
        std::vector<double> w = { 0.0, 1.0, -2.0, -3.0, 1.0 };  // col1=[1,-2] col2=[-3,1]
        phase_lock(w, 2);
        if (w[1] != -1.0 || w[2] != 2.0) { std::fprintf(stderr, "FAIL P2 col1 %f %f\n", w[1], w[2]); ok = false; }
        if (w[3] != 3.0 || w[4] != -1.0) { std::fprintf(stderr, "FAIL P2 col2 %f %f\n", w[3], w[4]); ok = false; }
        // all-zero column stays zero.
        std::vector<double> z = { 0.0, 0.0, 0.0, 0.0, 0.0 };
        phase_lock(z, 2);
        for (int i = 1; i <= 4; ++i)
            if (z[i] != 0.0) { std::fprintf(stderr, "FAIL P3\n"); ok = false; }
    }
    // ---- force() control flow ----
    {
        molkst_C::natoms = 3; molkst_C::numat = 3;
        molkst_C::id = 0; molkst_C::numcal = 1; molkst_C::last = 0;
        molkst_C::nvar = 0; molkst_C::ndep = 0; molkst_C::iflepo = 0;
        molkst_C::moperr = false; molkst_C::mozyme = false; molkst_C::uhf = false;
        molkst_C::keywrd = "";
        common_arrays_C::na.assign(8, 0);
        common_arrays_C::nb.assign(8, 0);
        common_arrays_C::nc.assign(8, 0);
        common_arrays_C::nat = { 0, 8, 1, 1 };
        common_arrays_C::labels = { 0, 8, 1, 1 };
        common_arrays_C::coord.assign(4, std::vector<double>(4, 0.0));
        // Water: O at origin, H1 along +x, H2 in the xy-plane, H-O-H ~104.5 deg.
        common_arrays_C::coord[1][2] = 0.9572;
        common_arrays_C::coord[1][3] = -0.239987;
        common_arrays_C::coord[2][3] = 0.926627;
        common_arrays_C::geo.assign(4, std::vector<double>(4, 0.0));
        common_arrays_C::loc.assign(3, std::vector<int>(10, 0));
        common_arrays_C::xparam.assign(32, 0.0);
        common_arrays_C::grad.assign(32, 0.0);
        force();
        // Fortran force() ends with nvar=0, ndep=ndeold; verify observable results.
        if (molkst_C::itemp_2 != 6) { std::fprintf(stderr, "FAIL F1 n_trivial=%d\n", molkst_C::itemp_2); ok = false; }
        if (std::fabs(molkst_C::gnorm - 0.3) > 1e-9) { std::fprintf(stderr, "FAIL F2 gnorm=%f\n", molkst_C::gnorm); ok = false; }
        // freq stub fills 500..1300 cm^-1; nvib=3 -> zpe = 1800*const.
        double constv = 0.5 * 6.0221367e23 * 6.6260755e-27 * 2.99792458e10 / (1.0e10 * 4.184);
        double expect = 1800.0 * constv;
        if (std::fabs(molkst_C::zpe - expect) > 1e-6) { std::fprintf(stderr, "FAIL F3 zpe=%f expect=%f\n", molkst_C::zpe, expect); ok = false; }
        if (molkst_C::itemp_2 != 6) { std::fprintf(stderr, "FAIL F4 n_trivial=%d\n", molkst_C::itemp_2); ok = false; }
        if (molkst_C::ndep != 0) { std::fprintf(stderr, "FAIL F5 ndep=%d\n", molkst_C::ndep); ok = false; }
    }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
