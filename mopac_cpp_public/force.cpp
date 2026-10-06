// force.cpp — C++ translation of MOPAC 2016 "force.F90".
// FORCE: force constants, vibrational frequencies, zero-point energy,
// isotope substitution, IRC/DRC driving, thermo interface.
#include "force.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "symmetry_C.h"
#include "parameters_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "to_screen_C.h"
#include "chanel_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace molkst_C {
extern int natoms, ndep, nvar, numat, id, numcal, last;
extern int itemp_1, itemp_2;
extern bool moperr, mozyme, uhf;
extern double gnorm, escf, zpe;
extern std::string keywrd;
extern int iflepo;
}
namespace common_arrays_C {
extern std::vector<double> xparam, grad, errfn, fmatrx, p, q, eigs, eigb;
extern std::vector<int> na, nb, nc, labels, nat, na_store;
extern std::vector<std::vector<int>> loc, lopt;
extern std::vector<std::vector<double>> geo, geoa, coord, c, cb;
}
namespace symmetry_C {
extern std::string name;
extern int igroup;
}
namespace parameters_C {
extern double tore[];
}
namespace funcon_C {
extern double fpc_10, fpc_6, fpc_8, fpc_9, a0, ev;
}
namespace to_screen_C {
extern std::vector<double> dipt, travel, freq, cnorml;
extern std::vector<std::vector<double>> redmas;
}
namespace chanel_C {
extern int iw;
}

extern void upcase(std::string&, int);
extern void mopend(const std::string&);
extern void to_screen(const std::string&);
extern void xyzint(double*, int, int*, int*, int*, double, double*);
extern void gmetry(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
extern double second(int);
extern void compfg(const std::vector<double>&, bool, double&, bool, std::vector<double>&, bool);
extern void axis(double&, double&, double&, double evec[4][4]);
extern double dipole(const std::vector<double>&, std::vector<std::vector<double>>&, std::vector<double>&, int);
extern void symtrz(double*, double*, int, int);
extern void vecprt(double*, int);
extern void write_trajectory(double*, int, double*, double, double, double, double);
extern void frame(std::vector<double>&, int, int);
extern void matou1(double*, double*, int, int&, int, int);
extern void matout(const double*, const double*, int, int&, int);
extern void drc(std::vector<double>&, const std::vector<double>&);
extern double reada(const std::string&, int);
extern void thermo(double, double, double, int, double, double*, int, double);
extern void intfc(double*, double*, double*, int*, int*, int*);
extern void fmat(std::vector<double>&, int&, double, double, std::vector<std::vector<double>>&, double, std::vector<double>&, bool);
extern void freqcy(std::vector<double>&, std::vector<double>&, std::vector<double>&, bool, std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&, bool);
extern void anavib(std::vector<double>&, std::vector<double>&, int, std::vector<std::vector<double>>&, std::vector<double>&, int, std::vector<double>&, std::vector<double>&);
extern void mullik();
extern void chrge(const std::vector<double>&, std::vector<double>&);
extern void rsp(double*, int, double*, double*);
extern void write_path_html();
extern void reverse_aux();

// BLAS DDOT (double precision inner product).
static double ddot(int n, const double* x, int incx, const double* y, int incy) {
    double s = 0.0;
    for (int i = 0, ix = 0, iy = 0; i < n; ++i, ix += incx, iy += incy)
        s += x[ix] * y[iy];
    return s;
}
// Fortran 1-based index(): position (1-based) or 0 when absent.
static int f_index(const std::string& s, const std::string& sub, int start1 = 1) {
    size_t p = s.find(sub, start1 - 1);
    return (p == std::string::npos) ? 0 : (int)p + 1;
}
// Copy a (3,n) row-major 1-based vector<vector<double>> into a flat
// column-major Fortran buffer (xyz[(atom-1)*3 + (x-1)]).
static void to_col3(std::vector<double>& flat, const std::vector<std::vector<double>>& a, int n) {
    flat.assign(3 * n + 1, 0.0);
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= 3; ++j) flat[(i - 1) * 3 + (j - 1)] = a[j][i];
}
static void from_col3(std::vector<std::vector<double>>& a, const std::vector<double>& flat, int n) {
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= 3; ++j) a[j][i] = flat[(i - 1) * 3 + (j - 1)];
}
// xyzint wrapper for (3,n) row-major arrays.
static void xyzint_vec(const std::vector<std::vector<double>>& coord, int numat,
                       std::vector<int>& na, std::vector<int>& nb, std::vector<int>& nc,
                       double degree, std::vector<std::vector<double>>& geo) {
    std::vector<double> xyz, out;
    to_col3(xyz, coord, numat);
    out.assign(3 * numat + 1, 0.0);
    xyzint(&xyz[1], numat, &na[1], &nb[1], &nc[1], degree, &out[1]);
    from_col3(geo, out, numat);
}
// Copy cnorml (packed n*n, column-major flat) into anavib's 2-D (n,n) view.
static void cnorml_to_2d(std::vector<std::vector<double>>& v, const std::vector<double>& cnorml, int n) {
    v.assign(n + 1, std::vector<double>(n + 1, 0.0));
    for (int col = 1; col <= n; ++col)
        for (int row = 1; row <= n; ++row)
            v[row][col] = cnorml[(col - 1) * n + (row - 1)];
}
static void cnorml_from_2d(std::vector<double>& cnorml, const std::vector<std::vector<double>>& v, int n) {
    for (int col = 1; col <= n; ++col)
        for (int row = 1; row <= n; ++row)
            cnorml[(col - 1) * n + (row - 1)] = v[row][col];
}

// Point-group rotational symmetry numbers (force.F90 data statement).
static const int irot[60] = {
    1, 1, 1, 2, 4, 2, 2, 4, 3, 4, 2, 6, 3, 3, 5, 6, 3, 7, 8, 4, 4,
    8, 4, 5, 6, 12, 4, 6, 6, 7, 7, 8, 16, 8, 8, 5, 10, 10, 10, 6, 12, 12,
    14, 14, 14, 8, 16, 12, 12, 24, 12, 24, 24, 24, 1, 2, 1, 1, 1, 1 };

// label98 block: thermo calculation using the vibrational frequencies.
static void force_thermo(std::vector<double>& store, int nvib, double a, double b,
                         double c, bool linear, double escf) {
    using namespace common_arrays_C;
    using namespace molkst_C;
    using namespace to_screen_C;
    using namespace symmetry_C;
    gmetry(geo, coord);
    double sym = irot[symmetry_C::igroup];
    std::fprintf(stdout, "\n SYMMETRY NUMBER FOR POINT-GROUP %s=%3d\n",
                 symmetry_C::name.c_str(), irot[symmetry_C::igroup]);
    int i = f_index(keywrd, " TRANS");
    int j;
    if (i != 0) {
        i = f_index(keywrd, " TRANS=");
        if (i != 0) {
            i = (int)(1 + reada(keywrd, i) + 0.5);
            j = nvib - i + 1;
            std::fprintf(stdout, "\n\n THE LOWEST%3d VIBRATIONS ARE NOT\n  TO BE USED IN THE THERMO CALCULATION\n", i - 1);
        } else {
            std::fprintf(stdout, "\n\n          SYSTEM IS A TRANSITION STATE\n");
            i = 2;
            j = nvib - 1;
        }
    } else {
        std::fprintf(stdout, "\n\n          SYSTEM IS A GROUND STATE\n");
        i = 1;
        j = nvib;
    }
    if (store.size() <= (size_t)nvib) store.assign(nvib + 1, 0.0);
    for (int idx = 1; idx <= nvib; ++idx) store[idx] = freq[idx];
    thermo(a, b, c, linear ? 1 : 0, sym, &store[i], j, escf);
}
void force() {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace funcon_C;
    using namespace to_screen_C;
    auto cleanup = [&]() { dipt.clear(); travel.clear(); freq.clear(); redmas.clear(); };

    double time2 = -1.0e9;
    // Save connectivity and geometry for restarts.
    std::vector<int> nar(natoms + 1, 0), nbr(natoms + 1, 0), ncr(natoms + 1, 0);
    for (int i = 1; i <= natoms; ++i) { nar[i] = na[i]; nbr[i] = nb[i]; ncr[i] = nc[i]; }
    geoa.assign(4, std::vector<double>(natoms + 1, 0.0));
    for (int i = 1; i <= natoms; ++i)
        for (int j = 1; j <= 3; ++j) geoa[j][i] = geo[j][i];
    bool ts = (f_index(keywrd, " FORCETS") + f_index(keywrd, " MINI") != 0);
    std::vector<std::vector<int>> locold(3, std::vector<int>(3 * natoms + 1, 0));
    int nvaold = 0, ndeold = 0;

    if (id == 0 && f_index(keywrd, " NOREOR") == 0 && !ts) {
        // Ensure XYZINT works correctly before the call to DRC.
        int l = 0;
        for (int i = 1; i <= natoms; ++i) {
            if (labels[i] == 99) continue;
            ++l;
            labels[l] = labels[i];
        }
        numat = l;
        natoms = numat;
        if (numat == 1) { cleanup(); return; }  // F90: goto 98 -> single atom, no force output
        xyzint_vec(coord, numat, na, nb, nc, 1.0, geo);
        na_store = na;
        gmetry(geo, coord);
        nvaold = nvar;
        for (int i = 1; i <= nvar; ++i) { locold[1][i] = loc[1][i]; locold[2][i] = loc[2][i]; }
    }
    // Unconditionally convert to Cartesian coordinates.
    if (!ts) nvar = 0;
    ndeold = ndep;
    ndep = 0;
    for (int i = 1; i <= numat + id; ++i) {
        for (int j = 1; j <= 3; ++j) {
            geo[j][i] = coord[j-1][i];
            if (!ts) {
                ++nvar;
                loc[1][nvar] = i;
                loc[2][nvar] = j;
                xparam[nvar] = geo[j][i];
            }
            labels[i] = nat[i];
        }
    }
    if (ts) {
        for (int i = 1; i <= nvar; ++i)
            xparam[i] = geo[loc[2][i]][loc[1][i]];
    }
    if (id > 0) {
        natoms = numat + id;
        nvar -= 3 * id;
        for (int i = numat + 1; i <= numat + id; ++i) labels[i] = 107;
    }
    std::fill(na.begin(), na.end(), 0);
    if (nvar > 8000) {
        std::fprintf(stdout, "    Insufficient memory to run FORCE\n");
        mopend("Insufficient memory to run FORCE");
        return;
    }
    int i = ts ? nvar : 3 * (numat + id);
    int j = 0;

    // Deallocate in case fewer variables were set in the input file.
    grad.clear(); errfn.clear();
    dipt.clear(); travel.clear(); freq.clear(); redmas.clear(); cnorml.clear(); fmatrx.clear();
    cnorml.resize((size_t)i * i + 1, 0.0);
    fmatrx.resize((size_t)i * (i + 1) / 2 + 1, 0.0);
    grad.resize(i + 1, 0.0);
    std::vector<double> ff((size_t)i * i + 1, 0.0);
    errfn.resize(i + 1, 0.0);
    std::vector<double> store((size_t)i * (i + 1) / 2 + 1, 0.0);
    std::vector<double> oldf((size_t)i * (i + 1) / 2 + 1, 0.0);
    dipt.resize(3 * numat + 1, 0.0);
    travel.resize(3 * numat + 1, 0.0);
    freq.resize(3 * numat + 1, 0.0);
    redmas.assign(3 * numat + 1, std::vector<double>(3, 0.0));
    if (i > 0 && cnorml.empty()) {
        std::fprintf(stdout, " Failed to allocate memory in FORCE\n");
        mopend("Failed to allocate memory in FORCE");
        return;
    }
    std::fill(errfn.begin(), errfn.end(), 0.0);

    double tscf = -1.0, tder = -1.0;
    bool prnt = (f_index(keywrd, " RC=") == 0);
    bool debug = (f_index(keywrd, " DFORCE") != 0);
    bool large = (f_index(keywrd, " LARGE") != 0);
    bool restrt = (f_index(keywrd, " RESTART") != 0);
    double time1 = second(1);
    if (!restrt) {
        compfg(xparam, true, escf, true, grad, false);
        if (moperr) { cleanup(); return; }
        std::fprintf(stdout, "\n\n          HEAT OF FORMATION =%17.6f KCALS/MOLE\n", escf);
        time2 = second(1);
        tscf = time2 - time1;
        compfg(xparam, true, escf, false, grad, true);
        double time3 = second(1);
        tder = time3 - time2;
        if (prnt && !ts) {
            std::fprintf(stdout, "\n\n          CARTESIAN COORDINATE DERIVATIVES\n\n   NUMBER  ATOM          X              Y              Z\n");
            int l = 0, iu = 0;
            for (i = 1; i <= natoms; ++i) {
                if (labels[i] == 99) continue;
                ++l;
                int il = iu + 1;
                iu = il + 2;
                if (labels[i] == 107) iu = il;
                std::fprintf(stdout, "%6d%6s%14.6f%13.6f%13.6f\n", l,
                             elemts_C::elemnt[labels[i]].c_str(),
                             grad[il], grad[il + 1], grad[il + 2]);
            }
            std::fprintf(stdout, "\n");
        }
        gnorm = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
        if (!mozyme && f_index(keywrd, " AUX") != 0) {
            chrge(p, q);
            for (i = 1; i <= numat; ++i) {
                int l = nat[i];
                q[i] = parameters_C::tore[l] - q[i];
            }
            double sum = dipole(p, coord, fmatrx, 1);
            (void)sum;
            // symtrz(ca, eigs, 1, .true.): adapt c (n,n) to flat column-major.
            int nsz = (int)c.size() - 1;
            if (nsz > 0) {
                std::vector<double> ca_flat((size_t)nsz * nsz + 1, 0.0);
                for (int cc = 1; cc <= nsz; ++cc)
                    for (int rr = 1; rr <= nsz; ++rr)
                        ca_flat[(cc - 1) * nsz + (rr - 1)] = c[rr][cc];
                symtrz(&ca_flat[1], &eigs[1], 1, 1);
                for (int cc = 1; cc <= nsz; ++cc)
                    for (int rr = 1; rr <= nsz; ++rr)
                        c[rr][cc] = ca_flat[(cc - 1) * nsz + (rr - 1)];
            }
            if (uhf) {
                int nsz2 = (int)cb.size() - 1;
                if (nsz2 > 0) {
                    std::vector<double> cb_flat((size_t)nsz2 * nsz2 + 1, 0.0);
                    for (int cc = 1; cc <= nsz2; ++cc)
                        for (int rr = 1; rr <= nsz2; ++rr)
                            cb_flat[(cc - 1) * nsz2 + (rr - 1)] = cb[rr][cc];
                    symtrz(&cb_flat[1], &eigb[1], 1, 1);
                    for (int cc = 1; cc <= nsz2; ++cc)
                        for (int rr = 1; rr <= nsz2; ++rr)
                            cb[rr][cc] = cb_flat[(cc - 1) * nsz2 + (rr - 1)];
                }
            }
            to_screen("To_file: Normal output");
        }
        std::fprintf(stdout, "\n\n          GRADIENT NORM =%12.5f\n", gnorm);
        if (gnorm > 10.0 && !ts) {
            if (f_index(keywrd, " LET ") != 0) {
                std::fprintf(stdout, "\n\n\n ** GRADIENT IS VERY LARGE, BUT SINCE \"LET\" IS USED, CALCULATION WILL CONTINUE\n");
            } else {
                std::fprintf(stdout, "\n\n\n ** GRADIENT IS TOO LARGE TO ALLOW FORCE MATRIX TO BE CALCULATED, (LIMIT=10) **\n\n");
                std::fprintf(stdout, "\n EITHER ADD 'LET' OR REDUCE GRADIENT\n");
                std::fprintf(stdout, " USING 'TS' OR OTHER GEOMETRY OPTIMIZER\n");
                to_screen("To_file: ERROR: GRADIENT IS TOO LARGE TO ALLOW FORCE MATRIX TO BE CALCULATED");
                mopend("Gradient in FORCE is too large");
                { cleanup(); return; }
            }
        }
    }
    if (f_index(keywrd, " THERMO") != 0 && gnorm > 1.0) {
        std::fprintf(stdout, "\n\n              **** WARNING ****\n\n          GRADIENT IS VERY LARGE FOR A THERMO CALCULATION\n          RESULTS ARE LIKELY TO BE INACCURATE IF THERE ARE\n");
        std::fprintf(stdout, "          ANY LOW-LYING VIBRATIONS (LESS THAN ABOUT 400CM-1)\n");
        std::fprintf(stdout, "          GRADIENT NORM SHOULD BE LESS THAN ABOUT 0.2 FOR THERMO\n          TO GIVE ACCURATE RESULTS\n");
    }
    if (!mozyme && !restrt) mullik();
    if (tscf > 0.0) {
        std::fprintf(stdout, "\n\n          TIME FOR SCF CALCULATION =%10.2f\n", tscf);
        std::fprintf(stdout, "           TIME FOR DERIVATIVES     =%10.2f\n", tder);
    }
    if (ndeold > 0)
        std::fprintf(stdout, "\n\n          SYMMETRY WAS SPECIFIED, BUT CANNOT BE USED HERE\n");
    double a = 0.0, b = 0.0, c = 0.0;
    double rot[4][4] = {};
    if (!ts) axis(a, b, c, rot);
    std::vector<std::vector<double>> store_coord(4, std::vector<double>(numat + 1, 0.0));
    for (i = 1; i <= numat; ++i)
        for (int k = 1; k <= 3; ++k) store_coord[k-1][i] = coord[k-1][i];
    for (i = 1; i <= numat; ++i)
        for (int k = 1; k <= 3; ++k) geo[k][i] = coord[k-1][i];
    if (rot[1][1] > 2.0) { cleanup(); return; }  // dummy use of rot
    int nvib = 3 * numat - 6;
    if (ts) {
        nvib = nvar;
    } else {
        if (std::fabs(c) < 1.0e-20) ++nvib;
        if (id != 0) nvib = 3 * numat - 3;
    }
    if (prnt) {
        std::fprintf(stdout, "\n          ORIENTATION OF MOLECULE IN FORCE CALCULATION\n");
        std::fprintf(stdout, "\n    NO.      ATOM         X         Y         Z\n");
    }
    int l = 0;
    if (!prnt) {
        for (i = 1; i <= natoms; ++i)
            if (labels[i] != 99) ++l;
    } else {
        for (i = 1; i <= natoms; ++i) {
            if (labels[i] == 99) continue;
            ++l;
            std::fprintf(stdout, "%6d%9s%3d%10.4f%10.4f%10.4f\n", l,
                         elemts_C::elemnt[labels[i]].c_str(), labels[i],
                         coord[0][l], coord[1][l], coord[2][l]);
        }
    }
    symtrz(&cnorml[1], &cnorml[1], 2, 0);
    std::vector<std::vector<double>> deldip(4, std::vector<double>(3 * numat + 1, 0.0));
    fmat(fmatrx, nvib, tscf, tder, deldip, escf, ff, ts);
    if (moperr) { cleanup(); return; }
    na[1] = 0;
    for (i = 1; i <= natoms; ++i) { na[i] = nar[i]; nb[i] = nbr[i]; nc[i] = ncr[i]; }
    for (i = 1; i <= natoms; ++i)
        for (int k = 1; k <= 3; ++k) geo[k][i] = geoa[k][i];
    if (nvib < 0) {
        ndep = ndeold;
        nvar = 0;
        iflepo = -1;
        { cleanup(); return; }
    }
    // Print the force matrix as an atom-atom matrix.
    int ij = 0, iu = 0;
    for (i = 1; i <= nvar / 3; ++i) {
        int il = iu + 1;
        iu = il + 2;
        int im1 = i - 1;
        int ju = 0;
        for (int jj2 = 1; jj2 <= im1; ++jj2) {
            int jl = ju + 1;
            ju = jl + 2;
            double sum = 0.0;
            for (int ii = il; ii <= iu; ++ii)
                for (int jj = jl; jj <= ju; ++jj)
                    sum += fmatrx[(ii * (ii - 1)) / 2 + jj] * fmatrx[(ii * (ii - 1)) / 2 + jj];
            ++ij;
            store[ij] = std::sqrt(sum);
        }
        ++ij;
        store[ij] = std::sqrt(fmatrx[((il + 0) * (il + 1)) / 2] * fmatrx[((il + 0) * (il + 1)) / 2] +
                              fmatrx[((il + 1) * (il + 2)) / 2] * fmatrx[((il + 1) * (il + 2)) / 2] +
                              fmatrx[((il + 2) * (il + 3)) / 2] * fmatrx[((il + 2) * (il + 3)) / 2] +
                              2.0 * (fmatrx[((il + 1) * (il + 2)) / 2 - 1] * fmatrx[((il + 1) * (il + 2)) / 2 - 1] +
                                     fmatrx[((il + 2) * (il + 3)) / 2 - 2] * fmatrx[((il + 2) * (il + 3)) / 2 - 2] +
                                     fmatrx[((il + 2) * (il + 3)) / 2 - 1] * fmatrx[((il + 2) * (il + 3)) / 2 - 1]));
    }
    if (debug) {
        std::fprintf(stdout, "\n\n          FULL FORCE MATRIX, INVOKED BY \"DFORCE\"\n");
        if (f_index(keywrd, " NOREOR") == 0) {
            std::fprintf(stdout, "\n Caution: NOREOR is NOT present, therefore system will be oriented\n");
            std::fprintf(stdout, " so that the moments of inertia are along the Cartesian axes.\n");
        }
        i = -nvar;
        vecprt(&fmatrx[1], i);
    }
    if (prnt) {
        double sum = 1.0e-21 * fpc_10 * a0 * a0 / (4.184 * ev * fpc_9);
        std::fprintf(stdout, "\n(To convert to Hartree/Bohr^2, multiply by (10^(-21) x N x a0^2)/(4.184 x 627.51) = %8.4f)\n", sum);
        std::fprintf(stdout, "\n\n          FORCE MATRIX IN MILLIDYNES/ANGSTROM\n");
        i = ts ? nvar / 3 : numat;
        vecprt(&store[1], i);
    }
    int l2 = (nvar * (nvar + 1)) / 2;
    for (i = 1; i <= l2; ++i) store[i] = fmatrx[i];
    if (prnt) axis(a, b, c, rot);
    // A molecule is linear if one or two of a,b,c are much smaller than the largest.
    bool linear = std::fabs((a / (a + b + c)) * (b / (a + b + c)) * (c / (a + b + c))) < 1.0e-10;
    if (prnt)
        std::fprintf(stdout, "\n\n          HEAT OF FORMATION =%17.6f KCALS/MOLE\n", escf);
    for (i = 1; i <= numat; ++i)
        for (int k = 1; k <= 3; ++k) coord[k-1][i] = store_coord[k-1][i];
    if (large) {
        frame(store, numat, 0);
        rsp(&store[1], nvar, &freq[1], &cnorml[1]);
        phase_lock(cnorml, nvar);
        for (i = nvib + 1; i <= nvar; ++i) {
            j = (int)((freq[i] + 50.0) * 0.01);
            freq[i] -= j * 100;
        }
        if (prnt) {
            std::fprintf(stdout, "\n\n          TRIVIAL VIBRATIONS, SHOULD BE ZERO\n");
            std::fprintf(stdout, "\n %9.4f=TX%9.4f=TY%9.4f=TZ%9.4f=RX%9.4f=RY%9.4f=RZ\n",
                         freq[nvib + 1], freq[nvib + 2], freq[nvib + 3],
                         freq[nvib + 4], freq[nvib + 5], freq[nvib + 6]);
            symtrz(&cnorml[1], &freq[1], 2, 1);
            std::fprintf(stdout, "\n\n      MOLECULAR POINT GROUP   :   %s\n", symmetry_C::name.c_str());
            std::fprintf(stdout, "\n\n          EIGENVECTORS\n");
            i = 3 * numat;
            matou1(&cnorml[1], &freq[1], nvib, i, nvib, 5);
            std::fprintf(stdout, "\n\n          FORCE CONSTANTS IN MILLIDYNES/ANGSTROM  (= 10**5 DYNES/CM)\n");
            int nprint = std::min(nvib, 8);
            for (int k = 1; k <= nprint; ++k) std::fprintf(stdout, "%8.5f", freq[k]);
            std::fprintf(stdout, "\n");
            std::fprintf(stdout, "\n\n          ASSOCIATED EIGENVECTORS\n");
            i = -nvar;
            matout(&cnorml[1], &freq[1], nvib, i, nvar);
        }
    }
    itemp_2 = nvar - nvib;
    freqcy(fmatrx, freq, travel, true, deldip, ff, oldf, ts);
    // Calculate zero-point energy.
    // const = 0.5*N*H*C/(1e10*4.184)
    double constv = 0.5 * fpc_10 * fpc_6 * fpc_8 / (1.0e10 * 4.184);
    double sum = 0.0;
    int ni = 0;
    for (i = 1; i <= nvib; ++i) {
        if (freq[i] > 0) sum += freq[i];
        else ++ni;
    }
    zpe = sum * constv;
    if (prnt) {
        std::fprintf(stdout, "\n          ZERO POINT ENERGY%13.3f KCAL/MOL\n", zpe);
        if (ni != 0) {
            std::fprintf(stdout, "\n\n NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT\n");
            std::fprintf(stdout, " ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED\n");
            std::fprintf(stdout, " DOES NOT INCLUDE THE%3d IMAGINARY FREQUENCIES\n", ni);
        }
    }
    std::vector<std::vector<double>> trdip(4, std::vector<double>(std::max(3 * numat, nvar) + 1, 0.0));
    double summ = 0.0;
    for (i = 1; i <= nvar; ++i) {
        double sum1 = 1.0e-20;
        for (j = 1; j <= nvar; ++j)
            sum1 += cnorml[j + (i - 1) * nvar] * cnorml[j + (i - 1) * nvar];
        sum1 = 1.0 / std::sqrt(sum1);
        std::fill(grad.begin(), grad.end(), 0.0);
        for (int k = 1; k <= 3; ++k) {
            sum = 0.0;
            for (j = 1; j <= nvar; ++j)
                sum += cnorml[j + (i - 1) * nvar] * deldip[k][j];
            summ += std::fabs(sum);
            trdip[k][i] = sum * sum1;
        }
        dipt[i] = std::sqrt(trdip[1][i] * trdip[1][i] + trdip[2][i] * trdip[2][i] + trdip[3][i] * trdip[3][i]);
    }
    if (prnt && large) {
        std::fprintf(stdout, "\n\n          FREQUENCIES, REDUCED MASSES AND VIBRATIONAL DIPOLES\n");
        int nto6 = nvar / 6;
        int nrem6 = nvar - nto6 * 6;
        int iinc1 = -5;
        if (nto6 >= 1) {
            for (i = 1; i <= nto6; ++i) {
                std::fprintf(stdout, "\n");
                iinc1 += 6;
                int iinc2 = iinc1 + 5;
                std::fprintf(stdout, "   I");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10d", k);
                std::fprintf(stdout, "\n");
                std::fprintf(stdout, " FREQ(I)");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.4f", freq[k]);
                std::fprintf(stdout, "\n");
                std::fprintf(stdout, " MASS(I)");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", redmas[k][1]);
                std::fprintf(stdout, "\n");
                std::fprintf(stdout, " DIPX(I)");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", trdip[1][k]);
                std::fprintf(stdout, "\n");
                std::fprintf(stdout, " DIPY(I)");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", trdip[2][k]);
                std::fprintf(stdout, "\n");
                std::fprintf(stdout, " DIPZ(I)");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", trdip[3][k]);
                std::fprintf(stdout, "\n");
                std::fprintf(stdout, " DIPT(I)");
                for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", dipt[k]);
                std::fprintf(stdout, "\n");
            }
        }
        if (nrem6 >= 1) {
            std::fprintf(stdout, "\n");
            iinc1 += 6;
            int iinc2 = iinc1 + (nrem6 - 1);
            std::fprintf(stdout, "   I");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10d", k);
            std::fprintf(stdout, "\n");
            std::fprintf(stdout, " FREQ(I)");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.4f", freq[k]);
            std::fprintf(stdout, "\n\n MASS(I)");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", redmas[k][1]);
            std::fprintf(stdout, "\n\n DIPX(I)");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", trdip[1][k]);
            std::fprintf(stdout, "\n DIPY(I)");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", trdip[2][k]);
            std::fprintf(stdout, "\n DIPZ(I)");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", trdip[3][k]);
            std::fprintf(stdout, "\n\n DIPT(I)");
            for (int k = iinc1; k <= iinc2; ++k) std::fprintf(stdout, "%10.5f", dipt[k]);
            std::fprintf(stdout, "\n");
        }
    }
    if (prnt) {
        std::fprintf(stdout, "\n\n          NORMAL COORDINATE ANALYSIS (Total motion = 1 Angstrom)\n");
        j = nvar;
        i = -nvar;
        matou1(&cnorml[1], &freq[1], nvib, j, nvib, 5);
    }
    // Carry out IRC if requested.
    if (f_index(keywrd, " IRC") + f_index(keywrd, " DRC") != 0) {
        if (f_index(keywrd, " HTML") != 0) write_path_html();
        for (i = 1; i <= nvar; ++i) { loc[1][i] = 0; loc[2][i] = 0; }
        nvar = nvaold;
        for (i = 1; i <= nvar; ++i) { loc[1][i] = locold[1][i]; loc[2][i] = locold[2][i]; }
        xyzint_vec(coord, numat, na, nb, nc, 1.0, geo);
        last = 1;
        for (i = 1; i <= numat; ++i)
            for (int k = 1; k <= 3; ++k) store_coord[k-1][i] = coord[k-1][i];
        for (i = 1; i <= numat; ++i)
            for (int k = 1; k <= 3; ++k) geo[k][i] = store_coord[k-1][i];
        std::fill(na.begin(), na.end(), 0);
        itemp_2 = 0;
        std::vector<double> velocity;
        int k = 0;
        if (ts) {
            velocity.assign(3 * numat + 1, 0.0);
            k = (f_index(keywrd, " FORCETS") != 0) ? 1 : 0;
            std::fill(velocity.begin(), velocity.end(), 0.0);
            j = 3;
            for (i = 1; i <= numat; ++i) {
                if (lopt[1][i] == k) {
                    velocity[i * 3 - 2] = cnorml[j - 2];
                    velocity[i * 3 - 1] = cnorml[j - 1];
                    velocity[i * 3 - 0] = cnorml[j - 0];
                    j += 3;
                }
            }
            std::fill(na_store.begin(), na_store.end(), 0);
            drc(velocity, freq);
        } else {
            velocity.assign(3 * numat + 1, 0.0);
            for (i = 1; i <= 3 * numat; ++i) velocity[i] = cnorml[i];
            drc(cnorml, freq);
            for (i = 1; i <= 3 * numat; ++i) cnorml[i] = velocity[i];
        }
        if (moperr) return;
        i = f_index(keywrd, " IRC");
        if (i != 0) {
            // Position of the first keyword after "IRC".
            int sp = f_index(keywrd.substr(i), " ");
            int start1 = i + sp - 1;
            std::string tail = keywrd.substr(start1 - 1);
            if (f_index(tail, "*") != 0) {
                // Double-sided IRC: first reverse the path already written.
                std::vector<double> geo_flat;
                to_col3(geo_flat, geo, numat);
                write_trajectory(&geo_flat[1], 2, nullptr, 0.0, 0.0, 0.0, 0.0);
                reverse_aux();
                last = 1;
                for (i = 1; i <= numat; ++i)
                    for (int kk = 1; kk <= 3; ++kk) geo[kk][i] = store_coord[kk-1][i];
                std::fill(na.begin(), na.end(), 0);
                itemp_1 = 0;
                ++numcal;
                for (i = 1; i <= nvar; ++i) { loc[1][i] = 0; loc[2][i] = 0; }
                nvar = nvaold;
                for (i = 1; i <= nvar; ++i) { loc[1][i] = locold[1][i]; loc[2][i] = locold[2][i]; }
                if (ts) {
                    std::fill(velocity.begin(), velocity.end(), 0.0);
                    j = 3;
                    for (i = 1; i <= numat; ++i) {
                        if (lopt[1][i] == k) {
                            velocity[i * 3 - 2] = -cnorml[j - 2];
                            velocity[i * 3 - 1] = -cnorml[j - 1];
                            velocity[i * 3 - 0] = -cnorml[j - 0];
                            j += 3;
                        }
                    }
                    drc(velocity, freq);
                } else {
                    for (size_t idx = 0; idx < cnorml.size(); ++idx) cnorml[idx] = -cnorml[idx];
                    drc(cnorml, freq);
                }
            }
            ndep = 0;
            for (i = 1; i <= numat; ++i) { lopt[1][i] = 1; lopt[2][i] = 1; }
        } else {
            ndep = ndeold;
            nvar = 0;
            for (i = 1; i <= natoms; ++i)
                for (int kk = 1; kk <= 3; ++kk) geo[kk][i] = geoa[kk][i];
        }
        { cleanup(); return; }
    }
    // Mass-weighted analysis.
    std::vector<double> trav2;
    for (i = 1; i <= 3 * numat; ++i)
        for (int k = 1; k <= 3; ++k) trav2.push_back(deldip[k][i]);
    freqcy(fmatrx, freq, trav2, false, deldip, ff, oldf, ts);
    std::fprintf(stdout, "\n\n          MASS-WEIGHTED COORDINATE ANALYSIS (NORMAL COORDINATES)\n");
    i = -nvar;
    j = nvar;
    matou1(&cnorml[1], &freq[1], nvib, j, nvib, 5);
    if (!ts) {
        std::vector<std::vector<double>> vibs;
        cnorml_to_2d(vibs, cnorml, nvar);
        anavib(freq, dipt, nvar, vibs, store, nvib, fmatrx, ff);
        cnorml_from_2d(cnorml, vibs, nvar);
    }
    if (f_index(keywrd, " THERMO") != 0)
        force_thermo(store, nvib, a, b, c, linear, escf);
    for (i = 1; i <= numat; ++i)
        for (int k = 1; k <= 3; ++k) coord[k-1][i] = store_coord[k-1][i];
    to_screen("To_file: Force output");
    for (i = 1; i <= (nvar * (nvar + 1)) / 2; ++i) fmatrx[i] = oldf[i] * 1.0e-5;
    if (!ts && numat > 1) {
        std::vector<double> geo_flat;
        to_col3(geo_flat, geoa, natoms);
        intfc(&oldf[1], &xparam[1], &geo_flat[1], &nar[1], &nbr[1], &ncr[1]);
    }
    na[1] = 0;
    nvar = 0;
    ndep = ndeold;
    cleanup();
}

// phase_lock: the phase of an eigenvector is set so that the largest
// coefficient is positive. vecs is a packed n*n array in column-major
// Fortran order: vecs[(col-1)*n + (row-1)].
void phase_lock(std::vector<double>& vecs, int n) {
    for (int i = 1; i <= n; ++i) {
        double asum = 0.0;
        double ssum = 0.0;
        for (int j = 1; j <= n; ++j) {
            if (std::fabs(vecs[(i - 1) * n + j]) > asum) {
                asum = std::fabs(vecs[(i - 1) * n + j]);
                ssum = vecs[(i - 1) * n + j];
            }
        }
        if (ssum < 0.0) {
            for (int j = 1; j <= n; ++j)
                vecs[(i - 1) * n + j] = -vecs[(i - 1) * n + j];
        }
    }
}
