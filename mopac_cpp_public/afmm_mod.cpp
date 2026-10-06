// afmm_mod.cpp — C++ translation of MOPAC 2016 "afmm_mod.F90" (~1389 lines).
//
// Adaptive fast multipole method.  Fortran 1-based indexing preserved for
// index/level arrays; multipole arrays use the RArr/CArr offset wrappers.

#include "afmm_mod.h"

#include <algorithm>
#include <cmath>
#include <complex>

#include "chanel_C.h"  // iw (debug output only; unused when DEBUG is false)

namespace afmm_mod {
namespace {

const bool DEBUG = false;
const bool PRINT = false;
const double small = 1.0e-3;
const double PI = 3.14159265358979323846;

// Active (currently-selected) tesselation arrays.
std::vector<Box> b;
std::vector<double> d;
std::vector<int> level, ind, ind1;

int nbox = 0, nlev = 0;
std::vector<int> points_index;

RArr pmn, y_norm, amn;
double fact[2 * P + 1] = {};

int max_tess = 0, current_tess = 0;
std::vector<Tess> tess;

// Internal helpers (not part of the public list).
void store_tesselation(int handle, int n);
void split_box(int parent, int level, double x0, double y0, double z0, double d_,
               const std::vector<std::vector<double>>& coord, int nc, int n,
               int b_pos, int ind_pos, int& m);
void split_index(int k, const std::vector<std::vector<double>>& coord, int nc,
                 double g, int ind_pos, int n, int& n1, int& n2);
void normal_multipoles(double x0, double y0, double z0, const int* ind, int n,
                       RArr& pmn, const RArr& y_norm,
                       const std::vector<std::vector<double>>& c, int nc,
                       const std::vector<double>& q, CArr& psi, int p);
void afmm_step1(const std::vector<std::vector<double>>& coord, int nc, int n,
                const std::vector<double>& q, const TrickyFn& tricky);
void afmm_step2();
void afmm_step35();

// (0,1)^m for integer m (Fortran (Cmplx(0,1))**m).
std::complex<double> imag_pow(int m) {
    switch (m & 3) {
        case 0: return {1, 0};
        case 1: return {0, 1};
        case 2: return {-1, 0};
        default: return {0, -1};
    }
}

}  // namespace

// ===========================================================================
void prepare_tesselations(int n1, int n2, std::vector<int>& out_ind, int n,
                          int& ierr) {
    if (!tess.empty()) {
        for (int j = 1; j <= max_tess; ++j) {
            if (!tess[j].b.empty()) {
                tess[j].b.clear();
                tess[j].d.clear();
                tess[j].level.clear();
                tess[j].ind.clear();
                tess[j].ind1.clear();
            }
        }
        tess.clear();
        points_index.clear();
    }
    tess.resize(n + 1);  // 1-based
    points_index.assign(std::max(n1, n2) + 1, 0);
    for (int j = 1; j <= n; ++j) {
        tess[j].b.clear();
        tess[j].d.clear();
        tess[j].level.clear();
        tess[j].ind.clear();
        tess[j].ind1.clear();
    }
    out_ind.assign(n + 1, 0);
    for (int i = 1; i <= n; ++i) out_ind[i] = i;
    max_tess = n;
    current_tess = 0;
    ierr = 0;
}

// ===========================================================================
void set_tesselation(int handle, int& ierr) {
    if (handle > 0 && handle <= max_tess) {
        ierr = 0;
        if (current_tess != handle) {
            current_tess = handle;
            b = tess[handle].b;
            d = tess[handle].d;
            level = tess[handle].level;
            ind = tess[handle].ind;
            ind1 = tess[handle].ind1;
            nbox = tess[handle].nbox;
            nlev = tess[handle].nlev;
        }
        for (int i = 1; i <= tess[handle].n; ++i) points_index[i] = ind[i];
    } else {
        ierr = -1;
    }
}

namespace {
void store_tesselation(int handle, int n) {
    current_tess = handle;
    tess[handle].b = b;
    tess[handle].d = d;
    tess[handle].level = level;
    tess[handle].ind = ind;
    tess[handle].ind1 = ind1;
    tess[handle].nbox = nbox;
    tess[handle].nlev = nlev;
    tess[handle].n = n;
}
}  // namespace

// ===========================================================================
void divide_box(const std::vector<std::vector<double>>& coord, int n1, int n,
                double min_d, int min_p, int max_tes, int handle, int& ierr) {
    double xmax, xmin, ymax, ymin, zmax, zmin;
    double dx, dy, dz, x0, y0, z0, t;
    int max_boxs, max_lev;
    int i, i0, i1, j, j0, k, m;
    bool division;

    if (n1 < 3) { ierr = -10; return; }

    xmax = xmin = coord[0][1];
    ymax = ymin = coord[1][1];
    zmax = zmin = coord[2][1];
    for (i = 2; i <= n; ++i) {
        if (xmax < coord[0][i]) xmax = coord[0][i];
        else if (xmin > coord[0][i]) xmin = coord[0][i];
        if (ymax < coord[1][i]) ymax = coord[1][i];
        else if (ymin > coord[1][i]) ymin = coord[1][i];
        if (zmax < coord[2][i]) zmax = coord[2][i];
        else if (zmin > coord[2][i]) zmin = coord[2][i];
    }
    dx = xmax - xmin; dy = ymax - ymin; dz = zmax - zmin;
    t = dx;
    if (t < dy) t = dy;
    if (t < dz) t = dz;
    x0 = (xmax + xmin) / 2.0;
    y0 = (ymax + ymin) / 2.0;
    z0 = (zmax + zmin) / 2.0;
    t = t / 2.0;
    xmin = x0 - t; xmax = x0 + t;
    ymin = y0 - t; ymax = y0 + t;
    zmin = z0 - t; zmax = z0 + t;

    if (max_tes < 3) {
        max_lev = (int)std::min(std::log(t / min_d) / std::log(2.0) + 1.0,
                                std::log((double)n / (double)min_p) / std::log(8.0) + 2.0);
    } else {
        max_lev = max_tes;
    }
    if (max_lev < 3) max_lev = 3;

    max_boxs = 0;
    j = 1;
    for (i = 1; i <= max_lev; ++i) {
        max_boxs = max_boxs + j * 8;
        j = j * 8;
    }

    b.assign(max_boxs + 1, Box{});  // 1-based
    d.assign(max_lev + 1, 0.0);
    level.assign(max_lev + 2, 0);
    ind.assign(n + 2, 0);
    ind1.assign(std::max(n + 1, 189 * 6) + 1, 0);

    for (i = 1; i <= n; ++i) ind[i] = i;

    split_box(0, 1, x0, y0, z0, t * 0.5, coord, n1, n, 0, 0, m);
    d[1] = t;
    level[1] = 1;
    level[2] = m + 1;
    nbox = m;

    for (i = 1; i <= max_lev - 1; ++i) {
        division = false;
        t = t / 2.0;
        for (j = level[i]; j <= level[i + 1] - 1; ++j) {
            if (b[j].n_p > min_p) {
                split_box(j, i + 1, b[j].x0, b[j].y0, b[j].z0, t * 0.5,
                          coord, n1, b[j].n_p, nbox, b[j].p_start - 1, m);
                division = true;
                b[j].n_childs = m;
                for (k = 1; k <= m; ++k) b[j].childs[k] = nbox + k;
                nbox = nbox + m;
            }
        }
        if (division) {
            nlev = i + 1;
            d[i + 1] = t;
            level[i + 2] = nbox + 1;
            for (j = level[i + 1]; j <= level[i + 2] - 1; ++j)
                b[j].p_start = b[j].p_start + b[b[j].parent].p_start - 1;
        } else {
            nlev = i;
            break;
        }
    }

    // Build L1 list.
    for (i = 1; i <= nlev; ++i) {
        t = d[i] + small;
        for (j = level[i]; j <= level[i + 1] - 1; ++j) {
            for (k = level[i]; k <= level[i + 1] - 1; ++k) {
                if (j == k) continue;
                dx = std::fabs(b[j].x0 - b[k].x0);
                dy = std::fabs(b[j].y0 - b[k].y0);
                dz = std::fabs(b[j].z0 - b[k].z0);
                if (dx < t && dy < t && dz < t) {
                    if (b[j].n_childs != 0) {
                        b[j].n_l1 = b[j].n_l1 + 1;
                        b[j].l1[b[j].n_l1] = k;
                    } else if (j > k || b[k].n_childs != 0) {
                        b[j].n_l1 = b[j].n_l1 + 1;
                        b[j].l1[b[j].n_l1] = k;
                    }
                }
            }
        }
    }
    // Build L2 list.
    for (i = level[2]; i <= nbox; ++i) {
        t = d[b[i].level] + small;
        i0 = b[i].parent;
        for (j = 1; j <= b[i0].n_l1; ++j) {
            j0 = b[i0].l1[j];
            for (k = 1; k <= b[j0].n_childs; ++k) {
                i1 = b[j0].childs[k];
                dx = std::fabs(b[i].x0 - b[i1].x0);
                dy = std::fabs(b[i].y0 - b[i1].y0);
                dz = std::fabs(b[i].z0 - b[i1].z0);
                if (!(dx < t && dy < t && dz < t)) {
                    b[i].n_l2 = b[i].n_l2 + 1;
                    b[i].l2[b[i].n_l2] = i1;
                }
            }
        }
    }

    ierr = 0;
    store_tesselation(handle, n);
}

// ===========================================================================
int count_short_ints(const std::vector<std::vector<double>>& coord, int nc,
                     bool compute, const SimulateDirIntFn& simulate_dir_int) {
    int n = 0;
    for (int i = 1; i <= nbox; ++i) {
        if (b[i].n_childs == 0) {
            simulate_dir_int(&points_index[b[i].p_start], b[i].n_p,
                             &points_index[b[i].p_start], b[i].n_p,
                             &coord[0][0], nc, true, compute, n);
            for (int j1 = 1; j1 <= b[i].n_l1; ++j1) {
                int j = b[i].l1[j1];
                simulate_dir_int(&points_index[b[i].p_start], b[i].n_p,
                                 &points_index[b[j].p_start], b[j].n_p,
                                 &coord[0][0], nc, false, compute, n);
            }
        }
    }
    return n;
}

// ===========================================================================
namespace {
void afmm_step1(const std::vector<std::vector<double>>& coord, int nc, int n,
                const std::vector<double>& q, const TrickyFn& tricky) {
    (void)n;
    for (int i = 1; i <= nbox; ++i) {
        b[i].phi = CArr{};
        b[i].psi = CArr{};
    }
    for (int i = 1; i <= nbox; ++i) {
        if (b[i].n_childs == 0) {
            if (tricky) {
                tricky(b[i].x0, b[i].y0, b[i].z0, &points_index[b[i].p_start],
                       b[i].n_p, pmn, y_norm, &coord[0][0], nc, q.data(),
                       b[i].phi, P);
            } else {
                normal_multipoles(b[i].x0, b[i].y0, b[i].z0,
                                  &points_index[b[i].p_start], b[i].n_p, pmn,
                                  y_norm, coord, nc, q, b[i].phi, P);
            }
        }
    }
}

void normal_multipoles(double x0, double y0, double z0, const int* ind_, int n,
                       RArr& pmn_, const RArr& ynorm_,
                       const std::vector<std::vector<double>>& c, int nc,
                       const std::vector<double>& q, CArr& psi, int p) {
    double dx, dy, dz, t, tt, r, cos_t, phi;
    std::complex<double> tc, tc1;
    int j, k, n1, m1;
    (void)nc;
    for (j = 1; j <= n; ++j) {
        k = ind_[j];
        dx = c[1][k] - x0;
        dy = c[2][k] - y0;
        dz = c[3][k] - z0;
        r = std::sqrt(dx * dx + dy * dy + dz * dz);
        cos_t = dz / r;
        phi = std::atan2(dy, dx);
        get_legendre(p, cos_t, pmn_);
        psi(0, 0) = psi(0, 0) + std::complex<double>(q[k], 0.0);
        t = 1.0;
        for (n1 = 1; n1 <= p; ++n1) {
            t = t * r;
            for (m1 = 0; m1 <= n1; ++m1) {
                tc = std::exp(std::complex<double>(0.0, -m1 * phi));
                tt = ynorm_(-m1, n1) * pmn_(-m1, n1) * t * q[k];
                tc1 = tc * tt;
                psi(m1, n1) = psi(m1, n1) + tc1;
            }
        }
        for (n1 = 1; n1 <= p; ++n1)
            for (m1 = 1; m1 <= n1; ++m1)
                psi(-m1, n1) = std::conj(psi(m1, n1));
    }
}

void afmm_step2() {
    double dx, dy, dz, tt, r, cos_t, phi;
    std::complex<double> tc, tc1;
    int i, j, k, j1, m, n1, m1, n2, m2;
    for (i = nlev - 1; i >= 2; --i) {
        for (j = level[i]; j <= level[i + 1] - 1; ++j) {
            for (k = 1; k <= b[j].n_childs; ++k) {
                j1 = b[j].childs[k];
                dx = b[j1].x0 - b[j].x0;
                dy = b[j1].y0 - b[j].y0;
                dz = b[j1].z0 - b[j].z0;
                r = std::sqrt(dx * dx + dy * dy + dz * dz);
                cos_t = dz / r;
                phi = std::atan2(dy, dx);
                get_legendre(P, cos_t, pmn);
                for (n1 = 0; n1 <= P; ++n1) {
                    for (m1 = 0; m1 <= n1; ++m1) {
                        tc = {0.0, 0.0};
                        for (n2 = 0; n2 <= n1; ++n2) {
                            for (m2 = -n2; m2 <= n2; ++m2) {
                                if (std::abs(m1 - m2) > n1 - n2) continue;
                                m = std::abs(m1) - std::abs(m2) - std::abs(m1 - m2);
                                tt = amn(m2, n2) * amn(m1 - m2, n1 - n2) /
                                         amn(m1, n1) * std::pow(r, n2) *
                                         y_norm(-m2, n2) * pmn(-m2, n2);
                                tc1 = std::exp(std::complex<double>(0.0, -m2 * phi)) *
                                      b[j1].phi(m1 - m2, n1 - n2) * imag_pow(m);
                                tc = tc + tt * tc1;
                            }
                        }
                        b[j].phi(m1, n1) = b[j].phi(m1, n1) + tc;
                    }
                }
            }
            for (n1 = 1; n1 <= P; ++n1)
                for (m1 = 1; m1 <= n1; ++m1)
                    b[j].phi(-m1, n1) = std::conj(b[j].phi(m1, n1));
        }
    }
}

void afmm_step35() {
    double dx, dy, dz, tt, r, cos_t, phi;
    std::complex<double> tc, tc1;
    int i, j, j1, m, n1, m1, n2, m2;
    // Step 4: move multipole expansions to boxes in l2 list.
    for (i = level[2]; i <= nbox; ++i) {
        for (j1 = 1; j1 <= b[i].n_l2; ++j1) {
            j = b[i].l2[j1];
            dx = b[i].x0 - b[j].x0;
            dy = b[i].y0 - b[j].y0;
            dz = b[i].z0 - b[j].z0;
            r = std::sqrt(dx * dx + dy * dy + dz * dz);
            cos_t = dz / r;
            if (std::abs(dx) < 1.0e-6 && std::abs(dy) < 1.0e-6) phi = 0.0;
            else phi = std::atan2(dy, dx);
            get_legendre(P, cos_t, pmn);
            for (n1 = 0; n1 <= P; ++n1) {
                for (m1 = 0; m1 <= n1; ++m1) {
                    tc = {0.0, 0.0};
                    for (n2 = 0; n2 <= P; ++n2) {
                        if (n1 + n2 > P) continue;
                        for (m2 = -n2; m2 <= n2; ++m2) {
                            if (std::abs(m2 - m1) > n1 + n2) continue;
                            m = std::abs(m1 - m2) - std::abs(m1) - std::abs(m2);
                            tt = amn(m2, n2) * amn(m1, n1) /
                                     (amn(m2 - m1, n1 + n2) * std::pow(r, n1 + n2 + 1)) *
                                     y_norm(m2 - m1, n1 + n2) * pmn(m2 - m1, n1 + n2) /
                                     ((n2 & 1) ? -1.0 : 1.0);
                            tc1 = std::exp(std::complex<double>(0.0, (m2 - m1) * phi)) *
                                  b[i].phi(m2, n2) * imag_pow(m);
                            tc = tc + tt * tc1;
                        }
                    }
                    b[j].psi(m1, n1) = b[j].psi(m1, n1) + tc;
                }
            }
            for (n1 = 1; n1 <= P; ++n1)
                for (m1 = 1; m1 <= n1; ++m1)
                    b[j].psi(-m1, n1) = std::conj(b[j].psi(m1, n1));
        }
    }
    // Step 5: shift local expansion to children.
    for (i = level[2]; i <= nbox; ++i) {
        for (j1 = 1; j1 <= b[i].n_childs; ++j1) {
            j = b[i].childs[j1];
            dx = b[i].x0 - b[j].x0;
            dy = b[i].y0 - b[j].y0;
            dz = b[i].z0 - b[j].z0;
            r = std::sqrt(dx * dx + dy * dy + dz * dz);
            cos_t = dz / r;
            phi = std::atan2(dy, dx);
            get_legendre(P, cos_t, pmn);
            for (n1 = 0; n1 <= P; ++n1) {
                for (m1 = -n1; m1 <= n1; ++m1) {
                    tc = {0.0, 0.0};
                    for (n2 = n1; n2 <= P; ++n2) {
                        for (m2 = -n2; m2 <= n2; ++m2) {
                            if (std::abs(m2 - m1) > n2 - n1) continue;
                            m = std::abs(m2) - std::abs(m2 - m1) - std::abs(m1);
                            tt = amn(m2 - m1, n2 - n1) * amn(m1, n1) / amn(m2, n2) *
                                     std::pow(r, n2 - n1) * y_norm(m2 - m1, n2 - n1) *
                                     pmn(m2 - m1, n2 - n1) /
                                     (((n2 + n1) & 1) ? -1.0 : 1.0);
                            tc1 = std::exp(std::complex<double>(0.0, (m2 - m1) * phi)) *
                                  b[i].psi(m2, n2) * imag_pow(m);
                            tc = tc + tt * tc1;
                        }
                    }
                    b[j].psi(m1, n1) = b[j].psi(m1, n1) + tc;
                }
            }
        }
    }
}
}  // namespace

// ===========================================================================
void afmm(const std::vector<std::vector<double>>& coord, int nc, int n,
          const std::vector<double>& q, std::vector<double>& res, int res_size,
          const DirIntFn& dir_int, const MultIntFn& mult_int,
          const TrickyFn& tricky) {
    double t, tt;
    int i, j, j1;
    if (tricky) afmm_step1(coord, nc, n, q, tricky);
    else afmm_step1(coord, nc, n, q, nullptr);
    afmm_step2();
    afmm_step35();
    res.assign(res_size + 1, 0.0);
    t = 0;
    tt = n;
    tt = tt * (tt - 1.0) * 0.5;
    for (i = 1; i <= nbox; ++i) {
        if (b[i].n_childs == 0) {
            dir_int(&points_index[b[i].p_start], b[i].n_p,
                    &points_index[b[i].p_start], b[i].n_p,
                    &coord[0][0], nc, q.data(), res.data(), true);
            t = t + (b[i].n_p * (b[i].n_p - 1)) / 2.0;
            for (j1 = 1; j1 <= b[i].n_l1; ++j1) {
                j = b[i].l1[j1];
                dir_int(&points_index[b[i].p_start], b[i].n_p,
                        &points_index[b[j].p_start], b[j].n_p,
                        &coord[0][0], nc, q.data(), res.data(), false);
                t = t + b[i].n_p * b[j].n_p;
            }
        }
    }
    for (i = 1; i <= nbox; ++i) {
        if (b[i].n_childs == 0) {
            mult_int(&points_index[b[i].p_start], b[i].n_p, b[i].psi, P,
                     b[i].x0, b[i].y0, b[i].z0, y_norm, pmn,
                     &coord[0][0], nc, q.data(), res.data());
        }
    }
}

// ===========================================================================
void simple_mm(const std::vector<std::vector<double>>& coord, int nc, int n,
               const std::vector<double>& q, std::vector<double>& res,
               const DirIntFn& dir_int, const FarIntFn& far_int,
               const TrickyFn& tricky) {
    double t, tt;
    int i, j, j1;
    res.assign(n + 1, 0.0);
    if (tricky) afmm_step1(coord, nc, n, q, tricky);
    else afmm_step1(coord, nc, n, q, nullptr);
    afmm_step2();
    for (i = level[2]; i <= nbox; ++i) {
        for (j1 = 1; j1 <= b[i].n_l2; ++j1) {
            j = b[i].l2[j1];
            far_int(&points_index[b[i].p_start], b[i].n_p, b[j].phi, P,
                    b[j].x0, b[j].y0, b[j].z0, y_norm, pmn,
                    &coord[0][0], nc, q.data(), res.data());
        }
    }
    t = 0;
    tt = n;
    tt = tt * (tt - 1.0) * 0.5;
    for (i = 1; i <= nbox; ++i) {
        if (b[i].n_childs == 0) {
            dir_int(&points_index[b[i].p_start], b[i].n_p,
                    &points_index[b[i].p_start], b[i].n_p,
                    &coord[0][0], nc, q.data(), res.data(), true);
            t = t + (b[i].n_p * (b[i].n_p - 1)) / 2.0;
            for (j1 = 1; j1 <= b[i].n_l1; ++j1) {
                j = b[i].l1[j1];
                dir_int(&points_index[b[i].p_start], b[i].n_p,
                        &points_index[b[j].p_start], b[j].n_p,
                        &coord[0][0], nc, q.data(), res.data(), false);
                t = t + b[i].n_p * b[j].n_p;
            }
        }
    }
}

// ===========================================================================
namespace {
void split_box(int parent, int level_, double x0, double y0, double z0, double d_,
               const std::vector<std::vector<double>>& coord, int nc, int n,
               int b_pos, int ind_pos, int& m) {
    int n1, n2, n3, n4, n5, n6;
    m = 0;
    split_index(1, coord, nc, x0, ind_pos, n, n1, n2);
    if (n1 > 0) {
        split_index(2, coord, nc, y0, ind_pos, n1, n3, n4);
        if (n3 > 0) {
            split_index(3, coord, nc, z0, ind_pos, n3, n5, n6);
            if (n5 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 - d_; b[b_pos + m].y0 = y0 - d_; b[b_pos + m].z0 = z0 - d_;
                b[b_pos + m].p_start = 1; b[b_pos + m].n_p = n5;
            }
            if (n6 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 - d_; b[b_pos + m].y0 = y0 - d_; b[b_pos + m].z0 = z0 + d_;
                b[b_pos + m].p_start = 1 + n5; b[b_pos + m].n_p = n6;
            }
        }
        if (n4 > 0) {
            split_index(3, coord, nc, z0, ind_pos + n3, n4, n5, n6);
            if (n5 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 - d_; b[b_pos + m].y0 = y0 + d_; b[b_pos + m].z0 = z0 - d_;
                b[b_pos + m].p_start = n3 + 1; b[b_pos + m].n_p = n5;
            }
            if (n6 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 - d_; b[b_pos + m].y0 = y0 + d_; b[b_pos + m].z0 = z0 + d_;
                b[b_pos + m].p_start = n3 + 1 + n5; b[b_pos + m].n_p = n6;
            }
        }
    }
    if (n2 > 0) {
        split_index(2, coord, nc, y0, ind_pos + n1, n2, n3, n4);
        if (n3 > 0) {
            split_index(3, coord, nc, z0, ind_pos + n1, n3, n5, n6);
            if (n5 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 + d_; b[b_pos + m].y0 = y0 - d_; b[b_pos + m].z0 = z0 - d_;
                b[b_pos + m].p_start = n1 + 1; b[b_pos + m].n_p = n5;
            }
            if (n6 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 + d_; b[b_pos + m].y0 = y0 - d_; b[b_pos + m].z0 = z0 + d_;
                b[b_pos + m].p_start = n1 + 1 + n5; b[b_pos + m].n_p = n6;
            }
        }
        if (n4 > 0) {
            split_index(3, coord, nc, z0, ind_pos + n1 + n3, n4, n5, n6);
            if (n5 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 + d_; b[b_pos + m].y0 = y0 + d_; b[b_pos + m].z0 = z0 - d_;
                b[b_pos + m].p_start = n1 + n3 + 1; b[b_pos + m].n_p = n5;
            }
            if (n6 > 0) {
                m = m + 1;
                b[b_pos + m].x0 = x0 + d_; b[b_pos + m].y0 = y0 + d_; b[b_pos + m].z0 = z0 + d_;
                b[b_pos + m].p_start = n1 + n3 + 1 + n5; b[b_pos + m].n_p = n6;
            }
        }
    }
    for (int i = b_pos + 1; i <= b_pos + m; ++i) {
        b[i].parent = parent;
        b[i].level = level_;
        b[i].n_childs = 0;
        b[i].n_l1 = 0;
        b[i].n_l2 = 0;
        b[i].phi = CArr{};
        b[i].psi = CArr{};
    }
}

void split_index(int k, const std::vector<std::vector<double>>& coord, int /*nc*/,
                 double g, int ind_pos, int n, int& n1, int& n2) {
    n1 = 0;
    n2 = n + 1;
    for (int i = 1; i <= n; ++i) {
        if (coord[k - 1][ind[ind_pos + i]] < g) {
            n1 = n1 + 1;
            ind1[n1] = ind[ind_pos + i];
        } else {
            n2 = n2 - 1;
            ind1[n2] = ind[ind_pos + i];
        }
    }
    n2 = n - n1;
    for (int i = 1; i <= n; ++i) ind[ind_pos + i] = ind1[i];
}
}  // namespace

// ===========================================================================
void get_legendre(int n_max, double& x, RArr& p_pmn) {
    double x2, t, s, s1;
    int n, m;
    x2 = x * x;
    if (x2 >= 1.0) {
        x = (x > 0) ? 1.0 : -1.0;
        x2 = 1.0;
        t = 0.0;
    } else {
        t = std::sqrt(1.0 - x2);
    }
    p_pmn(0, 0) = 1.0;
    p_pmn(0, 1) = x;
    p_pmn(1, 1) = -t;
    p_pmn(0, 2) = 0.5 * (3.0 * x2 - 1.0);
    p_pmn(1, 2) = -3.0 * x * t;
    p_pmn(2, 2) = 3.0 * t * t;
    s = p_pmn(2, 2);
    for (n = 3; n <= n_max; ++n) {
        s = s * (2 * n - 1) * t;
        p_pmn(n, n) = s * (((n & 1) ? -1.0 : 1.0));
    }
    for (n = 3; n <= n_max; ++n) {
        for (m = n - 1; m >= 0; --m) {
            s = 1.0 / (double)(n - m);
            s1 = x * (2.0 * n - 1.0) * p_pmn(m, n - 1);
            if (m <= n - 2) s1 = s1 - (double)(n + m - 1) * p_pmn(m, n - 2);
            p_pmn(m, n) = s1 * s;
        }
    }
    for (n = 1; n <= n_max; ++n)
        for (m = 1; m <= n; ++m) p_pmn(-m, n) = p_pmn(m, n);
}

// ===========================================================================
void afmm_ini() {
    double t = 1.0;
    fact[0] = t;
    fact[1] = t;
    for (int i = 2; i <= 2 * P; ++i) {
        t = t * i;
        fact[i] = t;
    }
    y_norm(0, 0) = 1.0;
    for (int i = 1; i <= P; ++i) {
        for (int j = 0; j <= i; ++j) {
            t = std::sqrt(fact[i - j] / fact[i + j]);
            y_norm(j, i) = t;
            y_norm(-j, i) = t;
        }
    }
    t = 1.0;
    amn(0, 0) = 1.0;
    for (int i = 1; i <= P; ++i) {
        t = -t;
        for (int j = -i; j <= i; ++j)
            amn(j, i) = t / std::sqrt(fact[i - j] * fact[i + j]);
    }
}

}  // namespace afmm_mod
