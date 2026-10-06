// afmm_mod.h — C++ translation of MOPAC 2016 "afmm_mod.F90".
//
// Adaptive fast multipole method (Cheng/Greengard/Rokhlin 1999), with the
// spherical-harmonic addition theorem.  Kept as a namespace mirroring the
// Fortran module.  The multipole arrays phi/psi are indexed m in [-p,p],
// n in [0,p]; a small offset-wrapper maps Fortran subscripts to storage.
#pragma once

#include <complex>
#include <functional>
#include <vector>

namespace afmm_mod {

constexpr int P = 3;  // multipole expansion level (Fortran parameter p)

// Real array indexed [-P..P][0..P] (Fortran pmn / y_norm / amn).
struct RArr {
    double d[2 * P + 1][P + 1];
    RArr() { for (auto& row : d) for (auto& x : row) x = 0.0; }
    double& operator()(int m, int n) { return d[m + P][n]; }
    double operator()(int m, int n) const { return d[m + P][n]; }
};

// Complex array indexed [-P..P][0..P] (Fortran phi / psi).
struct CArr {
    std::complex<double> d[2 * P + 1][P + 1];
    CArr() { for (auto& row : d) for (auto& x : row) x = {0.0, 0.0}; }
    std::complex<double>& operator()(int m, int n) { return d[m + P][n]; }
    std::complex<double> operator()(int m, int n) const { return d[m + P][n]; }
};

struct Box {
    int parent = 0, level = 0;
    int n_p = 0, p_start = 0;
    int n_childs = 0, n_l1 = 0, n_l2 = 0;
    int n_up = 0, n_down = 0, n_north = 0, n_south = 0, n_east = 0, n_west = 0;
    int childs[9] = {};   // Fortran dimension(8), 1-based
    int l1[29] = {};      // dimension(28)
    int l2[190] = {};     // dimension(189)
    double x0 = 0, y0 = 0, z0 = 0;
    CArr phi, psi;
};

struct Tess {
    std::vector<Box> b;
    std::vector<double> d;
    std::vector<int> level, ind, ind1;
    int nbox = 0, nlev = 0, n = 0;
};

// Callback procedure-argument types (raw pointers match Fortran assumed-size
// arrays; coord is column-major with leading dimension nc).
using SimulateDirIntFn = std::function<void(const int* ind1, int n1,
                                           const int* ind2, int n2,
                                           const double* c, int nc,
                                           bool same, bool compute, int& n)>;

using DirIntFn = std::function<void(const int* ind1, int n1,
                                     const int* ind2, int n2,
                                     const double* c, int nc,
                                     const double* q, double* r, bool same)>;

using MultIntFn = std::function<void(const int* ind, int n,
                                     const CArr& psi, int p,
                                     double x0, double y0, double z0,
                                     const RArr& y_norm, RArr& pmn,
                                     const double* c, int n3,
                                     const double* q, double* r)>;

using FarIntFn = std::function<void(const int* ind, int n,
                                    const CArr& phi, int p,
                                    double x0, double y0, double z0,
                                    const RArr& y_norm, RArr& pmn,
                                    const double* c, int nc,
                                    const double* q, double* r)>;

using TrickyFn = std::function<void(double x0, double y0, double z0,
                                    const int* ind, int n,
                                    RArr& pmn, const RArr& y_norm,
                                    const double* c, int nc, const double* q,
                                    CArr& psi, int p)>;

// ---- public API (Fortran "public" list) ----
void afmm_ini();
void prepare_tesselations(int n1, int n2, std::vector<int>& ind, int n, int& ierr);
void set_tesselation(int handle, int& ierr);
void divide_box(const std::vector<std::vector<double>>& coord, int n1, int n,
                double min_d, int min_p, int max_tes, int handle, int& ierr);
int count_short_ints(const std::vector<std::vector<double>>& coord, int nc,
                     bool compute, const SimulateDirIntFn& simulate_dir_int);
void afmm(const std::vector<std::vector<double>>& coord, int nc, int n,
          const std::vector<double>& q, std::vector<double>& res, int res_size,
          const DirIntFn& dir_int, const MultIntFn& mult_int,
          const TrickyFn& tricky = nullptr);
void simple_mm(const std::vector<std::vector<double>>& coord, int nc, int n,
               const std::vector<double>& q, std::vector<double>& res,
               const DirIntFn& dir_int, const FarIntFn& far_int,
               const TrickyFn& tricky = nullptr);
void get_legendre(int n_max, double& x, RArr& pmn);

}  // namespace afmm_mod
