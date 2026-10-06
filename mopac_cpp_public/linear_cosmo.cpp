// linear_cosmo.cpp — C++ translation of MOPAC 2016 "linear_cosmo.F90".
// Batch A: module data (cosmo_mini + linear_cosmo), ini_linear_cosmo,
// bpnew_vec, bz_vec, get_bvec, mult_triangle_vec, some_norm.
// F90 1-based logic maps directly onto the C++ vectors (0 slot unused).
#include "linear_cosmo.h"
#include "afmm_mod.h"
#include "cosmo.h"
#include "molkst_C.h"
#include "cosmo_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "overlaps_C.h"
#include "ijbo.h"
#include "mopend.h"
#include "mkl_bits.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>

using namespace molkst_C;
using namespace cosmo_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;

namespace linear_cosmo {

// ---- module cosmo_mini data ----
bool new_surface = false, new_iteration = false;
std::vector<double> a_block, r_vec, p_vec, q_vec, z_vec, a_diag, m_vec, a_part;
std::vector<double> tm;
std::vector<int> iblock_pos;

// ---- module linear_cosmo data ----
double c_proc = 0.0;
std::vector<double> rsc;
std::vector<int> nipsrs, nset, npoints, iatom_pos, ijbo_diag;
int max_block_size = 0, atom_handle = 0, surface_handle = 0;

// F90 ini_linear_cosmo @71: allocate Linear-Cosmo specific arrays, build
// iatom_pos / ijbo_diag, reset solv_energy, init AFMM.
void ini_linear_cosmo() {
    int maxrs = 70 * numat;
    if (!rsc.empty()) {
        rsc.clear(); tm.clear(); iatom_pos.clear(); nipsrs.clear();
        nset.clear(); npoints.clear(); ijbo_diag.clear();
    }
    rsc.assign(4 * maxrs, 0.0);
    tm.assign(3 * 3 * numat, 0.0);
    iatom_pos.assign(numat + 1, 0);
    nipsrs.assign(lenabc + 1, 0);
    nset.assign(1082 * numat + 1, 0);
    npoints.assign(numat + 2, 0);
    ijbo_diag.assign(numat + 1, 0);

    // F90 1-based: iatom_pos(1)=0; ijbo_diag(1)=ijbo(1,1)
    iatom_pos[1] = 0;
    ijbo_diag[1] = ijbo(1, 1);
    for (int i = 2; i <= numat; ++i) {
        int nnn = nlast[i - 1] - nfirst[i - 1] + 1;
        nnn = ((nnn + 1) * nnn) / 2;
        iatom_pos[i] = iatom_pos[i - 1] + nnn;
        ijbo_diag[i] = ijbo(i, i);
    }
    solv_energy = 0.0;
    afmm_mod::afmm_ini();
}

// F90 bpnew_vec @116: interaction of electronic densities with induced
// charges.  v(i) -= <density j at surface point i> over all atoms.
// p is the packed density matrix (1-based indexing via jj=ijbo_diag(j)).
void bpnew_vec(std::vector<double>& v) {
    double w[46];   // F90 w(45), 1-based
    if ((int)v.size() < nps + 1) v.assign(nps + 1, 0.0);
    for (int i = 1; i <= nps; ++i) v[i] = 0.0;
    for (int j = 1; j <= numat; ++j) {
        int nao = nlast[j] - nfirst[j];
        int jj = ijbo_diag[j];
        int ni = 0;
        double dip = 0.0, quad = 0.0;
        if (nao > 0) {
            ni = nat[j];
            dip = dd[ni] * a0;
            quad = (a0 * qq[ni]) * (a0 * qq[ni]);
        }
        for (int i = 1; i <= nps; ++i) {
            double dx = cosurf[1][i] - coord[0][j];
            double dy = cosurf[2][i] - coord[1][j];
            double dz = cosurf[3][i] - coord[2][j];
            double r = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
            if (nao == 0) {          // --- S-element
                v[i] -= p[jj + 1] * r;
                continue;
            }
            // --- sp element
            w[1] = r;
            dx *= r; dy *= r; dz *= r;
            double t = dip * r * r;
            w[2] = dx * t; w[4] = dy * t; w[7] = dz * t;
            t = quad * r * r * r;
            w[3] = r + (3.0 * dx * dx - 1.0) * t;
            w[6] = r + (3.0 * dy * dy - 1.0) * t;
            w[10] = r + (3.0 * dz * dz - 1.0) * t;
            w[5] = 3.0 * dx * dy * t;
            w[8] = 3.0 * dx * dz * t;
            w[9] = 3.0 * dy * dz * t;
            if (nao > 3) {           // --- d element
                for (int k = 11; k <= 45; ++k) w[k] = 0.0;
                w[15] = r; w[21] = r; w[28] = r; w[36] = r; w[45] = r;
            }
            int m = 0;
            t = 0.0;
            for (int k = 1; k <= nao + 1; ++k) {
                for (int l = 1; l <= k - 1; ++l) {
                    ++m;
                    t -= 2 * p[jj + m] * w[m];
                }
                ++m;
                t -= p[jj + m] * w[m];
            }
            v[i] += t;
        }
    }
}

// F90 bz_vec @208: interaction of atom cores with induced charges.
void bz_vec(std::vector<double>& v) {
    if ((int)v.size() < nps + 1) v.assign(nps + 1, 0.0);
    for (int i = 1; i <= nps; ++i) {
        double t = 0.0;
        for (int j = 1; j <= numat; ++j) {
            double dx = coord[0][j] - cosurf[1][i];
            double dy = coord[1][j] - cosurf[2][i];
            double dz = coord[2][j] - cosurf[3][i];
            t += tore[nat[j]] / std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        v[i] = t;
    }
}

// F90 get_bvec @239: multipole expansion of 1/r in the local frame of atom
// (x1 = surface point, x2 = atom center).  w(1..45) packed multipoles.
void get_bvec(const double* x1, const double* x2, int nao, int ni, double* w) {
    double dx[3];
    dx[0] = x1[0] - x2[0];
    dx[1] = x1[1] - x2[1];
    dx[2] = x1[2] - x2[2];
    double r = 1.0 / std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    w[1] = r;
    if (nao == 0) return;            // --- H
    // --- sp element
    double dip = dd[ni] * a0 * r * r;
    double quad = (a0 * qq[ni]) * (a0 * qq[ni]) * r * r * r;
    dx[0] *= r; dx[1] *= r; dx[2] *= r;
    w[2] = dx[0] * dip;
    w[4] = dx[1] * dip;
    w[7] = dx[2] * dip;
    w[3] = r + (3.0 * dx[0] * dx[0] - 1.0) * quad;
    w[6] = r + (3.0 * dx[1] * dx[1] - 1.0) * quad;
    w[10] = r + (3.0 * dx[2] * dx[2] - 1.0) * quad;
    w[5] = 3.0 * dx[0] * dx[1] * quad;
    w[8] = 3.0 * dx[0] * dx[2] * quad;
    w[9] = 3.0 * dx[1] * dx[2] * quad;
    if (nao < 4) return;             // --- no d element
    for (int i = 11; i <= 45; ++i) w[i] = 0.0;
    w[15] = r; w[21] = r; w[28] = r; w[36] = r; w[45] = r;
}

// F90 mult_triangle_vec @2507: y = T*x with T a full packed (lower) triangle
// (row-major packed order: (1,1),(2,1),(2,2),(3,1),...).
void mult_triangle_vec(const double* t, const double* x, int n, double* y) {
    int ii = 0;   // previous diagonal element
    for (int i = 1; i <= n; ++i) {
        double s = 0.0;
        int ij = ii;
        for (int j = 1; j <= i; ++j) { ++ij; s += t[ij] * x[j]; }
        ii = ij;
        for (int j = i + 1; j <= n; ++j) { ij = ij + j - 1; s += t[ij] * x[j]; }
        y[i] = s;
    }
}

// F90 some_norm @2229: max abs of vector (infinity norm).
double some_norm(const double* v, int n) {
    double t = std::fabs(v[1]);
    for (int i = 2; i <= n; ++i)
        if (t < std::fabs(v[i])) t = std::fabs(v[i]);
    return t;
}

// F90 coscavz @295: driver — build/refresh the SAS via coscanz, set up the
// AFMM tessellations and allocate the A-part / preconditioner storage.
void coscavz() {
    new_surface = true;
    coscanz();
    std::vector<int> ind(3, 0);
    int i = 0;
    afmm_mod::prepare_tesselations(numat, nps, ind, 2, i);
    if (i != 0) { mopend("Error in divide box by atoms "); return; }
    atom_handle = ind[1];
    surface_handle = ind[2];
    if (numat > nb2max || numat > nb1max) {
        afmm_mod::divide_box(coord, 3, numat, std::sqrt(overlaps_C::cutof2), 10, 0,
                             atom_handle, i);
    }
    if (i != 0) { mopend("Error in divide box by atoms "); return; }
    if (nps > na2max || nps > na1max) {
        afmm_mod::divide_box(cosurf, 4, nps, std::sqrt(disex2), 20, 0,
                             surface_handle, i);
    }
    if (i != 0) { mopend("Error in divide box by surface "); return; }
    c_proc = 0;

    if (compute_a_part) {
        if (!a_part.empty()) a_part.clear();
        int ndim = 0;
        if (nps > na2max || nps > na1max) {
            afmm_mod::set_tesselation(surface_handle, i);
            if (i != 0) { mopend("CosmoZ:  Tesselation ini error "); return; }
            ndim = afmm_mod::count_short_ints(cosurf, 4, false, simulate_aq_dir_int);
        } else {
            int maxrs = 60 * numat;
            simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                            nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc,
                            false, ndim);
        }
        a_part.assign(ndim, 0.0);
    }

    if (!m_vec.empty()) m_vec.clear();
    if (use_a_blocks) {
        if (!iblock_pos.empty()) { iblock_pos.clear(); a_block.clear(); }
        iblock_pos.assign(numat + 1, 0);
        int ndim = 0;
        max_block_size = 0;
        for (int i1 = 1; i1 <= numat; ++i1) {
            iblock_pos[i1] = ndim + 1;
            int j = npoints[i1 + 1] - npoints[i1];
            if (j > max_block_size) max_block_size = j;
            ndim += (j * (j + 1)) / 2;
        }
        m_vec.assign(ndim + 1, 0.0);
        a_block.assign((max_block_size + 1) * (max_block_size + 1), 0.0);
    } else {
        m_vec.assign(nps + 1, 0.0);
    }
}
// F90 coscavz(coord, nat): flat coordinate vector (length 3*nat, layout
// c[3*(i-1)+k] = coord(k,i)) -> stashed into global coord, then the no-arg
// driver runs.  Callers: prtpka (mozyme branch).
void coscavz(std::vector<double>& c, int nat) {
    coord.assign(4, std::vector<double>(nat + 1, 0.0));
    for (int i = 1; i <= nat; ++i)
        for (int k = 1; k <= 3; ++k) coord[k-1][i] = c[3 * (i - 1) + k];
    coscavz();
}


// F90 coscanz @449: construct the solvent-accessible surface (linear-scaling
// version of the COSMO coscan).  Box partition, atom assignment, neighbour
// search, per-atom transformation matrix tm, direction grid, segment packing.
void coscanz() {
    constexpr int MAXDIR = 1082;
    // F90 1-based flat helpers for the direction tables (column-major storage:
    // dirvec(ix, j) -> dirvec[(j-1)*4 + (ix-1)]).
    auto DV = [](int ix, int j) -> double& { return (&dirvec[0][0])[(j - 1) * 4 + (ix - 1)]; };
    auto DS = [](int ix, int j) -> double& { return (&dirsm[0][0])[(j - 1) * 4 + (ix - 1)]; };

    // --- make coordinates a bit asymmetric to avoid symmetry problems ---
    for (int i = 1; i <= numat; ++i)
        for (int j = 1; j <= 3; ++j)
            coord[j - 1][i] += std::cos(i * j * 0.1) * 3.0e-9;

    // --- determine how many boxes ---
    double xmax = coord[0][1], ymax = coord[1][1], zmax = coord[2][1];
    double xmin = coord[0][1], ymin = coord[1][1], zmin = coord[2][1];
    double r_max = srad[1];
    for (int i = 2; i <= numat; ++i) {
        if (coord[0][i] > xmax) xmax = coord[0][i];
        else if (coord[0][i] < xmin) xmin = coord[0][i];
        if (coord[1][i] > ymax) ymax = coord[1][i];
        else if (coord[1][i] < ymin) ymin = coord[1][i];
        if (coord[2][i] > zmax) zmax = coord[2][i];
        else if (coord[2][i] < zmin) zmin = coord[2][i];
        if (srad[i] > r_max) r_max = srad[i];
    }
    r_max = r_max + rsolv;
    r_max = r_max + r_max;          // cell size
    double o_2r = 1.0 / r_max;

    int nx = (int)((xmax - xmin) * o_2r);
    if (nx < 1) nx = 1;
    if (xmin + nx * r_max < xmax) nx = nx + 1;
    int ny = (int)((ymax - ymin) * o_2r);
    if (ny < 1) ny = 1;
    if (ymin + ny * r_max < ymax) ny = ny + 1;
    int nz = (int)((zmax - zmin) * o_2r);
    if (nz < 1) nz = 1;
    if (zmin + nz * r_max < zmax) nz = nz + 1;

    int maxrs = 100 * numat;
    // F90 atoms_in_box(0:nx*ny*nz), s_box(0:numat-1), box_number(numat), nn(3,numat)
    std::vector<int> atoms_in_box(nx * ny * nz + 1, 0), s_box(numat, 0),
                     box_number(numat + 1, 0), isort(maxrs + 1, 0),
                     ipsrs(maxrs + 1, 0), nipa(numat + 1, 0);
    std::vector<int> nn(3 * (numat + 1), 0);   // nn(k, i) -> nn[(i-1)*3+(k-1)]
    std::vector<char> din(std::max(numat, MAXDIR) + 1, 0);

    // --- assign atoms to boxes ---
    for (int i = 1; i <= numat; ++i) {
        int ix = (int)((coord[0][i] - xmin) * o_2r);
        if (ix < 0) ix = 0; else if (ix >= nx) ix = nx - 1;
        int iy = (int)((coord[1][i] - ymin) * o_2r);
        if (iy < 0) iy = 0; else if (iy >= ny) iy = ny - 1;
        int iz = (int)((coord[2][i] - zmin) * o_2r);
        if (iz < 0) iz = 0; else if (iz >= nz) iz = nz - 1;
        box_number[i] = ix + (iy + iz * ny) * nx;
        ++atoms_in_box[box_number[i]];
    }
    // prefix sums
    for (int i = 0; i < nx * ny * nz; ++i)
        atoms_in_box[i + 1] += atoms_in_box[i];
    for (int i = 1; i <= numat; ++i) {
        --atoms_in_box[box_number[i]];
        s_box[atoms_in_box[box_number[i]]] = i;
    }

    // --- intersect: count near-neighbour atom pairs (pass 1) ---
    int ilipa = 0;
    int k = -ny;
    for (int iz = 0; iz < nz; ++iz) {
        k += ny;
        for (int iy = 0; iy < ny; ++iy) {
            int l = (iy + k) * nx;
            for (int ix = 0; ix < nx; ++ix) {
                int i1 = atoms_in_box[l];
                ++l;
                int i2 = atoms_in_box[l];
                for (int ii = i1; ii < i2; ++ii) {
                    int i = s_box[ii];
                    double ri = srad[i];
                    double r = ri + rsolv;
                    double rr = r + rsolv;
                    double xa1 = coord[0][i], xa2 = coord[1][i], xa3 = coord[2][i];
                    int izfrom = std::max(iz - 1, 0), izto = std::min(iz + 2, nz) - 1;
                    int k1 = (izfrom - 1) * ny;
                    for (int iz1 = izfrom; iz1 <= izto; ++iz1) {
                        int iyfrom = std::max(iy - 1, 0), iyto = std::min(iy + 2, ny) - 1;
                        k1 += ny;
                        for (int iy1 = iyfrom; iy1 <= iyto; ++iy1) {
                            int l1 = (iy1 + k1) * nx;
                            int ixfrom = std::max(ix - 1, 0), ixto = std::min(ix + 2, nx) - 1;
                            l1 += ixfrom - 1;
                            for (int ix1 = ixfrom; ix1 <= ixto; ++ix1) {
                                ++l1;
                                for (int jj = atoms_in_box[l1]; jj < atoms_in_box[l1 + 1]; ++jj) {
                                    int j = s_box[jj];
                                    if (j != i) {
                                        double dist = (xa1 - coord[0][j]) * (xa1 - coord[0][j]) +
                                                     (xa2 - coord[1][j]) * (xa2 - coord[1][j]) +
                                                     (xa3 - coord[2][j]) * (xa3 - coord[2][j]);
                                        if (dist < (rr + srad[j]) * (rr + srad[j])) ++ilipa;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

        // pass 2: fill lipa/temp keeping atoms in original (s_box) order
    std::vector<int> temp(ilipa + 1, 0), lipa(ilipa + 1, 0);
    ilipa = 0;
    k = -ny;
    for (int iz = 0; iz < nz; ++iz) {
        k += ny;
        for (int iy = 0; iy < ny; ++iy) {
            int l = (iy + k) * nx;
            for (int ix = 0; ix < nx; ++ix) {
                int i1 = atoms_in_box[l];
                ++l;
                int i2 = atoms_in_box[l];
                for (int ii = i1; ii < i2; ++ii) {
                    int i = s_box[ii];
                    nipa[i] = 0;
                    double ri = srad[i];
                    double r = ri + rsolv;
                    double rr = r + rsolv;
                    double xa1 = coord[0][i], xa2 = coord[1][i], xa3 = coord[2][i];
                    int izfrom = std::max(iz - 1, 0), izto = std::min(iz + 2, nz) - 1;
                    int k1 = (izfrom - 1) * ny;
                    for (int iz1 = izfrom; iz1 <= izto; ++iz1) {
                        int iyfrom = std::max(iy - 1, 0), iyto = std::min(iy + 2, ny) - 1;
                        k1 += ny;
                        for (int iy1 = iyfrom; iy1 <= iyto; ++iy1) {
                            int l1 = (iy1 + k1) * nx;
                            int ixfrom = std::max(ix - 1, 0), ixto = std::min(ix + 2, nx) - 1;
                            l1 += ixfrom - 1;
                            for (int ix1 = ixfrom; ix1 <= ixto; ++ix1) {
                                ++l1;
                                for (int jj = atoms_in_box[l1]; jj < atoms_in_box[l1 + 1]; ++jj) {
                                    int j = s_box[jj];
                                    if (j != i) {
                                        double dist = (xa1 - coord[0][j]) * (xa1 - coord[0][j]) +
                                                     (xa2 - coord[1][j]) * (xa2 - coord[1][j]) +
                                                     (xa3 - coord[2][j]) * (xa3 - coord[2][j]);
                                        if (dist < (rr + srad[j]) * (rr + srad[j])) {
                                            // keep identical to original MOPAC COSMO:
                                            // insert j into lipa sorted descending w.r.t. current atom
                                            int i0 = ilipa + 1;
                                            for (int ik = ilipa; ik >= ilipa - nipa[i] + 1; --ik) {
                                                if (lipa[ik] > j) { lipa[i0] = lipa[ik]; --i0; }
                                                else break;
                                            }
                                            ++nipa[i];
                                            ++ilipa;
                                            lipa[i0] = j;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // --- reorder lipa to correspond to atoms in order 1..n ---
    box_number[1] = 0;
    for (int i = 2; i <= numat; ++i)
        box_number[i] = box_number[i - 1] + nipa[i - 1];
    int indi = 0;
    for (int ii = 0; ii < numat; ++ii) {
        int i = s_box[ii];
        for (int jj = 1; jj <= nipa[i]; ++jj)
            temp[box_number[i] + jj] = lipa[indi + jj];
        indi += nipa[i];
    }
    for (int i1 = 1; i1 <= indi; ++i1) lipa[i1] = temp[i1];

        // --- main loop ---
    int inset = 1;
    ilipa = 0;
    nps = 0;
    area = 0.0;
    cosvol = 0.0;

    double xa[4], xx[4], dirtm[3 * MAXDIR + 1];  // dirtm(ix,j) -> dirtm[(j-1)*3+(ix-1)]
    int iseg[MAXDIR + 1];
    for (int i = 1; i <= numat; ++i) {
        int nps0 = nps + 1;
        double ri = srad[i];
        double r = ri + rsolv;
        double rr = r + rsolv;
        double ri2 = ri * ri;
        xa[1] = coord[0][i]; xa[2] = coord[1][i]; xa[3] = coord[2][i];

        // --- search for 3 nearest neighbour atoms ---
        double dist1 = 1.e20, dist2 = 1.e20, dist3 = 1.e20;
        nn[(i - 1) * 3 + 0] = 0; nn[(i - 1) * 3 + 1] = 0; nn[(i - 1) * 3 + 2] = 0;
        for (int jj = 1; jj <= nipa[i]; ++jj) {
            int j = lipa[ilipa + jj];
            double dist = (xa[1] - coord[0][j]) * (xa[1] - coord[0][j]) +
                          (xa[2] - coord[1][j]) * (xa[2] - coord[1][j]) +
                          (xa[3] - coord[2][j]) * (xa[3] - coord[2][j]);
            if (dist + 0.05 < dist3) { dist3 = dist; nn[(i - 1) * 3 + 2] = j; }
            if (dist3 + 0.05 < dist2) {
                dist = dist2; dist2 = dist3; dist3 = dist;
                nn[(i - 1) * 3 + 2] = nn[(i - 1) * 3 + 1];
                nn[(i - 1) * 3 + 1] = j;
            }
            if (dist2 + 0.05 < dist1) {
                dist = dist1; dist1 = dist2; dist2 = dist;
                nn[(i - 1) * 3 + 1] = nn[(i - 1) * 3 + 0];
                nn[(i - 1) * 3 + 0] = j;
            }
        }
        ilipa += nipa[i];

        // --- build new transformation matrix ---
        // tm(ix,jy,i) -> tm[(i-1)*9 + (ix-1)*3 + (jy-1)]
        if (nn[(i - 1) * 3 + 0] == 0) {
            tm[(i - 1) * 9 + 0 * 3 + 0] = 1.0;
            tm[(i - 1) * 9 + 0 * 3 + 1] = 0.0;
            tm[(i - 1) * 9 + 0 * 3 + 2] = 0.0;
        } else {
            int n1 = nn[(i - 1) * 3 + 0];
            dist1 = (xa[1] - coord[0][n1]) * (xa[1] - coord[0][n1]) +
                    (xa[2] - coord[1][n1]) * (xa[2] - coord[1][n1]) +
                    (xa[3] - coord[2][n1]) * (xa[3] - coord[2][n1]);
            double d1 = 1.0 / std::sqrt(dist1);
            tm[(i - 1) * 9 + 0 * 3 + 0] = (coord[0][n1] - xa[1]) * d1;
            tm[(i - 1) * 9 + 0 * 3 + 1] = (coord[1][n1] - xa[2]) * d1;
            tm[(i - 1) * 9 + 0 * 3 + 2] = (coord[2][n1] - xa[3]) * d1;
        }
        for (;;) {
            if (nn[(i - 1) * 3 + 1] == 0) {
                tm[(i - 1) * 9 + 1 * 3 + 0] = -tm[(i - 1) * 9 + 0 * 3 + 1];
                tm[(i - 1) * 9 + 1 * 3 + 1] = tm[(i - 1) * 9 + 0 * 3 + 0];
                tm[(i - 1) * 9 + 1 * 3 + 2] = 0.0;
                break;
            } else {
                int n2 = nn[(i - 1) * 3 + 1];
                dist2 = (xa[1] - coord[0][n2]) * (xa[1] - coord[0][n2]) +
                        (xa[2] - coord[1][n2]) * (xa[2] - coord[1][n2]) +
                        (xa[3] - coord[2][n2]) * (xa[3] - coord[2][n2]);
                double d2 = 1.0 / std::sqrt(dist2);
                xx[1] = (coord[0][n2] - xa[1]) * d2;
                xx[2] = (coord[1][n2] - xa[2]) * d2;
                xx[3] = (coord[2][n2] - xa[3]) * d2;
                double sp = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 0] +
                            xx[2] * tm[(i - 1) * 9 + 0 * 3 + 1] +
                            xx[3] * tm[(i - 1) * 9 + 0 * 3 + 2];
                if (sp * sp > 0.99) {
                    nn[(i - 1) * 3 + 1] = nn[(i - 1) * 3 + 2];
                    nn[(i - 1) * 3 + 2] = 0;
                    dist2 = dist3;
                } else {
                    double sininv = 1.0 / std::sqrt(1.0 - sp * sp);
                    tm[(i - 1) * 9 + 1 * 3 + 0] = (xx[1] - sp * tm[(i - 1) * 9 + 0 * 3 + 0]) * sininv;
                    tm[(i - 1) * 9 + 1 * 3 + 1] = (xx[2] - sp * tm[(i - 1) * 9 + 0 * 3 + 1]) * sininv;
                    tm[(i - 1) * 9 + 1 * 3 + 2] = (xx[3] - sp * tm[(i - 1) * 9 + 0 * 3 + 2]) * sininv;
                    break;
                }
            }
        }
        tm[(i - 1) * 9 + 2 * 3 + 0] = tm[(i - 1) * 9 + 0 * 3 + 1] * tm[(i - 1) * 9 + 1 * 3 + 2] -
                                      tm[(i - 1) * 9 + 1 * 3 + 1] * tm[(i - 1) * 9 + 0 * 3 + 2];
        tm[(i - 1) * 9 + 2 * 3 + 1] = tm[(i - 1) * 9 + 0 * 3 + 2] * tm[(i - 1) * 9 + 1 * 3 + 0] -
                                      tm[(i - 1) * 9 + 1 * 3 + 2] * tm[(i - 1) * 9 + 0 * 3 + 0];
        tm[(i - 1) * 9 + 2 * 3 + 2] = tm[(i - 1) * 9 + 0 * 3 + 0] * tm[(i - 1) * 9 + 1 * 3 + 1] -
                                      tm[(i - 1) * 9 + 1 * 3 + 0] * tm[(i - 1) * 9 + 0 * 3 + 1];

        // --- transform dirvec according to tm ---
        for (int j = 1; j <= MAXDIR; ++j) {
            xx[1] = DV(1, j); xx[2] = DV(2, j); xx[3] = DV(3, j);
            dirtm[(j - 1) * 3 + 0] = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 0] +
                                     xx[2] * tm[(i - 1) * 9 + 1 * 3 + 0] +
                                     xx[3] * tm[(i - 1) * 9 + 2 * 3 + 0];
            dirtm[(j - 1) * 3 + 1] = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 1] +
                                     xx[2] * tm[(i - 1) * 9 + 1 * 3 + 1] +
                                     xx[3] * tm[(i - 1) * 9 + 2 * 3 + 1];
            dirtm[(j - 1) * 3 + 2] = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 2] +
                                     xx[2] * tm[(i - 1) * 9 + 1 * 3 + 2] +
                                     xx[3] * tm[(i - 1) * 9 + 2 * 3 + 2];
        }

        // --- find the points of the basic grid on the SAS ---
        int narea = 0;
        for (int j = 1; j <= MAXDIR; ++j) {
            din[j] = 0;
            xx[1] = xa[1] + dirtm[(j - 1) * 3 + 0] * r;
            xx[2] = xa[2] + dirtm[(j - 1) * 3 + 1] * r;
            xx[3] = xa[3] + dirtm[(j - 1) * 3 + 2] * r;
            // only try those atoms intersecting atom i
            bool hidden = false;
            for (int ik = ilipa - nipa[i] + 1; ik <= ilipa; ++ik) {
                int kk = lipa[ik];
                double dist = (xx[1] - coord[0][kk]) * (xx[1] - coord[0][kk]) +
                              (xx[2] - coord[1][kk]) * (xx[2] - coord[1][kk]) +
                              (xx[3] - coord[2][kk]) * (xx[3] - coord[2][kk]);
                dist = std::sqrt(dist) - rsolv - srad[kk];
                if (dist < 0) { hidden = true; break; }
            }
            if (hidden) continue;
            ++narea;
            cosvol += ri2 * DV(4, j) * (dirtm[(j - 1) * 3 + 0] * xa[1] +
                                        dirtm[(j - 1) * 3 + 1] * xa[2] +
                                        dirtm[(j - 1) * 3 + 2] * xa[3] + ri);
            area += ri2 * DV(4, j);
            din[j] = 1;
        }

        if (narea != 0) {
            // if hydrogen use the smaller set of points, else the larger
            int i0 = 1;
            if (nat[i] == 1) i0 = 2;
            int jmax = n0[i0];
            i0 = (i0 - 1) * n0[1];
            for (int j = 1; j <= jmax; ++j) {
                ++nps;
                if (nps > lenabc) {
                    mopend("NPS IS GREATER THAN LENABC-USE SMALLER NSPA");
                    return;
                } else {
                    iatsp[nps] = i;
                    xx[1] = DS(1, i0 + j);
                    xx[2] = DS(2, i0 + j);
                    xx[3] = DS(3, i0 + j);
                    cosurf[1][nps] = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 0] +
                                     xx[2] * tm[(i - 1) * 9 + 1 * 3 + 0] +
                                     xx[3] * tm[(i - 1) * 9 + 2 * 3 + 0];
                    cosurf[2][nps] = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 1] +
                                     xx[2] * tm[(i - 1) * 9 + 1 * 3 + 1] +
                                     xx[3] * tm[(i - 1) * 9 + 2 * 3 + 1];
                    cosurf[3][nps] = xx[1] * tm[(i - 1) * 9 + 0 * 3 + 2] +
                                     xx[2] * tm[(i - 1) * 9 + 1 * 3 + 2] +
                                     xx[3] * tm[(i - 1) * 9 + 2 * 3 + 2];
                }
            }

            int niter = 0;
            for (;;) {
                ++niter;
                for (int ips = nps0; ips <= nps; ++ips) {
                    nar_csm[ips] = 0;
                    phinet[ips][1] = 0.0;
                    phinet[ips][2] = 0.0;
                    phinet[ips][3] = 0.0;
                }
                for (int j = 1; j <= MAXDIR; ++j) {
                    if (din[j]) {
                        double spm = -1.0;
                        double x1 = dirtm[(j - 1) * 3 + 0];
                        double x2 = dirtm[(j - 1) * 3 + 1];
                        double x3 = dirtm[(j - 1) * 3 + 2];
                        int ipm = 0;
                        for (int ips = nps0; ips <= nps; ++ips) {
                            double sp = x1 * cosurf[1][ips] + x2 * cosurf[2][ips] + x3 * cosurf[3][ips];
                            if (sp >= spm) { spm = sp * (1.0 + 1.e-14); ipm = ips; }
                        }
                        iseg[j] = ipm;
                        ++nar_csm[ipm];
                        phinet[ipm][1] += dirtm[(j - 1) * 3 + 0] * DV(4, j);
                        phinet[ipm][2] += dirtm[(j - 1) * 3 + 1] * DV(4, j);
                        phinet[ipm][3] += dirtm[(j - 1) * 3 + 2] * DV(4, j);
                    }
                }
                int ips = nps0 - 1;
                for (;;) {
                    ++ips;
                    while (nar_csm[ips] == 0) {
                        niter = 1;
                        --nps;
                        if (ips > nps) break;
                        for (int jps = ips; jps <= nps; ++jps) {
                            nar_csm[jps] = nar_csm[jps + 1];
                            phinet[jps][1] = phinet[jps + 1][1];
                            phinet[jps][2] = phinet[jps + 1][2];
                            phinet[jps][3] = phinet[jps + 1][3];
                        }
                    }
                    if (ips > nps) break;
                    double dists = phinet[ips][1] * phinet[ips][1] +
                                   phinet[ips][2] * phinet[ips][2] +
                                   phinet[ips][3] * phinet[ips][3];
                    dists = std::max(dists, 1.e-20);
                    double dd = 1.0 / std::sqrt(dists);
                    cosurf[1][ips] = phinet[ips][1] * dd;
                    cosurf[2][ips] = phinet[ips][2] * dd;
                    cosurf[3][ips] = phinet[ips][3] * dd;
                    if (ips >= nps) break;
                }
                if (niter >= 2) break;
            }

            // all segments finally defined; close-pack the basic grid points
            for (int ips = nps0; ips <= nps; ++ips) {
                nsetf[ips] = inset;
                inset += nar_csm[ips];
                nar_csm[ips] = 0;
                cosurf[4][ips] = 0.0;
                cosurf[1][ips] = cosurf[1][ips] * ri + xa[1];
                cosurf[2][ips] = cosurf[2][ips] * ri + xa[2];
                cosurf[3][ips] = cosurf[3][ips] * ri + xa[3];
            }
            for (int j = 1; j <= MAXDIR; ++j) {
                if (din[j]) {
                    int ipm = iseg[j];
                    int nara = nar_csm[ipm];
                    nset[nsetf[ipm] + nara] = j;
                    nar_csm[ipm] = nara + 1;
                    cosurf[4][ipm] += DV(4, j) * ri2;
                }
            }
        }
    }

    // --- construction for a single atom ends here ---
    for (int i = 1; i <= numat; ++i) din[i] = 1;
    for (int j = 1; j <= nps; ++j) din[iatsp[j]] = 0;
    for (int i = 0; i < numat; ++i) s_box[i] = 0;
    for (int i = 1; i <= nps; ++i) {
        int j = iatsp[i] - 1;
        s_box[j] += (int)cosurf[4][i];
    }
    for (int i = 1; i <= numat; ++i)
        s_box[i - 1] = (int)(s_box[i - 1] / (srad[i] * srad[i]) * (srad[i] + rsolv) * (srad[i] + rsolv));

        // --- closure of the concave regions of the surface ---
    if (ioldcv == 0) {
        int idim = std::max(numat, MAXDIR);
        // surclo expects a planar [3][numat+1] coord buffer and double** rsc
        std::vector<double> coord_flat(3 * (numat + 1), 0.0);
        for (int a = 1; a <= numat; ++a)
            for (int ix = 0; ix < 3; ++ix)
                coord_flat[ix * (numat + 1) + a] = coord[ix][a];
        std::vector<std::array<double, 4>> rsc_local(maxrs + 1);
        std::vector<double*> rsc_ptrs(maxrs + 1);
        for (int i2 = 1; i2 <= maxrs; ++i2) rsc_ptrs[i2] = &rsc_local[i2][0];
        surclo(coord_flat.data(), nipa.data(), lipa.data(),
               reinterpret_cast<const bool*>(din.data()), idim, rsc_ptrs.data(),
               isort.data(), ipsrs.data(), nipsrs.data(), nat.data(),
               srad.data(), maxrs);
    }

        cosvol = cosvol / 3;

    std::fill(npoints.begin(), npoints.begin() + numat + 1, 0);
    for (int i = 1; i <= nps; ++i) ++npoints[iatsp[i]];
    int jj = npoints[1];
    npoints[1] = 1;
    for (int i = 2; i <= numat; ++i) {
        int kk = npoints[i];
        npoints[i] = npoints[i - 1] + jj;
        jj = kk;
    }
    npoints[numat + 1] = npoints[numat] + jj;

    std::fill(arat.begin(), arat.begin() + numat + 1, 0.0);
    double sumcf4=0; for (int i=1;i<=nps;++i) sumcf4+=cosurf[4][i];
    for (int i = 1; i <= nps; ++i) arat[iatsp[i]] += cosurf[4][i];

    // --- reallocate the iteration vectors ---
    if (!r_vec.empty()) {
        if (p_vec[1] < -1.e7 && z_vec[1] > 1.e7) return;   // dummy use of p_vec and z_vec
        r_vec.clear(); p_vec.clear(); q_vec.clear(); z_vec.clear(); a_diag.clear();
    }
    r_vec.assign(nps + 1, 0.0);
    p_vec.assign(nps + 1, 0.0);
    q_vec.assign(std::max(numat, nps) + 1, 0.0);
    z_vec.assign(nps + 1, 0.0);
    a_diag.assign(nps + 1, 0.0);
}

// ---- batch C: A-part machinery (aq_vec / simulate_aq_vec / aq_dir_int /
// ---- simulate_aq_dir_int / mfinel) --------------------------------------

// F90 mfinel (cosmo.F90 @1097): generates the list of all basic grid points
// and ring-segments belonging to segment ips, transformed into the atomic
// frame (finel[k][1..4][1..nfl], k=1..2 channel).
static void mfinel(int ips, int k, double finel[2][4][301], int& nfl,
                   int ioldcv, int maxrs, int lenabc, int numat,
                   const double* tm9, const double* x, double r) {
    auto DV = [](int ix, int j) -> double& { return (&dirvec[0][0])[(j - 1) * 4 + (ix - 1)]; };
    (void)maxrs; (void)lenabc; (void)numat;
    nfl = 0;
    int nari = nar_csm[ips];
    for (int l = nsetf[ips]; l <= nsetf[ips] + nari - 1; ++l) {
        int idir = nset[l];
        ++nfl;
        double y[3];
        for (int ix = 0; ix < 3; ++ix) y[ix] = DV(ix + 1, idir) * r;
        finel[k - 1][0][nfl] = y[0] * tm9[0] + y[1] * tm9[3] + y[2] * tm9[6] + x[0];
        finel[k - 1][1][nfl] = y[0] * tm9[1] + y[1] * tm9[4] + y[2] * tm9[7] + x[1];
        finel[k - 1][2][nfl] = y[0] * tm9[2] + y[1] * tm9[5] + y[2] * tm9[8] + x[2];
        finel[k - 1][3][nfl] = DV(4, idir) * r * r;
    }
    if (ioldcv == 1) return;
    // associated ring segments
    int irs0 = (ips > 1) ? nipsrs[ips - 1] + 1 : 1;
    int irs1 = nipsrs[ips];
    for (int irs = irs0; irs <= irs1; ++irs) {
        ++nfl;
        for (int ix = 0; ix < 4; ++ix)
            finel[k - 1][ix][nfl] = rsc[(irs - 1) * 4 + ix];
    }
}

// F90 aq_vec @1159: computes A*q product using the N**2 algorithm.
void aq_vec(const std::vector<std::vector<double>>& coord,
            const std::vector<double>& srad, int numat1,
            const std::vector<std::vector<double>>& cosurf, int nps1,
            const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
            const std::vector<int>& nset, std::vector<double>& rsc,
            std::vector<int>& nipsrs, const std::vector<int>& iatsp,
            const std::vector<double>& tm, int ioldcv, int maxrs, int lenabc,
            const std::vector<double>& a_diag, const std::vector<double>& q,
            std::vector<double>& v) {
    double finel[2][4][301] = {};
    int ijpos = 0;
    std::fill(v.begin(), v.begin() + nps1 + 1, 0.0);
    for (int ii = 1; ii <= nps1; ++ii) {
        int i = iatsp[ii];
        double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
        int nfl1 = 0;
        if (!compute_a_part) {
            double ri = srad[i];
            double xi[3] = {coord[0][i], coord[1][i], coord[2][i]};
            mfinel(ii, 1, finel, nfl1, ioldcv, maxrs, lenabc, numat1,
                   &tm[(i - 1) * 9], xi, ri);
        }
        for (int jj = 1; jj <= ii - 1; ++jj) {
            int j = iatsp[jj];
            double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
            double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                        (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                        (xb[2] - xa[2]) * (xb[2] - xa[2]);
            double aa;
            if (d2 > disex2) {
                aa = 1.0 / std::sqrt(d2);
            } else {
                if (compute_a_part) {
                    ++ijpos;
                    aa = a_part[ijpos];
                } else {
                    double xj[3] = {coord[0][j], coord[1][j], coord[2][j]};
                    double rj = srad[j];
                    int nfl2 = 0;
                    mfinel(jj, 2, finel, nfl2, ioldcv, maxrs, lenabc, numat1,
                           &tm[(j - 1) * 9], xj, rj);
                    aa = 0.0;
                    for (int k = 1; k <= nfl1; ++k) {
                        double x1 = finel[0][0][k], x2 = finel[0][1][k];
                        double x3 = finel[0][2][k], x4 = finel[0][3][k];
                        for (int l = 1; l <= nfl2; ++l) {
                            aa += x4 * finel[1][3][l] /
                                  std::sqrt((x1 - finel[1][0][l]) * (x1 - finel[1][0][l]) +
                                            (x2 - finel[1][1][l]) * (x2 - finel[1][1][l]) +
                                            (x3 - finel[1][2][l]) * (x3 - finel[1][2][l]));
                        }
                    }
                    aa /= (cosurf[3][ii] * cosurf[3][jj]);
                }
            }
            v[ii] += aa * q[jj];
            v[jj] += aa * q[ii];
        }
    }
    for (int i = 1; i <= nps1; ++i) v[i] += a_diag[i] * q[i];
}

// F90 simulate_aq_vec @1256: counts short-range elements in A matrix and,
// if requested (compute==true), fills a_part for those elements.
void simulate_aq_vec(const std::vector<std::vector<double>>& coord,
                     const std::vector<double>& srad, int numat1,
                     const std::vector<std::vector<double>>& cosurf, int nps1,
                     const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
                     const std::vector<int>& nset, std::vector<double>& rsc,
                     std::vector<int>& nipsrs, const std::vector<int>& iatsp,
                     const std::vector<double>& tm, int ioldcv, int maxrs,
                     int lenabc, bool compute, int& ndim) {
    double finel[2][4][301] = {};
    ndim = 0;
    for (int ii = 1; ii <= nps1; ++ii) {
        int i = iatsp[ii];
        double ri = srad[i];
        double xi[3] = {coord[0][i], coord[1][i], coord[2][i]};
        double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
        int nfl1 = 0;
        if (compute) {
            mfinel(ii, 1, finel, nfl1, ioldcv, maxrs, lenabc, numat1,
                   &tm[(i - 1) * 9], xi, ri);
        }
        for (int jj = 1; jj <= ii - 1; ++jj) {
            int j = iatsp[jj];
            double xj[3] = {coord[0][j], coord[1][j], coord[2][j]};
            double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
            double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                        (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                        (xb[2] - xa[2]) * (xb[2] - xa[2]);
            if (d2 <= disex2) {
                ++ndim;
                if (compute) {
                    double rj = srad[j];
                    int nfl2 = 0;
                    mfinel(jj, 2, finel, nfl2, ioldcv, maxrs, lenabc, numat1,
                           &tm[(j - 1) * 9], xj, rj);
                    double aa = 0.0;
                    for (int k = 1; k <= nfl1; ++k) {
                        double x1 = finel[0][0][k], x2 = finel[0][1][k];
                        double x3 = finel[0][2][k], x4 = finel[0][3][k];
                        for (int l = 1; l <= nfl2; ++l) {
                            aa += x4 * finel[1][3][l] /
                                  std::sqrt((x1 - finel[1][0][l]) * (x1 - finel[1][0][l]) +
                                            (x2 - finel[1][1][l]) * (x2 - finel[1][1][l]) +
                                            (x3 - finel[1][2][l]) * (x3 - finel[1][2][l]));
                        }
                    }
                    aa /= (cosurf[3][ii] * cosurf[3][jj]);
                    a_part[ndim] = aa;
                }
            }
        }
    }
}

// F90 aq_dir_int @1342: computes part of the A*q product for near surface
// elements given as direction lists (ind1/ind2). new_iteration resets ijpos.
void aq_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                const double* c, int nc, const double* q, double* r, bool same) {
    double finel[2][4][301] = {};
    static int ijpos = 0;
    double ri_dummy = c[0];  // dummy use of c (F90: ri = c(1,1))
    (void)ri_dummy; (void)nc;
    if (new_iteration) {
        ijpos = 0;
        new_iteration = false;
    }
    int maxrs = 60 * numat;
    int nfl1 = 0, nfl2 = 0;
    if (same) {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            int i = iatsp[ii];
            double ri = srad[i];
            double xi[3] = {coord[0][i], coord[1][i], coord[2][i]};
            double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
            if (!compute_a_part) {
                mfinel(ii, 1, finel, nfl1, ioldcv, maxrs, lenabc, numat,
                       &tm[(i - 1) * 9], xi, ri);
            }
            for (int j3 = 1; j3 <= i3 - 1; ++j3) {
                int jj = ind1[j3];
                int j = iatsp[jj];
                double xj[3] = {coord[0][j], coord[1][j], coord[2][j]};
                double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                double aa;
                if (d2 > disex2) {
                    aa = 1.0 / std::sqrt(d2);
                } else {
                    if (compute_a_part) {
                        ++ijpos;
                        aa = a_part[ijpos];
                    } else {
                        double rj = srad[j];
                        mfinel(jj, 2, finel, nfl2, ioldcv, maxrs, lenabc, numat,
                               &tm[(j - 1) * 9], xj, rj);
                        aa = 0.0;
                        for (int k = 1; k <= nfl1; ++k) {
                            double x1 = finel[0][0][k], x2 = finel[0][1][k];
                            double x3 = finel[0][2][k], x4 = finel[0][3][k];
                            for (int l = 1; l <= nfl2; ++l) {
                                aa += x4 * finel[1][3][l] /
                                      std::sqrt((x1 - finel[1][0][l]) * (x1 - finel[1][0][l]) +
                                                (x2 - finel[1][1][l]) * (x2 - finel[1][1][l]) +
                                                (x3 - finel[1][2][l]) * (x3 - finel[1][2][l]));
                            }
                        }
                        aa /= (cosurf[3][ii] * cosurf[3][jj]);
                    }
                }
                r[ii] += aa * q[jj];
                r[jj] += aa * q[ii];
            }
        }
    } else {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            int i = iatsp[ii];
            double ri = srad[i];
            double xi[3] = {coord[0][i], coord[1][i], coord[2][i]};
            double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
            if (!compute_a_part) {
                mfinel(ii, 1, finel, nfl1, ioldcv, maxrs, lenabc, numat,
                       &tm[(i - 1) * 9], xi, ri);
            }
            for (int j3 = 1; j3 <= n2; ++j3) {
                int jj = ind2[j3];
                int j = iatsp[jj];
                double xj[3] = {coord[0][j], coord[1][j], coord[2][j]};
                double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                double aa;
                if (d2 > disex2) {
                    aa = 1.0 / std::sqrt(d2);
                } else {
                    if (compute_a_part) {
                        ++ijpos;
                        aa = a_part[ijpos];
                    } else {
                        double rj = srad[j];
                        mfinel(jj, 2, finel, nfl2, ioldcv, maxrs, lenabc, numat,
                               &tm[(j - 1) * 9], xj, rj);
                        aa = 0.0;
                        for (int k = 1; k <= nfl1; ++k) {
                            double x1 = finel[0][0][k], x2 = finel[0][1][k];
                            double x3 = finel[0][2][k], x4 = finel[0][3][k];
                            for (int l = 1; l <= nfl2; ++l) {
                                aa += x4 * finel[1][3][l] /
                                      std::sqrt((x1 - finel[1][0][l]) * (x1 - finel[1][0][l]) +
                                                (x2 - finel[1][1][l]) * (x2 - finel[1][1][l]) +
                                                (x3 - finel[1][2][l]) * (x3 - finel[1][2][l]));
                            }
                        }
                        aa /= (cosurf[3][ii] * cosurf[3][jj]);
                    }
                }
                r[ii] += aa * q[jj];
                r[jj] += aa * q[ii];
            }
        }
    }
}

// F90 simulate_aq_dir_int @1493: counts / fills short-range A elements for
// direction lists. npos is inout (may accumulate over multiple calls).
void simulate_aq_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                         const double* c, int nc, bool same, bool compute, int& npos) {
    double finel[2][4][301] = {};
    double ri_dummy = c[0];
    (void)ri_dummy; (void)nc;
    int maxrs = 60 * numat;
    int nfl1 = 0, nfl2 = 0;
    if (same) {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            int i = iatsp[ii];
            double ri = srad[i];
            double xi[3] = {coord[0][i], coord[1][i], coord[2][i]};
            double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
            if (compute) {
                mfinel(ii, 1, finel, nfl1, ioldcv, maxrs, lenabc, numat,
                       &tm[(i - 1) * 9], xi, ri);
            }
            for (int j3 = 1; j3 <= i3 - 1; ++j3) {
                int jj = ind1[j3];
                int j = iatsp[jj];
                double xj[3] = {coord[0][j], coord[1][j], coord[2][j]};
                double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                if (d2 <= disex2) {
                    ++npos;
                    if (compute) {
                        double rj = srad[j];
                        mfinel(jj, 2, finel, nfl2, ioldcv, maxrs, lenabc, numat,
                               &tm[(j - 1) * 9], xj, rj);
                        double aa = 0.0;
                        for (int k = 1; k <= nfl1; ++k) {
                            double x1 = finel[0][0][k], x2 = finel[0][1][k];
                            double x3 = finel[0][2][k], x4 = finel[0][3][k];
                            for (int l = 1; l <= nfl2; ++l) {
                                aa += x4 * finel[1][3][l] /
                                      std::sqrt((x1 - finel[1][0][l]) * (x1 - finel[1][0][l]) +
                                                (x2 - finel[1][1][l]) * (x2 - finel[1][1][l]) +
                                                (x3 - finel[1][2][l]) * (x3 - finel[1][2][l]));
                            }
                        }
                        aa /= (cosurf[3][ii] * cosurf[3][jj]);
                        a_part[npos] = aa;
                    }
                }
            }
        }
    } else {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            int i = iatsp[ii];
            double ri = srad[i];
            double xi[3] = {coord[0][i], coord[1][i], coord[2][i]};
            double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
            if (compute) {
                mfinel(ii, 1, finel, nfl1, ioldcv, maxrs, lenabc, numat,
                       &tm[(i - 1) * 9], xi, ri);
            }
            for (int j3 = 1; j3 <= n2; ++j3) {
                int jj = ind2[j3];
                int j = iatsp[jj];
                double xj[3] = {coord[0][j], coord[1][j], coord[2][j]};
                double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                if (d2 <= disex2) {
                    ++npos;
                    if (compute) {
                        double rj = srad[j];
                        mfinel(jj, 2, finel, nfl2, ioldcv, maxrs, lenabc, numat,
                               &tm[(j - 1) * 9], xj, rj);
                        double aa = 0.0;
                        for (int k = 1; k <= nfl1; ++k) {
                            double x1 = finel[0][0][k], x2 = finel[0][1][k];
                            double x3 = finel[0][2][k], x4 = finel[0][3][k];
                            for (int l = 1; l <= nfl2; ++l) {
                                aa += x4 * finel[1][3][l] /
                                      std::sqrt((x1 - finel[1][0][l]) * (x1 - finel[1][0][l]) +
                                                (x2 - finel[1][1][l]) * (x2 - finel[1][1][l]) +
                                                (x3 - finel[1][2][l]) * (x3 - finel[1][2][l]));
                            }
                        }
                        aa /= (cosurf[3][ii] * cosurf[3][jj]);
                        a_part[npos] = aa;
                    }
                }
            }
        }
    }
}

// ============================ batch D: cgm_solve family ======================
// F90 precond_aq_dir_int @1629: fills the block-diagonal preconditioner m_vec
// (same atom pairs only) from the a_part list; npos is in/out so consecutive
// direction-list passes keep one running counter.
void precond_aq_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                        const double* c, int nc, bool same, bool compute, int& npos) {
    (void)nc; (void)compute; (void)c;   // dummy use of compute/c (F90: if (compute) aa=c(1,1))
    std::fill(m_vec.begin(), m_vec.end(), 0.0);
    if (same) {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            int i = iatsp[ii];
            double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
            for (int j3 = 1; j3 <= i3 - 1; ++j3) {
                int jj = ind1[j3];
                int j = iatsp[jj];
                double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                double aa;
                if (d2 > disex2) {
                    aa = 1.0 / std::sqrt(d2);
                } else {
                    ++npos;
                    aa = a_part[npos];
                }
                if (i == j) {   // same atom
                    int ips = ii - npoints[i] + 1;
                    int jps = jj - npoints[i] + 1;
                    if (ips >= jps)
                        m_vec[iblock_pos[i] - 1 + (ips * (ips - 1)) / 2 + jps] = aa;
                    else
                        m_vec[iblock_pos[i] - 1 + (jps * (jps - 1)) / 2 + ips] = aa;
                }
            }
        }
    } else {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            int i = iatsp[ii];
            double xa[3] = {cosurf[1][ii], cosurf[2][ii], cosurf[3][ii]};
            for (int j3 = 1; j3 <= n2; ++j3) {
                int jj = ind2[j3];
                int j = iatsp[jj];
                double xb[3] = {cosurf[1][jj], cosurf[2][jj], cosurf[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                double aa;
                if (d2 > disex2) {
                    aa = 1.0 / std::sqrt(d2);
                } else {
                    ++npos;
                    aa = a_part[npos];
                }
                if (i == j) {   // same atom
                    int ips = ii - npoints[i] + 1;
                    int jps = jj - npoints[i] + 1;
                    if (ips >= jps)
                        m_vec[iblock_pos[i] - 1 + (ips * (ips - 1)) / 2 + jps] = aa;
                    else
                        m_vec[iblock_pos[i] - 1 + (jps * (jps - 1)) / 2 + ips] = aa;
                }
            }
        }
    }
}

// F90 aq_mult_int @1724: far-surface A*q part using local expansions.
void aq_mult_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                 double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                 afmm_mod::RArr& pmn, const double* c, int n3,
                 const double* q, double* r) {
    (void)q; (void)n3;
    for (int i = 1; i <= n; ++i) {
        int ip = ind[i];
        double dx = c[(ip - 1) * n3 + 0] - x0;
        double dy = c[(ip - 1) * n3 + 1] - y0;
        double dz = c[(ip - 1) * n3 + 2] - z0;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        double cos_t = dz / d;
        double phi = std::atan2(dy, dx);
        afmm_mod::get_legendre(p, cos_t, pmn);
        double t = 0.0;
        for (int j = 0; j <= p; ++j) {
            std::complex<double> tc(0.0, 0.0);
            for (int k = 1; k <= j; ++k) {
                tc += y_norm(k, j) * pmn(k, j) * std::pow(d, j) * psi(k, j) *
                      std::exp(std::complex<double>(0.0, k * phi));
            }
            t += 2.0 * std::real(tc) + pmn(0, j) * std::pow(d, j) * std::real(psi(0, j));
        }
        r[ip] += t;
    }
}

// F90 aq_far_int @1777: far-surface A*q part using multipole expansions.
void aq_far_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                afmm_mod::RArr& pmn, const double* c, int nc,
                const double* q, double* r) {
    (void)q; (void)nc;
    for (int i = 1; i <= n; ++i) {
        int ip = ind[i];
        double dx = c[(ip - 1) * nc + 0] - x0;
        double dy = c[(ip - 1) * nc + 1] - y0;
        double dz = c[(ip - 1) * nc + 2] - z0;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        double cos_t = dz / d;
        double phi = std::atan2(dy, dx);
        afmm_mod::get_legendre(p, cos_t, pmn);
        double t = 0.0;
        for (int j = 0; j <= p; ++j) {
            std::complex<double> tc(0.0, 0.0);
            for (int k = 1; k <= j; ++k) {
                tc += y_norm(k, j) * pmn(k, j) / std::pow(d, j + 1) * psi(k, j) *
                      std::exp(std::complex<double>(0.0, k * phi));
            }
            t += 2.0 * std::real(tc) +
                 pmn(0, j) / std::pow(d, j + 1) * std::real(psi(0, j));
        }
        r[ip] += t;
    }
}

// F90 amat_diag @2420: diagonal elements of the A matrix (self-interaction of
// each surface segment; fdiag=2.1*sqrt(pi) is the diagonal singularity term).
void amat_diag(const std::vector<std::vector<double>>& coord1,
               const std::vector<double>& srad1, int numat1,
               const std::vector<std::vector<double>>& cosurf1, int nps1,
               const std::vector<int>& nar_csm1, const std::vector<int>& nsetf1,
               const std::vector<int>& nset1, std::vector<double>& rsc1,
               std::vector<int>& nipsrs1, const std::vector<int>& iatsp1,
               const std::vector<double>& tm1, int ioldcv1, int maxrs1,
               int lenabc1, std::vector<double>& m) {
    double finel[2][4][301] = {};
    double fdiag = 2.1 * std::sqrt(pi);
    for (int ii = 1; ii <= nps1; ++ii) {
        int i = iatsp1[ii];
        double ri = srad1[i];
        double xi[3] = {coord1[0][i], coord1[1][i], coord1[2][i]};
        int nfl1 = 0;
        mfinel(ii, 1, finel, nfl1, ioldcv1, maxrs1, lenabc1, numat1,
               &tm1[(i - 1) * 9], xi, ri);
        double aa = 0.0;
        for (int k = 1; k <= nfl1; ++k) {
            double f4k = finel[0][3][k];
            aa += fdiag * std::sqrt(f4k * f4k * f4k);
            double x1 = finel[0][0][k], x2 = finel[0][1][k];
            double x3 = finel[0][2][k], x4 = f4k;
            for (int l = 1; l <= k - 1; ++l) {
                aa += 2.0 * x4 * finel[0][3][l] /
                      std::sqrt((x1 - finel[0][0][l]) * (x1 - finel[0][0][l]) +
                                (x2 - finel[0][1][l]) * (x2 - finel[0][1][l]) +
                                (x3 - finel[0][2][l]) * (x3 - finel[0][2][l]));
            }
        }
        m[ii] = aa / (cosurf1[3][ii] * cosurf1[3][ii]);
    }
}

// F90 precondition @2259: builds the block-diagonal preconditioner m.
// Without a_blocks: m = 1/a_diag. With blocks: fill the same-atom blocks
// (packed lower triangle), then invert each block via dpotrf/dpotri ('U').
void precondition(const std::vector<std::vector<double>>& cosurf1, int nps1,
                  const std::vector<int>& iatsp1, int numat1,
                  const std::vector<double>& a_diag1, std::vector<double>& m) {
    if (!use_a_blocks) {
        for (int i = 1; i <= nps1; ++i) m[i] = 1.0 / a_diag1[i];
        return;
    }
    if (nps1 > na2max || nps1 > na1max) {
        int i = afmm_mod::count_short_ints(cosurf1, 4, false, precond_aq_dir_int);
        for (i = 1; i <= nps1; ++i) {
            int j = iatsp1[i];
            int ii = iblock_pos[j] - 1;
            int jj = i - npoints[j] + 1;
            m[ii + (jj * (jj + 1)) / 2] = a_diag1[i];
        }
    } else {
        int ijpos = 0;
        for (int ii = 1; ii <= nps1; ++ii) {
            int i = iatsp1[ii];
            double xa[3] = {cosurf1[1][ii], cosurf1[2][ii], cosurf1[3][ii]};
            for (int jj = 1; jj <= ii; ++jj) {
                int j = iatsp1[jj];
                double xb[3] = {cosurf1[1][jj], cosurf1[2][jj], cosurf1[3][jj]};
                double d2 = (xb[0] - xa[0]) * (xb[0] - xa[0]) +
                            (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                            (xb[2] - xa[2]) * (xb[2] - xa[2]);
                if (d2 <= disex2 && ii != jj) ++ijpos;
                if (i == j) {   // same atom
                    double aa;
                    if (ii == jj)
                        aa = a_diag1[jj];
                    else if (d2 > disex2)
                        aa = 1.0 / std::sqrt(d2);
                    else
                        aa = a_part[ijpos];
                    int ips = ii - npoints[i] + 1;
                    int jps = jj - npoints[i] + 1;
                    m[iblock_pos[i] - 1 + (ips * (ips - 1)) / 2 + jps] = aa;
                }
            }
        }
    }
    // invert the per-atom blocks (column-major upper triangle, lda = max_block_size)
    for (int i = 1; i <= numat1; ++i) {
        int ii = npoints[i + 1] - npoints[i];
        if (ii > 0) {
            int ip = iblock_pos[i];
            for (int j = 1; j <= ii; ++j)
                for (int jj = 1; jj <= j; ++jj) {
                    a_block[(jj - 1) + (j - 1) * max_block_size] = m[ip];
                    ++ip;
                }
            int ierr = 0;
            dpotrf('U', ii, a_block.data(), max_block_size, ierr);
            dpotri('U', ii, a_block.data(), max_block_size, ierr);
            if (ierr != 0) {
                std::fprintf(stderr, " Error in Matrix invert %d = %d\n", i, ierr);
                mopend(" Internal error");
                return;
            }
            ip = iblock_pos[i];
            for (int j = 1; j <= ii; ++j)
                for (int jj = 1; jj <= j; ++jj) {
                    m[ip] = a_block[(jj - 1) + (j - 1) * max_block_size];
                    ++ip;
                }
        }
    }
}

// F90 precondition_solve @2479: z = M^-1 * r. With blocks: per-atom
// mult_triangle_vec on packed (lower) triangle M-block.
void precondition_solve(const std::vector<double>& m, std::vector<double>& z,
                        const std::vector<double>& r, int n) {
    if (use_a_blocks) {
        int ipos = 1;
        for (int i = 1; i <= numat; ++i) {
            int ndim = npoints[i + 1] - npoints[i];
            if (ndim > 0) {
                mult_triangle_vec(&m[iblock_pos[i] - 1], &r[ipos - 1], ndim, &z[ipos - 1]);
                ipos += ndim;
            }
        }
    } else {
        for (int i = 1; i <= n; ++i) z[i] = r[i] * m[i];
    }
}

// F90 cgm_solve @2075: conjugate-gradient solver with block preconditioner.
void cgm_solve(std::vector<double>& x, int n, const std::vector<double>& b,
               std::vector<double>& r, std::vector<double>& p, std::vector<double>& q,
               std::vector<double>& z, std::vector<double>& m,
               bool new_x, bool new_points) {
    const int max_iters = 100;
    const double stop_tol = 1e-6, start_tol = 1e-2;
    static double current_tol = 0.0;
    int i = 0;
    new_iteration = true;
    int maxrs = 60 * numat;
    if (new_points) {
        current_tol = start_tol;
        amat_diag(coord, srad, numat, cosurf, n, nar_csm, nsetf, nset, rsc,
                  nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, a_diag);
    }
    int ierr = 0;
    if (n > na1max || n > na2max) {
        afmm_mod::set_tesselation(surface_handle, ierr);
        if (ierr != 0) { mopend("Internal error in cgm_solve"); return; }
        if (compute_a_part && new_points)
            i = afmm_mod::count_short_ints(cosurf, 4, true, simulate_aq_dir_int);
    } else {
        if (compute_a_part && new_points)
            simulate_aq_vec(coord, srad, numat, cosurf, n, nar_csm, nsetf, nset,
                            rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, true, i);
    }
    if (new_points) precondition(cosurf, n, iatsp, numat, a_diag, m);
    if (new_x) precondition_solve(m, x, b, n);
    // q = A*x
    if (n > na2max) {
        afmm_mod::afmm(cosurf, 4, n, x, q, n, aq_dir_int, aq_mult_int);
        for (i = 1; i <= n; ++i) q[i] += a_diag[i] * x[i];
    } else if (n > na1max) {
        afmm_mod::simple_mm(cosurf, 4, n, x, q, aq_dir_int, aq_far_int);
        for (i = 1; i <= n; ++i) q[i] += a_diag[i] * x[i];
    } else {
        aq_vec(coord, srad, numat, cosurf, n, nar_csm, nsetf, nset, rsc, nipsrs,
               iatsp, tm, ioldcv, maxrs, lenabc, a_diag, x, q);
    }
    for (i = 1; i <= n; ++i) r[i] = b[i] - q[i];
    double t;
    if (c_proc < 0.2) t = start_tol;
    else if (c_proc < 0.4) t = start_tol * 0.1;
    else if (c_proc < 0.6) t = start_tol * 0.01;
    else if (c_proc < 0.8) t = start_tol * 0.001;
    else t = stop_tol;
    t = std::min(t, current_tol);
    current_tol = t;
    double rho = 0.0, rho_old = 0.0, alpha = 0.0, beta = 0.0, tt = 0.0;
    for (i = 1; i <= max_iters; ++i) {
        new_iteration = true;
        precondition_solve(m, z, r, n);
        rho = ddot(n, &r[1], 1, &z[1], 1);
        if (i == 1) {
            for (int j = 1; j <= n; ++j) p[j] = z[j];
        } else {
            beta = rho / rho_old;
            for (int j = 1; j <= n; ++j) p[j] = z[j] + beta * p[j];
        }
        if (n > na2max) {
            afmm_mod::afmm(cosurf, 4, n, p, q, n, aq_dir_int, aq_mult_int);
            for (int j = 1; j <= n; ++j) q[j] += a_diag[j] * p[j];
        } else if (n > na1max) {
            afmm_mod::simple_mm(cosurf, 4, n, p, q, aq_dir_int, aq_far_int);
            for (int j = 1; j <= n; ++j) q[j] += a_diag[j] * p[j];
        } else {
            aq_vec(coord, srad, numat, cosurf, n, nar_csm, nsetf, nset, rsc,
                   nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, a_diag, p, q);
        }
        alpha = rho / ddot(n, &p[1], 1, &q[1], 1);
        for (int j = 1; j <= n; ++j) { x[j] += alpha * p[j]; r[j] -= alpha * q[j]; }
        tt = some_norm(&r[0], n);   // some_norm is 1-based
        if (tt < t) break;
        rho_old = rho;
    }
    if (tt > t)
        std::fprintf(stderr, " Warning! cgm_solve finished with error = %g\n", tt);
}


// F90 bz_far_int @1828: far-surface B*q part using multipole expansions.
// ind lists ATOMS; each atom contributes its surface points
// npoints(ia) .. npoints(ia+1)-1.
void bz_far_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                afmm_mod::RArr& pmn, const double* c, int nc,
                const double* q, double* r) {
    (void)q; (void)c; (void)nc;
    for (int i = 1; i <= n; ++i) {
        int ia = ind[i];
        for (int jj = npoints[ia]; jj <= npoints[ia + 1] - 1; ++jj) {
            double dx = cosurf[1][jj] - x0;
            double dy = cosurf[2][jj] - y0;
            double dz = cosurf[3][jj] - z0;
            double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            double cos_t = dz / d;
            double phi = std::atan2(dy, dx);
            afmm_mod::get_legendre(p, cos_t, pmn);
            double t = 0.0;
            for (int j = 0; j <= p; ++j) {
                std::complex<double> tc(0.0, 0.0);
                for (int k = 1; k <= j; ++k) {
                    tc += y_norm(k, j) * pmn(k, j) / std::pow(d, j + 1) * psi(k, j) *
                          std::exp(std::complex<double>(0.0, k * phi));
                }
                t += 2.0 * std::real(tc) +
                     pmn(0, j) / std::pow(d, j + 1) * std::real(psi(0, j));
            }
            r[jj] += t;
        }
    }
}

// F90 bz_mult_int @1885: near-surface B*q part using local expansions.
void bz_mult_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                 double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                 afmm_mod::RArr& pmn, const double* c, int nc,
                 const double* q, double* r) {
    (void)q; (void)c; (void)nc;
    for (int i = 1; i <= n; ++i) {
        int ia = ind[i];
        for (int jj = npoints[ia]; jj <= npoints[ia + 1] - 1; ++jj) {
            double dx = cosurf[1][jj] - x0;
            double dy = cosurf[2][jj] - y0;
            double dz = cosurf[3][jj] - z0;
            double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            double cos_t = dz / d;
            double phi = std::atan2(dy, dx);
            afmm_mod::get_legendre(p, cos_t, pmn);
            double t = 0.0;
            for (int j = 0; j <= p; ++j) {
                std::complex<double> tc(0.0, 0.0);
                for (int k = 1; k <= j; ++k) {
                    tc += y_norm(k, j) * pmn(k, j) * std::pow(d, j) * psi(k, j) *
                          std::exp(std::complex<double>(0.0, k * phi));
                }
                t += 2.0 * std::real(tc) +
                     pmn(0, j) * std::pow(d, j) * std::real(psi(0, j));
            }
            r[jj] += t;
        }
    }
}

// F90 bp_dir_int @1951: near-surface B*q part (density/charge interaction).
// w = get_bvec multipoles; over the packed triangle jj=ijbo_diag(j):
// off-diagonal w contribute -2*p*w, diagonal -p*w; plus nuclear tore*w(1).
void bp_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                const double* c, int nc, const double* q, double* r, bool same) {
    double w[46];   // F90 w(45), 1-based
    (void)q; (void)nc;
    if (same) {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            for (int i = npoints[ii]; i <= npoints[ii + 1] - 1; ++i) {
                double t = 0.0;
                double xa[3] = {cosurf[1][i], cosurf[2][i], cosurf[3][i]};
                for (int j3 = 1; j3 <= n1; ++j3) {
                    int j = ind1[j3];
                    int nao = nlast[j] - nfirst[j];
                    get_bvec(xa, &c[(j - 1) * nc], nao, nat[j], w);
                    int jj = ijbo_diag[j];
                    int m = 0;
                    for (int k = 1; k <= nao + 1; ++k) {
                        for (int l = 1; l <= k - 1; ++l) {
                            ++m;
                            t -= 2 * p[jj + m] * w[m];
                        }
                        ++m;
                        t -= p[jj + m] * w[m];
                    }
                    t += tore[nat[j]] * w[1];
                }
                r[i] += t;
            }
        }
    } else {
        // --- (1) points in box 1 with atoms in box 2 ---
        for (int i3 = 1; i3 <= n1; ++i3) {
            int ii = ind1[i3];
            for (int i = npoints[ii]; i <= npoints[ii + 1] - 1; ++i) {
                double t = 0.0;
                double xa[3] = {cosurf[1][i], cosurf[2][i], cosurf[3][i]};
                for (int j3 = 1; j3 <= n2; ++j3) {
                    int j = ind2[j3];
                    int nao = nlast[j] - nfirst[j];
                    get_bvec(xa, &c[(j - 1) * nc], nao, nat[j], w);
                    int jj = ijbo_diag[j];
                    int m = 0;
                    for (int k = 1; k <= nao + 1; ++k) {
                        for (int l = 1; l <= k - 1; ++l) {
                            ++m;
                            t -= 2 * p[jj + m] * w[m];
                        }
                        ++m;
                        t -= p[jj + m] * w[m];
                    }
                    t += tore[nat[j]] * w[1];
                }
                r[i] += t;
            }
        }
        // --- (2) points in box 2 with atoms in box 1 ---
        for (int i3 = 1; i3 <= n2; ++i3) {
            int ii = ind2[i3];
            for (int i = npoints[ii]; i <= npoints[ii + 1] - 1; ++i) {
                double t = 0.0;
                double xa[3] = {cosurf[1][i], cosurf[2][i], cosurf[3][i]};
                for (int j3 = 1; j3 <= n1; ++j3) {
                    int j = ind1[j3];
                    int nao = nlast[j] - nfirst[j];
                    get_bvec(xa, &c[(j - 1) * nc], nao, nat[j], w);
                    int jj = ijbo_diag[j];
                    int m = 0;
                    for (int k = 1; k <= nao + 1; ++k) {
                        for (int l = 1; l <= k - 1; ++l) {
                            ++m;
                            t -= 2 * p[jj + m] * w[m];
                        }
                        ++m;
                        t -= p[jj + m] * w[m];
                    }
                    t += tore[nat[j]] * w[1];
                }
                r[i] += t;
            }
        }
    }
}


// F90 addnucz @2544: zero the channel-1 phi/qscnet/qdenet, then load the
// nuclear charges tore(nat(i)) into qdenet(idenat(i),1).
void addnucz(std::vector<std::vector<double>>& phinet,
             std::vector<std::vector<double>>& qscnet,
             std::vector<std::vector<double>>& qdenet) {
    for (int i = 1; i <= nps; ++i) { phinet[i][1] = 0.0; qscnet[i][1] = 0.0; }
    for (int i = 1; i <= lm61; ++i) qdenet[i][1] = 0.0;
    for (int i = 1; i <= numat; ++i) qdenet[idenat[i]][1] = tore[nat[i]];
}

// F90 addfckz @2573: add the surface-charge correction to the Fock matrix and
// the solvation energy.  N^2 branch (numat <= nb1max): bz_vec+bpnew_vec build
// the potential, cgm_solve inverts A to get the induced charges, then the
// density interaction updates f and solv_energy.
void addfckz() {
    const double fcon = a0 * ev;
    for (int i = 1; i <= numat; ++i) qscat[i] = 0.0;
    for (int i = 1; i <= nps; ++i) phinet[i][2] = 0.0;
    if (numat > nb1max || numat > nb2max) {
        // large system: net charges q_vec(i) = -trace-density + nuclear
        static const int idiag[10] = {0, 1, 3, 6, 10, 15, 21, 28, 36, 45};
        for (int i = 1; i <= numat; ++i) {
            int ii = ijbo_diag[i];
            int j = nlast[i] - nfirst[i] + 1;
            double s1 = 0.0;
            for (int k = 1; k <= j; ++k) s1 += p[ii + idiag[k]];
            q_vec[i] = -s1 + tore[nat[i]];
        }
        int ierr = 0;
        afmm_mod::set_tesselation(atom_handle, ierr);
        if (ierr == -1000) return;   // dummy use of ierr
    }

    std::vector<double> phi2(nps + 1, 0.0);
    if (numat > nb2max) {
        afmm_mod::afmm(coord, 3, numat, q_vec, phi2, nps, bp_dir_int, bz_mult_int);
    } else if (numat > nb1max) {
        afmm_mod::simple_mm(coord, 3, numat, q_vec, phi2, bp_dir_int, bz_far_int);
    } else {
        std::vector<double> v1(nps + 1, 0.0), v2(nps + 1, 0.0);
        bz_vec(v1);
        bpnew_vec(v2);
        for (int i = 1; i <= nps; ++i) {
            phi2[i] = v1[i] + v2[i];
            phinet[i][1] = 0.0;
        }
    }
    for (int i = 1; i <= nps; ++i) {
        phinet[i][2] = phi2[i];
        phinet[i][3] = phi2[i];
    }

    std::vector<double> qs2(nps + 1, 0.0);
    for (int i = 1; i <= nps; ++i) qs2[i] = qscnet[i][2];
    if (new_surface) {
        cgm_solve(qs2, nps, phi2, r_vec, p_vec, q_vec, z_vec, m_vec, true, true);
        new_surface = false;
    } else {
        for (int i = 1; i <= nps; ++i) qs2[i] /= (-fepsi);
        cgm_solve(qs2, nps, phi2, r_vec, p_vec, q_vec, z_vec, m_vec, false, false);
    }
    for (int i = 1; i <= nps; ++i) qscnet[i][2] = qs2[i];

    ediel = 0.0;
    double s1 = 0.0;
    for (int i = 1; i <= nps; ++i) {
        int iat = iatsp[i];
        qscnet[i][2] = -fepsi * qscnet[i][2];
        double qsc3 = qscnet[i][1] + qscnet[i][2];
        qscnet[i][3] = qsc3;
        ediel += qsc3 * phinet[i][3];
        s1 += qscnet[i][1];
        qscat[iat] += qsc3;
    }
    ediel *= fcon / 2;

    if (numat > nb1max || numat > nb2max) {
        int ierr = 0;
        afmm_mod::set_tesselation(atom_handle, ierr);
        for (int i = 1; i <= lm61; ++i) qdenet[i][2] = 0.0;
    }

    if (numat > nb2max) {
        std::vector<double> qd2(lm61 + 1, 0.0);
        afmm_mod::afmm(coord, 3, numat, qs2, qd2, lm61, fock_dir_int,
                       fock_mult_int, sphere_f_multipoles);
        for (int i = 1; i <= lm61; ++i) {
            qdenet[i][2] = qd2[i];
            f[ipiden[i]] -= qdenet[i][2] * fcon;
        }
        s1 = 0.0;
        for (int i = 1; i <= lm61; ++i) s1 += qdenet[i][2] * p[ipiden[i]] * gden[i];
        s1 = -s1 * fcon;
    } else if (numat > nb1max) {
        std::vector<double> qd2(lm61 + 1, 0.0);
        afmm_mod::simple_mm(coord, 3, numat, qs2, qd2, fock_dir_int,
                            fock_far_int, sphere_f_multipoles);
        for (int i = 1; i <= lm61; ++i) {
            qdenet[i][2] = qd2[i];
            f[ipiden[i]] -= qdenet[i][2] * fcon;
        }
        s1 = 0.0;
        for (int i = 1; i <= lm61; ++i) s1 += qdenet[i][2] * p[ipiden[i]] * gden[i];
        s1 = -s1 * fcon;
    } else {
        s1 = 0.0;
        int im = 0;
        for (int i = 1; i <= numat; ++i) {
            int nao = nlast[i] - nfirst[i];
            int ii = ijbo_diag[i];
            if (nao == 0) {   // --- S element
                double v1 = 0.0;
                for (int j = 1; j <= nps; ++j) {
                    double dx = cosurf[1][j] - coord[0][i];
                    double dy = cosurf[2][j] - coord[1][i];
                    double dz = cosurf[3][j] - coord[2][i];
                    v1 += qscnet[j][2] / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
                v1 = -v1 * fcon;
                f[ii + 1] += v1;
                ++im;
                s1 += v1 * p[ii + 1] * gden[im];
            } else {
                int nnn = ((nao + 2) * (nao + 1)) / 2;
                std::vector<double> v(nnn + 1, 0.0), w(46, 0.0);
                for (int j = 1; j <= nps; ++j) {
                    double xa[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
                    double xj[3] = {coord[0][i], coord[1][i], coord[2][i]};
                    get_bvec(xa, xj, nao, nat[i], w.data());
                    for (int k = 1; k <= nnn; ++k) v[k] += w[k] * qscnet[j][2];
                }
                for (int k = 1; k <= nnn; ++k) v[k] = -v[k] * fcon;
                for (int j = 1; j <= nnn; ++j) {
                    f[ii + j] += v[j];
                    ++im;
                    s1 += v[j] * p[ii + j] * gden[im];
                }
            }
        }
    }

    solv_energy = s1 * 0.5 + ediel;

    // --- qdenet from the density matrix ---
    for (int i = 1; i <= lm61; ++i) {
        qdenet[i][2] = gden[i] * p[ipiden[i]];
        qdenet[i][3] = qdenet[i][2] + qdenet[i][1];
    }
}

// F90 fock_dir_int @2747: near-charge correction to the Fock matrix.
// For each atom: v(1:nnn) += get_bvec multipoles * qscnet(point,2), then
// r(iatom_pos(i)+j) += v(j).  same -> atoms/points in one box (n1 x n1);
// !same -> bidirectional two-box interaction.
void fock_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                  const double* c, int nc, const double* q, double* r, bool same) {
    double w[46];   // F90 w(45), 1-based
    (void)q; (void)nc;
    if (same) {
        for (int i3 = 1; i3 <= n1; ++i3) {
            int i = ind1[i3];
            int nao = nlast[i] - nfirst[i];
            int nnn = ((nao + 2) * (nao + 1)) / 2;
            double v[46] = {0.0};
            int ii = iatom_pos[i];
            for (int j3 = 1; j3 <= n1; ++j3) {
                int jj = ind1[j3];
                for (int j = npoints[jj]; j <= npoints[jj + 1] - 1; ++j) {
                    double xa[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
                    get_bvec(xa, &c[(i - 1) * nc], nao, nat[i], w);
                    for (int k = 1; k <= nnn; ++k) v[k] += w[k] * qscnet[j][2];
                }
            }
            for (int j = 1; j <= nnn; ++j) r[ii + j] += v[j];
        }
    } else {
        // --- (1) atoms in box 1 with points in box 2 ---
        for (int i3 = 1; i3 <= n1; ++i3) {
            int i = ind1[i3];
            int nao = nlast[i] - nfirst[i];
            int nnn = ((nao + 2) * (nao + 1)) / 2;
            double v[46] = {0.0};
            int ii = iatom_pos[i];
            for (int j3 = 1; j3 <= n2; ++j3) {
                int jj = ind2[j3];
                for (int j = npoints[jj]; j <= npoints[jj + 1] - 1; ++j) {
                    double xa[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
                    get_bvec(xa, &c[(i - 1) * nc], nao, nat[i], w);
                    for (int k = 1; k <= nnn; ++k) v[k] += w[k] * qscnet[j][2];
                }
            }
            for (int j = 1; j <= nnn; ++j) r[ii + j] += v[j];
        }
        // --- (2) atoms in box 2 with points in box 1 ---
        for (int i3 = 1; i3 <= n2; ++i3) {
            int i = ind2[i3];
            int nao = nlast[i] - nfirst[i];
            int nnn = ((nao + 2) * (nao + 1)) / 2;
            double v[46] = {0.0};
            int ii = iatom_pos[i];
            for (int j3 = 1; j3 <= n1; ++j3) {
                int jj = ind1[j3];
                for (int j = npoints[jj]; j <= npoints[jj + 1] - 1; ++j) {
                    double xa[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
                    get_bvec(xa, &c[(i - 1) * nc], nao, nat[i], w);
                    for (int k = 1; k <= nnn; ++k) v[k] += w[k] * qscnet[j][2];
                }
            }
            for (int j = 1; j <= nnn; ++j) r[ii + j] += v[j];
        }
    }
}

// F90 fock_mult_int @2850: far-charge Fock correction via local expansions
// (d^j multipoles).  r(iatom_pos(i)+idiag(1:nnn)) += t.
void fock_mult_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                   double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                   afmm_mod::RArr& pmn, const double* c, int nc,
                   const double* q, double* r) {
    static const int idiag[10] = {0, 1, 3, 6, 10, 15, 21, 28, 36, 45};
    (void)q; (void)nc;
    for (int i3 = 1; i3 <= n; ++i3) {
        int i = ind[i3];
        int ii = iatom_pos[i];
        int nnn = nlast[i] - nfirst[i] + 1;
        double dx = c[(i - 1) * nc] - x0;
        double dy = c[(i - 1) * nc + 1] - y0;
        double dz = c[(i - 1) * nc + 2] - z0;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        double cos_t = dz / d;
        double phi = std::atan2(dy, dx);
        afmm_mod::get_legendre(p, cos_t, pmn);
        double t = 0.0;
        for (int j = 0; j <= p; ++j) {
            std::complex<double> tc(0.0, 0.0);
            for (int k = 1; k <= j; ++k) {
                tc += y_norm(k, j) * pmn(k, j) * std::pow(d, j) * psi(k, j) *
                      std::exp(std::complex<double>(0.0, k * phi));
            }
            t += 2.0 * std::real(tc) +
                 pmn(0, j) * std::pow(d, j) * std::real(psi(0, j));
        }
        for (int j = 1; j <= nnn; ++j) r[ii + idiag[j]] += t;
    }
}

// F90 fock_far_int @2923: far-charge Fock correction via multipole expansions
// (1/d^(j+1) multipoles).  Same accumulator pattern as fock_mult_int.
void fock_far_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                  double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                  afmm_mod::RArr& pmn, const double* c, int nc,
                  const double* q, double* r) {
    static const int idiag[10] = {0, 1, 3, 6, 10, 15, 21, 28, 36, 45};
    (void)q; (void)nc;
    for (int i3 = 1; i3 <= n; ++i3) {
        int i = ind[i3];
        int ii = iatom_pos[i];
        int nnn = nlast[i] - nfirst[i] + 1;
        double dx = c[(i - 1) * nc] - x0;
        double dy = c[(i - 1) * nc + 1] - y0;
        double dz = c[(i - 1) * nc + 2] - z0;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        double cos_t = dz / d;
        double phi = std::atan2(dy, dx);
        afmm_mod::get_legendre(p, cos_t, pmn);
        double t = 0.0;
        for (int j = 0; j <= p; ++j) {
            std::complex<double> tc(0.0, 0.0);
            for (int k = 1; k <= j; ++k) {
                tc += y_norm(k, j) * pmn(k, j) / std::pow(d, j + 1) * psi(k, j) *
                      std::exp(std::complex<double>(0.0, k * phi));
            }
            t += 2.0 * std::real(tc) +
                 pmn(0, j) / std::pow(d, j + 1) * std::real(psi(0, j));
        }
        for (int j = 1; j <= nnn; ++j) r[ii + idiag[j]] += t;
    }
}

// F90 sphere_f_multipoles @2997: multipole expansion of surface charges
// (TrickyFn).  Accumulates psi over each surface point of each box atom.
void sphere_f_multipoles(double x0, double y0, double z0, const int* ind, int n,
                         afmm_mod::RArr& pmn, const afmm_mod::RArr& y_norm,
                         const double* c, int nc, const double* q,
                         afmm_mod::CArr& psi, int p) {
    (void)c; (void)nc; (void)q;
    for (int j = 1; j <= n; ++j) {
        int j1 = ind[j];
        for (int k = npoints[j1]; k <= npoints[j1 + 1] - 1; ++k) {
            double dx = cosurf[1][k] - x0;
            double dy = cosurf[2][k] - y0;
            double dz = cosurf[3][k] - z0;
            double qq = qscnet[k][2];
            double rr = std::sqrt(dx * dx + dy * dy + dz * dz);
            double cos_t = dz / rr;
            double phi = std::atan2(dy, dx);
            afmm_mod::get_legendre(p, cos_t, pmn);
            psi(0, 0) += std::complex<double>(qq, 0.0);
            double t = 1.0;
            for (int n1 = 1; n1 <= p; ++n1) {
                t *= rr;
                for (int m1 = 0; m1 <= n1; ++m1) {
                    std::complex<double> tc =
                        std::exp(std::complex<double>(0.0, -m1 * phi));
                    double tt = y_norm(-m1, n1) * pmn(-m1, n1) * t * qq;
                    psi(m1, n1) += tc * tt;
                }
            }
            for (int n1 = 1; n1 <= p; ++n1)
                for (int m1 = 1; m1 <= n1; ++m1)
                    psi(-m1, n1) = std::conj(psi(m1, n1));
        }
    }
}

// F90 am1dft_solve @3073: wrapper for cgm_solve used by writmn (AM1DFT).
void am1dft_solve(std::vector<double>& qsct, std::vector<double>& phit, int nps) {
    cgm_solve(qsct, nps, phit, r_vec, p_vec, q_vec, z_vec, m_vec, false, false);
}

}  // namespace linear_cosmo
