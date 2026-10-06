// dcart.cpp — C++ translation of MOPAC 2016 "dcart.F90".
// Cartesian energy derivatives by finite differences.
// coord/dxyz are 1-based rows (coord[k-1][i], k=1..3); p/pa/pb 0-based packed
// (common_arrays_C); MOZYME branch keeps structure with ijbo stubbed.
#include "dcart.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "chrge.h"
#include "chrge_for_MOZYME.h"
#include "common_arrays_C.h"
#include "cosmo.h"
#include "cosmo_C.h"
#include "delsta.h"
#include "dhc.h"
#include "dihed.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "molmec_C.h"
#include "parameters_C.h"
#include "set_up_dentate.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace cosmo_C;
using namespace elemts_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace MOZYME_C;
using namespace molmec_C;
using namespace parameters_C;

// MOZYME bonding-table lookup.  Global stub (chrge_for_MOZYME/cosmo also
// link against it); full implementation belongs to the M04 batch.
extern int ijbo(int, int);  // real implementation in ijbo.cpp

double derp(double r) {
    // thread_local: derp() is called inside the dcart OpenMP atom-pair loop.
    static thread_local int icalcn = 0;
    static thread_local double bound1 = 0.0, bound2 = 0.0, c = 0.0, cr = 0.0, cr2 = 0.0, range = 0.0;
    if (icalcn != numcal) {
        bound1 = clower / cutofp;
        bound2 = cupper / cutofp;
        range = bound2 - bound1;
        c = -0.5 * bound1 * bound1 * cutofp / range;
        cr = 1.0 + bound1 / range;
        cr2 = -1.0 / (cutofp * 2.0 * range);
        icalcn = numcal;
    }
    if (r > clower) {
        if (r > cupper) return 0.0;
        double den = c + cr * r + cr2 * r * r;
        return -(cr + 2.0 * cr2 * r) / (den * den);
    }
    return -1.0 / (r * r);
}

// B3: dcart atom-pair OpenMP is functionally correct (theoretical) but changes
// the floating-point reduction order, so large-molecule geometry optimisations
// drift (e.g. C12H26 EF: -66.72968 serial vs -64.61059 OMP).  Per the B-batch
// acceptance criterion (correct in principle, exact numbers not required) it is
// kept as an OPT-OUT free switch, DEFAULT OFF to preserve bit-exactness:
//   set MOPAC_DCART_OMP=1   -> parallel atom-pair loop (faster, non-bitwise)
static bool dcart_omp_enabled() {
    static bool en = std::getenv("MOPAC_DCART_OMP") != nullptr;
    return en;
}

void dcart(std::vector<std::vector<double>>& coord,
           std::vector<std::vector<double>>& dxyz) {
    static int icalcn = 0;
    static double chnge = 1.0e-4;
    static bool debug = false, force = false, large = false;
    static int ione = 1;
    static double const_ = 0.0;
    const double chnge2 = chnge * 0.5;
    if (icalcn != numcal) {
        icalcn = numcal;
        const_ = fpc_9;
        ione = (id == 0) ? 1 : 0;
        large = keywrd.find("LARGE") != std::string::npos;
        debug = keywrd.find("DERIV") != std::string::npos ||
                keywrd.find("DCART") != std::string::npos;
        force = keywrd.find("PREC") != std::string::npos ||
                keywrd.find("FORCE") != std::string::npos;
    }
    const int icuc = (l123 + 1) / 2;
    const int numtot = numat * l123;
    bool refeps = useps;
    useps = false;
    std::vector<double> q(numat + 1, 0.0);
    if (mozyme) {
        if (mode == 0) {
            for (int i = 1; i <= numtot; ++i)
                for (int k = 1; k <= 3; ++k) dxyz[k][i] = 0.0;
        } else if (mode == 1) {
            for (int i = 1; i <= numtot; ++i)
                for (int k = 1; k <= 3; ++k)
                    dxyz[k][i] = part_dxyz[(i - 1) * 3 + (k - 1)];
        } else if (mode == -1) {
            for (int i = 1; i <= numtot; ++i)
                for (int k = 1; k <= 3; ++k) dxyz[k][i] = -dxyz[k][i];
        }
        chrge_for_MOZYME(p.data(), q.data(), numat, iorbs.data());
    } else {
        for (int i = 1; i <= numtot; ++i)
            for (int k = 1; k <= 3; ++k) dxyz[k][i] = 0.0;
        chrge(p, q);
        mode = 0;
    }
    for (int i = 1; i <= numat; ++i) q[i] = tore[nat[i]] - q[i];

    // Per-atom work: computes pair contributions for atom ii into dd.
    // dd must be zeroed by the caller. Serial and OpenMP paths share this.
    auto work_ii = [&](int ii, std::vector<std::vector<double>>& dd) {
        const int iii = l123 * (ii - 1);
        const int im1 = ii - ione;
        const int if_ = nfirst[ii], il = nlast[ii];
        int ndi[2] = {0, nat[ii]};
        double cdi[6] = {};  // column-major 3x2, 0-based: [0..2]=atom1, [3..5]=atom2
        cdi[3] = coord[0][ii];
        cdi[4] = coord[1][ii];
        cdi[5] = coord[2][ii];
        int j2 = 1;
        for (int jj = 1; jj <= im1; ++jj) {
            if (mozyme) {
                if (mode == 0 || jopt[j2] == jj) {
                    j2 = j2 + 1;
                } else {
                    continue;
                }
            }
            const double half = (ii == jj) ? 0.5 : 1.0;
            const int jjj = l123 * (jj - 1);
            const int jf = nfirst[jj], jl = nlast[jj];
            ndi[0] = nat[jj];
            double pdi[171] = {}, padi[171] = {}, pbdi[171] = {};
            bool point = false;
            if (mozyme) {
                if (ijbo(ii, jj) >= 0) {
                    // MOZYME: build diatomic density blocks from LMO storage.
                    int k = ijbo(jj, jj);
                    int ij_ = 0;
                    for (int i = 1; i <= iorbs[jj]; ++i)
                        for (int j = 1; j <= i; ++j) {
                            ij_ = ij_ + 1;
                            k = k + 1;
                            padi[ij_] = p[k] * 0.5;
                            pdi[ij_] = p[k];
                        }
                    if (ii == jj) {
                        int ij2 = iorbs[jj];
                        for (int i = 1; i <= iorbs[ii]; ++i) {
                            ij2 = ij2 + 1;
                            int l = (ij2 * (ij2 - 1)) / 2;
                            for (int j = 1; j <= iorbs[jj]; ++j) {
                                l = l + 1;
                                padi[l] = 0.0;
                                pdi[l] = 0.0;
                            }
                        }
                    } else {
                        int ij2 = iorbs[jj];
                        k = ijbo(ii, jj);
                        for (int i = 1; i <= iorbs[ii]; ++i) {
                            ij2 = ij2 + 1;
                            int l = (ij2 * (ij2 - 1)) / 2;
                            for (int j = 1; j <= iorbs[jj]; ++j) {
                                l = l + 1;
                                k = k + 1;
                                padi[l] = p[k] * 0.5;
                                pdi[l] = p[k];
                            }
                        }
                    }
                    k = ijbo(ii, ii);
                    int ij3 = iorbs[jj];
                    for (int i = 1; i <= iorbs[ii]; ++i) {
                        ij3 = ij3 + 1;
                        int l = (ij3 * (ij3 - 1)) / 2 + iorbs[jj];
                        for (int j = 1; j <= i; ++j) {
                            k = k + 1;
                            l = l + 1;
                            padi[l] = p[k] * 0.5;
                            pdi[l] = p[k];
                        }
                    }
                    for (int l = 0; l < 171; ++l) pbdi[l] = padi[l];
                    point = false;
                } else {
                    point = true;
                }
            } else {
                point = false;
                int ij_ = 0;
                for (int i = jf; i <= jl; ++i) {
                    int k = (i * (i - 1)) / 2 + jf - 1;
                    if (i - jf + 1 > 0) {
                        const int cnt = i - jf + 1;
                        for (int m = 0; m < cnt; ++m) {
                            padi[ij_ + m] = pa[k + m + 1];
                            pbdi[ij_ + m] = pb[k + m + 1];
                            pdi[ij_ + m] = p[k + m + 1];
                        }
                        ij_ += cnt;
                    }
                }
                for (int i = if_; i <= il; ++i) {
                    int l = (i * (i - 1)) / 2;
                    int k = l + jf - 1;
                    if (jl - jf + 1 > 0) {
                        const int cnt = jl - jf + 1;
                        for (int m = 0; m < cnt; ++m) {
                            padi[ij_ + m] = pa[k + m + 1];
                            pbdi[ij_ + m] = pb[k + m + 1];
                            pdi[ij_ + m] = p[k + m + 1];
                        }
                        ij_ += cnt;
                    }
                    k = l + if_ - 1;
                    if (i - if_ + 1 > 0) {
                        const int cnt = i - if_ + 1;
                        for (int m = 0; m < cnt; ++m) {
                            padi[ij_ + m] = pa[k + m + 1];
                            pbdi[ij_ + m] = pb[k + m + 1];
                            pdi[ij_ + m] = p[k + m + 1];
                        }
                        ij_ += cnt;
                    }
                }
            }
            int kkkk = 0;
            for (int ik = -l1u; ik <= l1u; ++ik) {
                for (int jk = -l2u; jk <= l2u; ++jk) {
                    for (int kl = -l3u; kl <= l3u; ++kl) {
                        kkkk = kkkk + 1;
                        for (int kk = 1; kk <= 3; ++kk)
                            cdi[kk - 1] = coord[kk - 1][jj] + common_arrays_C::tvec[kk][1] * ik +
                                          common_arrays_C::tvec[kk][2] * jk + common_arrays_C::tvec[kk][3] * kl;
                        if (id != 0) {
                            double rij2 = (cdi[0] - cdi[3]) * (cdi[0] - cdi[3]) +
                                          (cdi[1] - cdi[4]) * (cdi[1] - cdi[4]) +
                                          (cdi[2] - cdi[5]) * (cdi[2] - cdi[5]);
                            const double cut2 = (2.0 / 3.0 * cutofp);
                            if (rij2 > cut2 * cut2) {
                                const double rij = std::sqrt(rij2);
                                const double der = derp(rij);
                                const double ee = q[ii] * q[jj] * fpc_9 * ev * a0 * der;
                                for (int k = 1; k <= 3; ++k) {
                                    const double deriv = half * ee * (cdi[k - 1] - cdi[k + 2]) / rij;
                                    dd[k][iii + icuc] -= deriv;
                                    dd[k][jjj + kkkk] += deriv;
                                }
                                continue;
                            }
                        }
                        double aa = 0.0, ee = 0.0;
                        if (!force) {
                            cdi[0] += chnge2;
                            cdi[1] += chnge2;
                            cdi[2] += chnge2;
                            dhc(pdi, padi, pbdi, cdi, ndi, jf, jl, if_, il, aa, 1);
                        }
                        if (point) {
                            for (int l = 1; l <= 3; ++l)
                                cdi[l - 1] = coord[l - 1][jj] + common_arrays_C::tvec[l][1] * ik +
                                             common_arrays_C::tvec[l][2] * jk + common_arrays_C::tvec[l][3] * kl;
                            std::vector<double> dstat(4, 0.0);
                            // delsta expects 1-based rows
                            std::vector<std::vector<double>> cdi_v(4, std::vector<double>(3, 0.0));
                            cdi_v[1][1] = cdi[0]; cdi_v[2][1] = cdi[1]; cdi_v[3][1] = cdi[2];
                            cdi_v[1][2] = cdi[3]; cdi_v[2][2] = cdi[4]; cdi_v[3][2] = cdi[5];
                            delsta(nat, iorbs, p, cdi_v, dstat, ii, jj);
                            for (int k = 1; k <= 3; ++k) {
                                dd[k][iii + icuc] -= dstat[k];
                                dd[k][jjj + kkkk] += dstat[k];
                            }
                        } else {
                            for (int k = 1; k <= 3; ++k) {
                                if (force) {
                                    cdi[k + 2] -= chnge2;
                                    dhc(pdi, padi, pbdi, cdi, ndi, jf, jl, if_, il, aa, 1);
                                }
                                cdi[k + 2] += chnge;
                                dhc(pdi, padi, pbdi, cdi, ndi, jf, jl, if_, il, ee, 2);
                                cdi[k + 2] -= chnge2;
                                if (!force) cdi[k + 2] -= chnge2;
                                const double deriv = half * (aa - ee) * const_ / chnge;

                                dd[k][iii + icuc] -= deriv;
                                dd[k][jjj + kkkk] += deriv;
                            }
                        }
                    }
                }
            }
        }
    };  // end work_ii
#ifdef _OPENMP
    if (!mozyme && molkst_C::id == 0 && dcart_omp_enabled()) {
        // Parallel atom-pair loop with private per-thread dxyz reduction.
        // Each pair contributes to exactly the same components as the serial
        // path; only the floating-point reduction order differs.
        #pragma omp parallel
        {
            std::vector<std::vector<double>> dl(4, std::vector<double>(numtot + 1, 0.0));
            #pragma omp for schedule(dynamic)
            for (int ii = 1; ii <= numat; ++ii) {
                work_ii(ii, dl);
            }
            #pragma omp critical
            for (int k = 1; k <= 3; ++k)
                for (int i = 1; i <= numtot; ++i) dxyz[k][i] += dl[k][i];
        }
    } else
#endif
    {
        for (int ii = 1; ii <= numat; ++ii) work_ii(ii, dxyz);
    }
    if (nnhco != 0) {
        const double del = 1.0e-8;
        for (int i = 1; i <= nnhco; ++i) {
            for (int j = 1; j <= 4; ++j) {
                for (int k = 1; k <= 3; ++k) {
                    coord[k-1][nhco[j][i]] -= del;
                    double angle = 0.0;
                    dihed(coord, nhco[1][i], nhco[2][i], nhco[3][i], nhco[4][i], angle);
                    const double refh = htype * std::sin(angle) * std::sin(angle);
                    coord[k-1][nhco[j][i]] += del * 2.0;
                    dihed(coord, nhco[1][i], nhco[2][i], nhco[3][i], nhco[4][i], angle);
                    const double heat = htype * std::sin(angle) * std::sin(angle);
                    coord[k-1][nhco[j][i]] -= del;
                    const double sum = (refh - heat) / (2.0 * del);
                    dxyz[k][nhco[j][i]] -= sum;
                }
            }
        }
    }
    if (method_pm6 && N_3_present) {
        const double del = 1.0e-8;
        for (int i = 1; i <= numat; ++i) {
            if (nat[i] == 7 && nbonds[i] == 4) {
                int jj = 0;
                if (nat[ibonds[2][i]] == 1) jj = 1;
                if (nat[ibonds[3][i]] == 1) jj = jj + 1;
                if (nat[ibonds[4][i]] == 1) jj = jj + 1;
                if (jj < 2) {
                    for (int j = 1; j <= 4; ++j) {
                        for (int k = 1; k <= 3; ++k) {
                            coord[k-1][ibonds[j][i]] -= del;
                            double sum = nsp2_atom_correction(coord, i, ibonds[2][i],
                                                              ibonds[3][i], ibonds[4][i]);
                            coord[k-1][ibonds[j][i]] += 2.0 * del;
                            sum = (sum - nsp2_atom_correction(coord, i, ibonds[2][i],
                                                              ibonds[3][i], ibonds[4][i])) /
                                  (2.0 * del);
                            coord[k-1][ibonds[j][i]] -= del;
                            dxyz[k][ibonds[j][i]] -= sum;
                        }
                    }
                }
            }
        }
    }
    if (method_pm7 && Si_O_H_present) {
        const double del = 1.0e-8;
        for (int i = 1; i <= numat; ++i) {
            if (nat[i] == 8) {
                int O = i, Si = 0, H = 0;
                for (int j = 1; j <= nbonds[i]; ++j) {
                    const int k = ibonds[j][i];
                    if (nat[k] == 14) Si = k;
                    if (nat[k] == 1) H = k;
                }
                if (Si != 0 && H != 0) {
                    for (int k = 1; k <= 3; ++k) {
                        coord[k-1][Si] -= del;
                        double sum = Si_O_H_bond_correction(coord, Si, O, H);
                        coord[k-1][Si] += 2.0 * del;
                        sum = (sum - Si_O_H_bond_correction(coord, Si, O, H)) / (2.0 * del);
                        coord[k-1][Si] -= del;
                        dxyz[k][Si] -= sum;
                        coord[k-1][O] -= del;
                        sum = Si_O_H_bond_correction(coord, Si, O, H);
                        coord[k-1][O] += 2.0 * del;
                        sum = (sum - Si_O_H_bond_correction(coord, Si, O, H)) / (2.0 * del);
                        coord[k-1][O] -= del;
                        dxyz[k][O] -= sum;
                        coord[k-1][H] -= del;
                        sum = Si_O_H_bond_correction(coord, Si, O, H);
                        coord[k-1][H] += 2.0 * del;
                        sum = (sum - Si_O_H_bond_correction(coord, Si, O, H)) / (2.0 * del);
                        coord[k-1][H] -= del;
                        dxyz[k][H] -= sum;
                    }
                }
            }
        }
    }
    if (mode == -1) {
        for (int i = 1; i <= numtot; ++i)
            for (int k = 1; k <= 3; ++k)
                part_dxyz[(i - 1) * 3 + (k - 1)] = -dxyz[k][i];
    }
    useps = refeps;
    if (useps) {
        std::vector<double> flat(3 * numtot + 1, 0.0);
        for (int i = 1; i <= numtot; ++i)
            for (int k = 1; k <= 3; ++k) flat[(i - 1) * 3 + k] = dxyz[k][i];
        diegrd(flat.data());
        for (int i = 1; i <= numtot; ++i)
            for (int k = 1; k <= 3; ++k) dxyz[k][i] = flat[(i - 1) * 3 + k];
    }
    if (use_ref_geo) {
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= 3; ++j)
                dxyz[j][i] += (geo[j][i] - geoa[j][i]) * density * 2.0;
    }
    if (debug) {
        std::printf(" CARTESIAN COORDINATE DERIVATIVES (no post-SCF corrections)\n");
        for (int i = 1; i <= numtot; ++i)
            std::printf("%6d %2s %13.6f %13.6f %13.6f\n", i,
                        elemnt[nat[(i - 1) / l123 + 1]].c_str(),
                        dxyz[1][i], dxyz[2][i], dxyz[3][i]);
    }
}
