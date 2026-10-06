// esp.cpp — C++ translation of MOPAC 2016 "esp.F90".
// ESP (electrostatic potential) point-charge fitting: Williams grid (pdgrid),
// Connolly surface (surfac), nuclear+electronic potential (elesn), least-
// squares fit (espfit).  Column-major Fortran arrays are kept 1-based:
//   a(j,k) -> a[j][k]  (outer index = first Fortran dimension).
#include "esp.h"

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "densit.h"
#include "dist2.h"
#include "elemts_C.h"
#include "esp_C.h"
#include "funcon_C.h"
#include "genun.h"
#include "gmetry.h"
#include "molkst_C.h"
#include "mopend.h"
#include "mult.h"
#include "overlaps_C.h"
#include "parameters_C.h"
#include "reada.h"
#include "rsp.h"
#include "second.h"
#include "setupg.h"
#include "to_screen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>

using esp_C::nesp;
using esp_C::idip;
using esp_C::iz;
using esp_C::ipx;
using esp_C::isc;
using esp_C::is_esp;
using esp_C::icd;
using esp_C::ipe;
using esp_C::npr;
using esp_C::ic;
using esp_C::ip;
using esp_C::ncc;
using esp_C::dens;
using esp_C::scale;
using esp_C::cf;
using esp_C::rms;
using esp_C::rrms;
using esp_C::dx;
using esp_C::dy;
using esp_C::dz;
using esp_C::den;
using esp_C::fv;
using esp_C::dex;
using esp_C::tf;
using esp_C::fac;
using esp_C::ind;
using esp_C::itemp;
using esp_C::ird;
using esp_C::indc;
using esp_C::iam;
using esp_C::cc;
using esp_C::ex;
using esp_C::temp;
using esp_C::rad;
using esp_C::es;
using esp_C::b_esp;
using esp_C::esp_array;
using esp_C::cesp;
using esp_C::cespml;
using esp_C::al;
using esp_C::qesp;
using esp_C::td;
using esp_C::dx_array;
using esp_C::dy_array;
using esp_C::dz_array;
using esp_C::ptd;
using esp_C::pexs;
using esp_C::pce;
using esp_C::pexpn;
using esp_C::pewcx;
using esp_C::pewcy;
using esp_C::pewcz;
using esp_C::pf0;
using esp_C::pf1;
using esp_C::pf2;
using esp_C::fc;
using esp_C::qsc;
using esp_C::cequiv;
using esp_C::cespm;
using esp_C::cen;
using esp_C::co;
using esp_C::potpt;
using esp_C::ovl;
using esp_C::cespm2;
using esp_C::a2;
using esp_C::exs;
using esp_C::ce;
using esp_C::expn;
using esp_C::ewcx;
using esp_C::ewcy;
using esp_C::ewcz;
using esp_C::f0;
using esp_C::f1;
using esp_C::u_esp;
using esp_C::rnai;
using esp_C::rnai1;
using esp_C::rnai2;
using esp_C::espi;
using esp_C::exsr;

using common_arrays_C::geo;
using common_arrays_C::coord;
using common_arrays_C::labels;
using common_arrays_C::nat;
using common_arrays_C::q;
using common_arrays_C::c;
using molkst_C::keywrd;
using molkst_C::moperr;
using molkst_C::numat;
using molkst_C::natoms;
using molkst_C::norbs;
using molkst_C::nopen;
using molkst_C::fract;
using molkst_C::nclose;
using molkst_C::ux;
using molkst_C::uy;
using molkst_C::uz;
using molkst_C::method_mndo;
using chanel_C::iw;
using chanel_C::esp_fn;
using chanel_C::esr_fn;
using chanel_C::restart_fn;
using funcon_C::a0;
using funcon_C::fpc_1;
using funcon_C::fpc_8;
using parameters_C::tore;
using parameters_C::zs;
using parameters_C::zp;
using overlaps_C::allc;
using overlaps_C::allz;
using elemts_C::elemnt;

namespace {

// Partial-pivoting Gaussian elimination in place, solving A*x = b.
// A is n x n (a[row][col]); b is overwritten with the solution.
// LAPACK dgesv replacement (espfit.F90 dependency; no LAPACK in this project).
void dgesv_local(std::vector<std::vector<double>>& a, int n,
                 std::vector<double>& b) {
    for (int k = 1; k <= n; ++k) {
        int p = k;
        double mx = std::fabs(a[k][k]);
        for (int i = k + 1; i <= n; ++i) {
            if (std::fabs(a[i][k]) > mx) {
                mx = std::fabs(a[i][k]);
                p = i;
            }
        }
        if (p != k) {
            for (int j = 1; j <= n; ++j) std::swap(a[p][j], a[k][j]);
            std::swap(b[p], b[k]);
        }
        for (int i = k + 1; i <= n; ++i) {
            double f = a[i][k] / a[k][k];
            a[i][k] = f;
            for (int j = k + 1; j <= n; ++j) a[i][j] -= f * a[k][j];
        }
    }
    // Forward substitution (L y = b), then back substitution (U x = y).
    for (int i = 2; i <= n; ++i)
        for (int k = 1; k < i; ++k) b[i] -= a[i][k] * b[k];
    for (int i = n; i >= 1; --i) {
        for (int k = i + 1; k <= n; ++k) b[i] -= a[i][k] * b[k];
        b[i] /= a[i][i];
    }
}
}  // namespace

// ---------------------------------------------------------------------
// esp() — main driver.
// ---------------------------------------------------------------------
void esp() {
    int n, i;
    double scincr, time1;
    time1 = second(1);
    if (keywrd.find(" CUBE") != std::string::npos ||
        keywrd.find(" ESPGRID") != std::string::npos) {
        // new_esp()  — CUBE/ESPGRID branch handled by new_esp.cpp.
        return;
    }
    std::string::size_type pos;
    if ((pos = keywrd.find("SCALE=")) != std::string::npos)
        scale = reada(keywrd, (int)pos);
    else
        scale = 1.4;
    if ((pos = keywrd.find("DEN=")) != std::string::npos)
        den = reada(keywrd, (int)pos);
    else
        den = 1.0;
    if ((pos = keywrd.find("SCINCR=")) != std::string::npos)
        scincr = reada(keywrd, (int)pos);
    else
        scincr = 0.20;
    if ((pos = keywrd.find("NSURF=")) != std::string::npos)
        n = (int)std::lround(reada(keywrd, (int)pos));
    else
        n = 4;
    time1 = second(1);
    setup_esp(1);
    nesp = 0;
    if (keywrd.find("WILLIAMS") != std::string::npos) {
        pdgrid();
        if (moperr) return;
    } else {
        for (i = 1; i <= n; ++i) {
            surfac();
            if (moperr) return;
            scale = scale + scincr;
        }
    }
    potcal();
    time1 = second(1) - time1;
    // write(iw,20) 'TIME TO CALCULATE ESP:', time1, ' SECONDS'
    return;
}

// ---------------------------------------------------------------------
// pdgrid() — Williams surface on a cubic grid.
// ---------------------------------------------------------------------
void pdgrid() {
    int icntr, i, j, ia, npnt, l, jz;
    double vderw[54] = {0.0};
    double dist[101] = {0.0};
    double xmin[4] = {0.0}, xmax[4] = {0.0};
    double shell, grid, closer, vdmax, xstart, ystart, zstart, zgrid, ygrid,
        xgrid;
    vderw[1] = 2.4;
    vderw[5] = 3.0;
    vderw[6] = 2.9;
    vderw[7] = 2.7;
    vderw[8] = 2.6;
    vderw[9] = 2.55;
    vderw[15] = 3.1;
    vderw[16] = 3.05;
    vderw[17] = 3.0;
    vderw[35] = 3.15;
    vderw[53] = 3.35;
    shell = 1.2;
    nesp = 0;
    grid = 0.8;
    closer = 0.0;
    // CONVERT INTERNAL TO CARTESIAN COORDINATES
    gmetry(geo, coord);
    // STRIP COORDINATES AND ATOM LABEL FOR DUMMIES (I.E. 99)
    icntr = 0;
    if (co.size() < (size_t)natoms + 1)
        co.assign(4, std::vector<double>(natoms + 1, 0.0));
    for (i = 1; i <= natoms; ++i) {
        co[1][i] = coord[0][i];
        co[2][i] = coord[1][i];
        co[3][i] = coord[2][i];
        if (labels[i] == 99) continue;
        ++icntr;
        nat[icntr] = labels[i];
    }
    numat = icntr;
    for (i = 1; i <= numat; ++i) {
        j = nat[i];
        if (vderw[j] == 0.0) {
            mopend("VAN DER WAALS' RADIUS NOT DEFINED FOR ATOM IN WILLIAMS "
                   "SURFACE ROUTINE PDGRID!");
            return;
        }
    }
    // NOW CREATE LIMITS FOR A BOX
    for (int k = 1; k <= 3; ++k) {
        xmin[k] = 100000.0;
        xmax[k] = -100000.0;
    }
    for (ia = 1; ia <= numat; ++ia) {
        for (int k = 1; k <= 3; ++k) {
            if (co[k][ia] - xmin[k] < 0.0) xmin[k] = co[k][ia];
            if (co[k][ia] - xmax[k] > 0.0) xmax[k] = co[k][ia];
        }
    }
    // ADD (OR SUBTRACT) THE MAXIMUM VDERW PLUS SHELL
    vdmax = 0.0;
    for (i = 1; i <= 53; ++i)
        if (vderw[i] > vdmax) vdmax = vderw[i];
    for (int k = 1; k <= 3; ++k) {
        xmin[k] -= vdmax + shell;
        xmax[k] += vdmax + shell;
    }
    // STEP GRID BACK FROM ZERO TO FIND STARTING POINTS
    xstart = -grid;
    while (xstart > xmin[1]) xstart -= grid;
    ystart = -grid;
    while (ystart > xmin[2]) ystart -= grid;
    zstart = -grid;
    while (zstart > xmin[3]) zstart -= grid;
    npnt = 0;
    zgrid = zstart;
    for (;;) {
        ygrid = ystart;
        for (;;) {
            xgrid = xstart;
            for (;;) {
                for (l = 1; l <= numat; ++l) {
                    jz = nat[l];
                    dist[l] = std::sqrt((co[1][l] - xgrid) * (co[1][l] - xgrid) +
                                        (co[2][l] - ygrid) * (co[2][l] - ygrid) +
                                        (co[3][l] - zgrid) * (co[3][l] - zgrid));
                    // REJECT GRID POINT IF ANY ATOM IS TOO CLOSE
                    if (dist[l] < vderw[jz] - closer) goto g220;
                }
                // BUT AT LEAST ONE ATOM MUST BE CLOSE ENOUGH
                for (l = 1; l <= numat; ++l) {
                    jz = nat[l];
                    if (dist[l] > vderw[jz] + shell) continue;
                    goto g210;
                }
                goto g220;
            g210:
                ++npnt;
                ++nesp;
                potpt[1][nesp] = xgrid;
                potpt[2][nesp] = ygrid;
                potpt[3][nesp] = zgrid;
            g220:;
                xgrid += grid;
                if (xgrid <= xmax[1]) continue;
                break;
            }
            ygrid += grid;
            if (ygrid <= xmax[2]) continue;
            break;
        }
        zgrid += grid;
        if (zgrid <= xmax[3]) continue;
        break;
    }
}

// ---------------------------------------------------------------------
// surfac() — Connolly solvent-accessible surface.
// ---------------------------------------------------------------------
void surfac() {
    std::vector<int> ias(numat + 1, 0);
    int icntr, i, ipoint, iatom, nnbr, jatom, ncon;
    double vander[101] = {0.0};
    double pi, rw, ri, d2;
    bool si;
    const double vander_init[100] = {
        1.20, 1.20, 1.37, 1.45, 1.45, 1.50, 1.50, 1.40, 1.35, 1.30,
        1.57, 1.36, 1.24, 1.17, 1.80, 1.75, 1.70, 0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
        0,    0,    0,    0,    2.30, 0,    0,    0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0};
    for (i = 1; i <= 100; ++i) vander[i] = vander_init[i - 1];
    // INSERT VAN DER WAAL RADII FOR ZINC
    vander[30] = 1.00;
    pi = 4.0 * std::atan(1.0);
    // CONVERT INTERNAL TO CARTESIAN COORDINATES
    gmetry(geo, coord);
    // STRIP COORDINATES AND ATOM LABEL FOR DUMMIES (I.E. 99)
    icntr = 0;
    for (i = 1; i <= numat; ++i) {
        co[1][i] = coord[0][i];
        co[2][i] = coord[1][i];
        co[3][i] = coord[2][i];
    }
    for (i = 1; i <= natoms; ++i) {
        if (labels[i] == 99) continue;
        ++icntr;
        nat[icntr] = labels[i];
    }
    // ONLY VAN DER WAALS' TYPE SURFACE IS GENERATED
    rw = 0.0;
    numat = icntr;
    dens = den;
    for (i = 1; i <= numat; ++i) {
        ipoint = nat[i];
        rad[i] = vander[ipoint] * scale;
        ias[i] = 2;
    }
    // BIG LOOP FOR EACH ATOM
    for (iatom = 1; iatom <= numat; ++iatom) {
        if (ias[iatom] == 0) continue;
        ri = rad[iatom];
        si = ias[iatom] == 2;
        std::array<double, 3> ci = {co[1][iatom], co[2][iatom], co[3][iatom]};
        // GATHER THE NEIGHBORING ATOMS OF IATOM
        nnbr = 0;
        std::vector<std::vector<double>> cnbr(4, std::vector<double>(201, 0.0));
        std::vector<double> rnbr(201, 0.0);
        for (jatom = 1; jatom <= numat; ++jatom) {
            if (iatom == jatom || ias[jatom] == 0) continue;
            std::array<double, 3> cj = {co[1][jatom], co[2][jatom], co[3][jatom]};
            d2 = dist2(ci, cj);
            if (d2 >= (2.0 * rw + ri + rad[jatom]) * (2.0 * rw + ri + rad[jatom]))
                continue;
            ++nnbr;
            if (nnbr > 200) {
                mopend("ERROR.  TOO MANY NEIGHBORS");
                return;
            }
            cnbr[1][nnbr] = co[1][jatom];
            cnbr[2][nnbr] = co[2][jatom];
            cnbr[3][nnbr] = co[3][jatom];
            rnbr[nnbr] = rad[jatom];
        }
        // CONTACT SURFACE
        if (!si) continue;
        ncon = (int)((4.0 * pi * ri * ri) * den);
        ncon = std::min(1000, ncon);
        if (ncon == 0) {
            mopend("VECTOR LENGTH OF ZERO IN SURFAC");
            return;
        }
        std::vector<std::vector<double>> con(4, std::vector<double>(1001, 0.0));
        genun(con, ncon);
        // CONTACT PROBE PLACEMENT LOOP
        for (i = 1; i <= ncon; ++i) {
            std::vector<double> cw(4, 0.0);
            cw[1] = ci[0] + (ri + rw) * con[1][i];
            cw[2] = ci[1] + (ri + rw) * con[2][i];
            cw[3] = ci[2] + (ri + rw) * con[3][i];
            if (collid(rw, cw, cnbr, rnbr, nnbr, 1)) continue;
            ++nesp;
            potpt[1][nesp] = ci[0] + ri * con[1][i];
            potpt[2][nesp] = ci[1] + ri * con[2][i];
            potpt[3][nesp] = ci[2] + ri * con[3][i];
        }
    }
}

// ---------------------------------------------------------------------
// elesn() — electronic contribution to the ESP (STO-6G auxiliary basis).
// ---------------------------------------------------------------------
void elesn() {
    int j1, i, m, k, j, nqn, ipc, in, l;
    double norm, pi, t, fval, sum, ra, rij;
    bool potwrt, sto3g;
    pi = 4.0 * std::atan(1.0);
    for (i = -1; i <= 10; ++i) dex[i + 1] = dex2(i);
    if (fac[4] > 1.0) {
        for (i = 1; i <= 7; ++i) fac[i] = 1.0 / fac[i];
    }
    for (m = 0; m <= 8; ++m) {
        k = 1;
        fv[m][1] = 1.0 / (2.0 * (double)m + 1.0);
        t = 0.0;
        for (j1 = 1; j1 <= (int)(41.0 / 0.05); ++j1) {
            t += 0.05;
            ++k;
            fsub(m, t, fval);
            fv[m][k] = fval;
        }
        t += 0.05;
    }
    // LOAD BASIS FUNCTIONS INTO ARRAYS
    sto3g = keywrd.find("STO3G") != std::string::npos;
    if (sto3g) {
        icd = 3;
        setup3();
    } else {
        icd = 6;
        setupg();
    }
    ncc = 0;
    npr = 0;
    for (i = 1; i <= numat; ++i) {
        if (nat[i] <= 2) {
            for (k = 1; k <= icd; ++k) {
                cc[npr + k] = allc[k][1][1];
                ex[npr + k] = allz[k][1][1] * zs[1] * zs[1];
                cen[npr + k][1] = co[1][i] / a0;
                cen[npr + k][2] = co[2][i] / a0;
                cen[npr + k][3] = co[3][i] / a0;
                iam[npr + k][1] = 0;
                iam[npr + k][2] = 0;
                fc[npr + k] = (double)i;
            }
            ++ncc;
            npr += icd;
        } else {
            // DETERMINE PRINCIPAL QUANTUM NUMBER(NQN) OF ORBITALS TO BE USED
            nqn = 2;
            if (nat[i] > 10 && nat[i] <= 18) nqn = 3;
            if (nat[i] > 18 && nat[i] <= 36) nqn = 4;
            if (nat[i] > 36 && nat[i] <= 54) nqn = 5;
            for (k = 1; k <= icd; ++k) {
                cc[npr + k] = allc[k][nqn][1];
                ex[npr + k] = allz[k][nqn][1] * zs[nat[i]] * zs[nat[i]];
                cen[npr + k][1] = co[1][i] / a0;
                cen[npr + k][2] = co[2][i] / a0;
                cen[npr + k][3] = co[3][i] / a0;
                iam[npr + k][1] = 0;
                iam[npr + k][2] = 0;
            }
            ++ncc;
            npr += icd;
            for (k = 1; k <= 3; ++k) {
                for (j = 1; j <= icd; ++j) {
                    cc[npr + j] = allc[j][nqn][2];
                    ex[npr + j] = allz[j][nqn][2] * zp[nat[i]] * zp[nat[i]];
                    cen[npr + j][1] = co[1][i] / a0;
                    cen[npr + j][2] = co[2][i] / a0;
                    cen[npr + j][3] = co[3][i] / a0;
                    iam[npr + j][1] = 1;
                    iam[npr + j][2] = k;
                }
                ++ncc;
                npr += icd;
            }
        }
    }
    // CALCULATE NORMALIZATION CONSTANTS AND INCLUDE THEM IN THE CONTRACTION
    // COEFFICIENTS
    for (i = 1; i <= npr; ++i) {
        norm = std::pow(2.0 * ex[i] / pi, 0.75) *
               std::pow(4.0 * ex[i], iam[i][1] / 2) /
               std::sqrt(dex[2 * iam[i][1] - 1 + 1]);
        cc[i] *= norm;
    }
    // PERFORM SORT OF PRIMITIVES BY ANGULAR MOMENTUM
    is_esp = 0;
    isc = 0;
    for (i = 1; i <= npr; ++i) {
        if (iam[i][1] != 0) continue;
        ++is_esp;
        ind[is_esp] = i;
    }
    ip = is_esp;
    for (i = 1; i <= npr; ++i) {
        if (iam[i][1] != 1 || iam[i][2] != 1) continue;
        ++ip;
        ind[ip] = i;
    }
    for (i = 1; i <= npr; ++i) {
        if (iam[i][1] != 1 || iam[i][2] != 2) continue;
        ++ip;
        ind[ip] = i;
    }
    for (i = 1; i <= npr; ++i) {
        if (iam[i][1] != 1 || iam[i][2] != 3) continue;
        ++ip;
        ind[ip] = i;
    }
    for (i = 1; i <= ncc; ++i) {
        in = i * icd - icd + 1;
        if (iam[in][1] != 0) continue;
        ++isc;
        indc[isc] = i;
    }
    ipc = isc;
    for (i = 1; i <= ncc; ++i) {
        in = i * icd - icd + 1;
        if (iam[in][1] != 1 || iam[in][2] != 1) continue;
        ++ipc;
        indc[ipc] = i;
    }
    for (i = 1; i <= ncc; ++i) {
        in = i * icd - icd + 1;
        if (iam[in][1] != 1 || iam[in][2] != 2) continue;
        ++ipc;
        indc[ipc] = i;
    }
    for (i = 1; i <= ncc; ++i) {
        in = i * icd - icd + 1;
        if (iam[in][1] != 1 || iam[in][2] != 3) continue;
        ++ipc;
        indc[ipc] = i;
    }
    for (i = 1; i <= npr; ++i) temp[i] = cc[ind[i]];
    for (i = 1; i <= npr; ++i) cc[i] = temp[i];
    for (i = 1; i <= npr; ++i) temp[i] = ex[ind[i]];
    for (i = 1; i <= npr; ++i) ex[i] = temp[i];
    for (i = 1; i <= npr; ++i) temp[i] = cen[ind[i]][1];
    for (i = 1; i <= npr; ++i) cen[i][1] = temp[i];
    for (i = 1; i <= npr; ++i) temp[i] = cen[ind[i]][2];
    for (i = 1; i <= npr; ++i) cen[i][2] = temp[i];
    for (i = 1; i <= npr; ++i) temp[i] = cen[ind[i]][3];
    for (i = 1; i <= npr; ++i) cen[i][3] = temp[i];
    for (i = 1; i <= npr; ++i) itemp[i] = iam[ind[i]][1];
    for (i = 1; i <= npr; ++i) iam[i][1] = itemp[i];
    for (i = 1; i <= npr; ++i) itemp[i] = iam[ind[i]][2];
    for (i = 1; i <= npr; ++i) iam[i][2] = itemp[i];
    // CALCULATE OVERLAP MATRIX OF STO-6G FUNCTIONS
    for (j = 1; j <= ncc; ++j) ovlp(j);
    for (j = 1; j <= ncc; ++j)
        for (k = 1; k <= ncc; ++k) cespm2[indc[j]][indc[k]] = ovl[j][k];
    for (j = 1; j <= ncc; ++j)
        for (k = 1; k <= ncc; ++k) ovl[j][k] = cespm2[j][k];
    l = 0;
    for (i = 1; i <= ncc; ++i) {
        for (k = 1; k <= i; ++k) {
            cesp[l + k - 1] = ovl[i][k];  // 0-based packed
        }
        l += i;
    }
    // DEORTHOGONALIZE THE COEFFICIENTS AND REFORM THE DENSITY MATRIX
    rsp(cesp.data(), ncc, temp.data(), cespml.data());
    for (i = 1; i <= ncc; ++i) {
        for (j = 1; j <= i; ++j) {
            sum = 0.0;
            for (k = 1; k <= ncc; ++k) {
                sum += cespml[(k - 1) * ncc + (i - 1)] / std::sqrt(temp[k - 1]) *
                       cespml[(k - 1) * ncc + (j - 1)];
            }
            cesp[(j - 1) * ncc + (i - 1)] = sum;
            cesp[(i - 1) * ncc + (j - 1)] = sum;
        }
    }
    // mult(c, cesp, cespml, ncc)
    {
        std::vector<double> cflat((size_t)ncc * ncc + 1, 0.0);
        for (i = 1; i <= ncc; ++i)
            for (j = 1; j <= ncc; ++j)
                cflat[(j - 1) * ncc + (i - 1)] = c[i][j];
        mult(cflat.data(), cesp.data(), cespml.data(), ncc);
    }
    // densit(cespml, ncc, ncc, nclose, 2.d0, nopen, fract, cesp, 2)
    // densit's p argument is 1-based packed; cesp here is 0-based, so copy
    // through a temporary 1-based buffer.
    {
        std::vector<std::vector<double>> cm(ncc + 1,
                                            std::vector<double>(ncc + 1, 0.0));
        for (i = 1; i <= ncc; ++i)
            for (j = 1; j <= ncc; ++j)
                cm[i][j] = cespml[(j - 1) * ncc + (i - 1)];
        std::vector<double> p1(ncc * (ncc + 1) / 2 + 1, 0.0);
        densit(cm, ncc, ncc, nclose, 2.0, nopen, fract, p1, 2);
        for (i = 1; i <= ncc; ++i)
            for (j = 1; j <= i; ++j) {
                int l1 = i * (i - 1) / 2 + j;  // 1-based packed index
                cesp[l1 - 1] = p1[l1];
            }
    }
    // NOW CALCULATE THE ELECTRONIC CONTRIBUTION TO THE ESP
    l = 0;
    for (i = 1; i <= ncc; ++i) {
        for (j = 1; j <= i; ++j) {
            ++l;
            cespm[i][j] = cesp[l - 1];  // cesp is 0-based packed
            cespm[j][i] = cesp[l - 1];
        }
    }
    ipx = (npr - is_esp) / 3;
    ipe = is_esp + ipx;
    for (i = 1; i <= nesp; ++i) es[i] = 0.0;
    naicas();
    naicap();
    // CALCULATE TOTAL ESP AND FORM ARRAYS FOR ESPFIT
    for (i = 1; i <= numat + 4; ++i) b_esp[i] = 0.0;
    for (i = 1; i <= nesp; ++i) {
        esp_array[i] = 0.0;
        for (j = 1; j <= numat; ++j) {
            ra = std::sqrt((co[1][j] - potpt[1][i]) * (co[1][j] - potpt[1][i]) +
                           (co[2][j] - potpt[2][i]) * (co[2][j] - potpt[2][i]) +
                           (co[3][j] - potpt[3][i]) * (co[3][j] - potpt[3][i]));
            esp_array[i] += tore[nat[j]] / (ra / a0);
        }
        esp_array[i] -= es[i];
        for (j = 1; j <= numat; ++j) {
            rij = std::sqrt((co[1][j] - potpt[1][i]) * (co[1][j] - potpt[1][i]) +
                            (co[2][j] - potpt[2][i]) * (co[2][j] - potpt[2][i]) +
                            (co[3][j] - potpt[3][i]) * (co[3][j] - potpt[3][i])) /
                  a0;
            b_esp[j] += esp_array[i] * 1.0 / rij;
        }
    }
    // IF REQUESTED WRITE OUT ELECTRIC POTENTIAL DATA TO UNIT IESP
    potwrt = keywrd.find("POTWRT") != std::string::npos;
    if (potwrt) {
        // open(iesp, file=esp_fn, status='UNKNOWN', position='asis')
        // write(iesp,'(I5)') nesp
        // do i=1,nesp: write(iesp,420) esp_array(i),potpt(1,i)/a0,...
        // close(unit=iesp)
    }
}

// ---------------------------------------------------------------------
// espfit() — least-squares fit of ESP to point charges (dgesv).
// ---------------------------------------------------------------------
void espfit() {
    int k, j, i;
    double rik, rij, espc, atmp;
    int istart, iend;
    cf = fpc_1 * a0 * fpc_8 * 1.0e-10;
    std::vector<std::vector<double>> a(numat + 5,
                                       std::vector<double>(numat + 5, 0.0));
    // THE FOLLOWING SETS UP THE LINEAR EQUATION A*Q=B
    // SET UP THE A(J,K) ARRAY
    istart = 1;
    iend = nesp;
    for (k = 1; k <= numat; ++k) {
        for (j = 1; j <= numat; ++j) {
            atmp = 0.0;
            for (i = istart; i <= iend; ++i) {
                rik = (co[1][k] - potpt[1][i]) * (co[1][k] - potpt[1][i]) +
                      (co[2][k] - potpt[2][i]) * (co[2][k] - potpt[2][i]) +
                      (co[3][k] - potpt[3][i]) * (co[3][k] - potpt[3][i]);
                rij = (co[1][j] - potpt[1][i]) * (co[1][j] - potpt[1][i]) +
                      (co[2][j] - potpt[2][i]) * (co[2][j] - potpt[2][i]) +
                      (co[3][j] - potpt[3][i]) * (co[3][j] - potpt[3][i]);
                atmp += 1.0 / std::sqrt(rij * rik);
            }
            a[j][k] = atmp * a0 * a0;
        }
    }
    for (k = 1; k <= numat; ++k) a[numat + 1][k] = 1.0;  // row
    for (j = numat + 1; j <= numat + 4; ++j)             // corner
        for (k = numat + 1; k <= numat + 4; ++k) a[j][k] = 0.0;
    for (j = 1; j <= numat; ++j) a[j][numat + 1] = 1.0;  // column
    if (idip == 1) {
        for (j = 1; j <= numat; ++j) {
            a[numat + 2][j] = co[1][j] / a0;
            a[numat + 3][j] = co[2][j] / a0;
            a[numat + 4][j] = co[3][j] / a0;
            a[j][numat + 2] = co[1][j] / a0;
            a[j][numat + 3] = co[2][j] / a0;
            a[j][numat + 4] = co[3][j] / a0;
        }
        for (j = numat + 2; j <= numat + 4; ++j)
            for (k = numat + 2; k <= numat + 4; ++k) a[j][k] = 0.0;
    } else {
        for (j = 1; j <= numat + 4; ++j)
            for (k = numat + 2; k <= numat + 4; ++k) a[j][k] = 0.0;
        for (j = numat + 2; j <= numat + 4; ++j)
            for (k = 1; k <= numat + 1; ++k) a[j][k] = 0.0;
    }
    // INSERT CHARGE AND DIPOLAR (IF DESIRED) CONSTRAINTS
    b_esp[numat + 1] = (double)iz;
    b_esp[numat + 2] = dx / cf;
    b_esp[numat + 3] = dy / cf;
    b_esp[numat + 4] = dz / cf;
    q.assign(numat + 5, 0.0);
    for (i = 1; i <= numat + 4; ++i) q[i] = b_esp[i];
    if (idip == 1)
        dgesv_local(a, numat + 4, q);
    else
        dgesv_local(a, numat + 1, q);
    // CALCULATE ROOT MEAN SQUARE FITS AND RELATIVE ROOT MEAN SQUARE FITS
    rms = 0.0;
    rrms = 0.0;
    for (i = 1; i <= nesp; ++i) {
        espc = 0.0;
        for (j = 1; j <= numat; ++j) {
            rij = std::sqrt((co[1][j] - potpt[1][i]) * (co[1][j] - potpt[1][i]) +
                            (co[2][j] - potpt[2][i]) * (co[2][j] - potpt[2][i]) +
                            (co[3][j] - potpt[3][i]) * (co[3][j] - potpt[3][i])) /
                  a0;
            espc += q[j] / rij;
        }
        rms += (espc - esp_array[i]) * (espc - esp_array[i]);
        rrms += esp_array[i] * esp_array[i];
    }
    rms = std::sqrt(rms / nesp);
    rrms = rms / std::sqrt(rrms / nesp);
    rms *= 627.51;
}

// ---------------------------------------------------------------------
// naicas() — (S|S) and (S|P) ESP integrals over ESP points.
// ---------------------------------------------------------------------
void naicas() {
    int ipx2, i, jstart, np, ic, ipr, istart, j, iesp, iref, k, ips, jps;
    double pi, potp1, potp2, potp3, ref, res, term, f, ts, term1, fi, ts1;
    pi = 4.0 * std::atan(1.0);
    ipx2 = 2 * ipx;
    // IF THIS IS A RESTART RUN, READ IN RESTART INFO
    jstart = 1;
    np = is_esp + 1;
    for (ic = jstart; ic <= isc; ++ic) {
        ipr = ic * icd - icd + 1;
        istart = ipr;
        for (i = istart; i <= ipe; ++i) {
            dx_array[i] = cen[ipr][1] - cen[i][1];
            dy_array[i] = cen[ipr][2] - cen[i][2];
            dz_array[i] = cen[ipr][3] - cen[i][3];
            td[i] = dx_array[i] * dx_array[i] + dy_array[i] * dy_array[i] +
                    dz_array[i] * dz_array[i];
        }
        // CALCULATE EXPONENT SUM
        for (i = istart; i <= ipe; ++i) {
            for (j = 1; j <= icd; ++j) {
                exsr[i][j] = ex[ipr + j - 1] + ex[i];
                exs[i][j] = 1.0 / exsr[i][j];
                ce[i][j] = ex[ipr + j - 1] * ex[i] * exs[i][j];
                expn[i][j] = std::exp(-ce[i][j] * td[i]);
            }
        }
        // CALCULATE EXPONENT WEIGHTED CENTERS
        for (i = istart; i <= ipe; ++i) {
            for (j = 1; j <= icd; ++j) {
                ewcx[i][j] = (ex[i] * cen[i][1] +
                              ex[ipr + j - 1] * cen[ipr + j - 1][1]) *
                             exs[i][j];
                ewcy[i][j] = (ex[i] * cen[i][2] +
                              ex[ipr + j - 1] * cen[ipr + j - 1][2]) *
                             exs[i][j];
                ewcz[i][j] = (ex[i] * cen[i][3] +
                              ex[ipr + j - 1] * cen[ipr + j - 1][3]) *
                             exs[i][j];
            }
        }
        // BEGIN LOOP OVER ESP POINTS
        for (iesp = 1; iesp <= nesp; ++iesp) {
            potp1 = potpt[1][iesp] / a0;
            potp2 = potpt[2][iesp] / a0;
            potp3 = potpt[3][iesp] / a0;
            // BEGIN LOOP OVER COMPONENTS OF CONTRACTED FUNCTION IC
            for (j = 1; j <= icd; ++j) {
                // CALCULATE DISTANCE BETWEEN EXPONENT WEIGHTED AND PROBE POINT
                for (i = istart; i <= ipe; ++i) {
                    u_esp[i][j] =
                        ((ewcx[i][j] - potp1) * (ewcx[i][j] - potp1) +
                         (ewcy[i][j] - potp2) * (ewcy[i][j] - potp2) +
                         (ewcz[i][j] - potp3) * (ewcz[i][j] - potp3)) *
                        exsr[i][j];
                    rnai[i][j] = std::sqrt(pi / u_esp[i][j]);
                }
                // CALCULATE ESP INTEGRALS
                for (i = istart; i <= ipe; ++i) {
                    if (u_esp[i][j] <= tf[0]) {
                        iref = (int)std::lround(u_esp[i][j] * 20.0);
                        ref = 0.05 * iref;
                        res = u_esp[i][j] - ref;
                        term = 1.0;
                        f0[i][j] = 0.0;
                        for (k = 0; k <= 6; ++k) {
                            f = fv[k][iref + 1];
                            ts = f * term * fac[k];
                            term = -term * res;
                            f0[i][j] += ts;
                        }
                    } else {
                        f0[i][j] = rnai[i][j] * 0.5;
                    }
                }
                for (i = np; i <= ipe; ++i) {
                    if (u_esp[i][j] <= tf[1]) {
                        iref = (int)std::lround(u_esp[i][j] * 20.0);
                        ref = 0.05 * iref;
                        res = u_esp[i][j] - ref;
                        term1 = 1.0;
                        f1[i][j] = 0.0;
                        for (k = 0; k <= 6; ++k) {
                            fi = fv[k + 1][iref + 1];
                            ts1 = fi * term1 * fac[k];
                            term1 = -term1 * res;
                            f1[i][j] += ts1;
                        }
                    } else {
                        f1[i][j] = rnai[i][j] * 0.25 / u_esp[i][j];
                    }
                }
                for (i = istart; i <= is_esp; ++i)
                    u_esp[i][j] =
                        2.0 * pi * exs[i][j] * expn[i][j] * f0[i][j];
                np = is_esp + 1;
                for (i = np; i <= ipe; ++i) {
                    rnai[i][j] = 2.0 * pi * exs[i][j] * expn[i][j] * f0[i][j];
                    rnai1[i][j] =
                        2.0 * pi * exs[i][j] * expn[i][j] * f1[i][j];
                }
                // CALCULATE (S||P) ESP INTEGRALS
                if (iam[ipr][1] != 0 || is_esp == ip) continue;
                for (i = np; i <= ipe; ++i)
                    u_esp[i][j] = (ewcx[i][j] - cen[i][1]) * rnai[i][j] -
                                  (ewcx[i][j] - potp1) * rnai1[i][j];
                for (i = ipe + 1 - ipx; i <= ipe; ++i)
                    u_esp[i + ipx][j] =
                        (ewcy[i][j] - cen[i][2]) * rnai[i + ipx][j] -
                        (ewcy[i][j] - potp2) * rnai1[i + ipx][j];
                for (i = ipe + 1 + ipx - ipx2; i <= npr - ipx2; ++i)
                    u_esp[i + ipx2][j] =
                        (ewcz[i][j] - cen[i][3]) * rnai[i + ipx2][j] -
                        (ewcz[i][j] - potp3) * rnai1[i + ipx2][j];
            }
            ips = ic * icd - icd + 1;
            for (i = ic; i <= ncc; ++i) {
                jps = i * icd - icd + 1;
                espi[i][ic] = 0.0;
                for (j = jps; j <= jps + icd - 1; ++j) {
                    double ssum = 0.0;
                    for (k = 1; k <= icd; ++k)
                        ssum += cc[j] * cc[ips + k - 1] * u_esp[j][k];
                    espi[i][ic] += ssum;
                }
                es[iesp] += 2.0 * cespm[indc[i]][indc[ic]] * espi[i][ic];
            }
            es[iesp] -= cespm[indc[ic]][indc[ic]] * espi[ic][ic];

        }
        // WRITE OUT RESTART INFORMATION
        // open(unit=iesr, file=restart_fn, status='UNKNOWN', form='UNFORMATTED')
        // iesps = 0
        // write(iesr) ic, iesps
        // do i=1,nesp: write(iesr) es(i)
        // close(iesr)
    }
}

// ---------------------------------------------------------------------
// ovlp(ic) — overlap integrals of the STO-6G contracted function ic.
// ---------------------------------------------------------------------
void ovlp(int ic) {
    int ipr, istart, i, j, np, ips, jps, k;
    double pi;
    pi = 4.0 * std::atan(1.0);
    ipr = ic * icd - icd + 1;
    istart = ipr;
    for (i = istart; i <= npr; ++i) {
        dx_array[i] = cen[ipr][1] - cen[i][1];
        dy_array[i] = cen[ipr][2] - cen[i][2];
        dz_array[i] = cen[ipr][3] - cen[i][3];
        td[i] = dx_array[i] * dx_array[i] + dy_array[i] * dy_array[i] +
                dz_array[i] * dz_array[i];
    }
    // CALCULATE EXPONENT SUM
    for (i = istart; i <= npr; ++i) {
        for (j = 1; j <= icd; ++j) {
            exs[i][j] = 1.0 / (ex[ipr + j - 1] + ex[i]);
            ce[i][j] = ex[ipr + j - 1] * ex[i] * exs[i][j];
            // CALCULATE EXPONENT WEIGHTED CENTERS
            ewcx[i][j] = (ex[i] * cen[i][1] + ex[ipr + j - 1] * cen[ipr + j - 1][1]) *
                         exs[i][j];
            ewcy[i][j] = (ex[i] * cen[i][2] + ex[ipr + j - 1] * cen[ipr + j - 1][2]) *
                         exs[i][j];
            ewcz[i][j] = (ex[i] * cen[i][3] + ex[ipr + j - 1] * cen[ipr + j - 1][3]) *
                         exs[i][j];
        }
    }
    for (i = 1; i <= npr; ++i) {
        for (j = 1; j <= icd; ++j) {
            expn[i][j] = std::exp(-ce[i][j] * td[i]);
            rnai[i][j] = std::pow(pi * exs[i][j], 1.5) * expn[i][j];
            expn[i][j] = rnai[i][j];
        }
    }
    // CALCULATE (S||P) ESP INTEGRALS
    if (iam[ipr][1] == 0 && is_esp != ip) {
        np = is_esp + 1;
        for (i = np; i <= npr; ++i) {
            for (j = 1; j <= icd; ++j) {
                switch (iam[i][2]) {
                    case 2:
                        rnai[i][j] = (ewcy[i][j] - cen[i][2]) * expn[i][j];
                        break;
                    case 3:
                        rnai[i][j] = (ewcz[i][j] - cen[i][3]) * expn[i][j];
                        break;
                    default:
                        rnai[i][j] = (ewcx[i][j] - cen[i][1]) * expn[i][j];
                        break;
                }
            }
        }
    }
    // CALCULATE (P||S) ESP INTEGRALS
    if (iam[ipr][1] == 1 && is_esp != ip) {
        np = is_esp + 1;
        for (i = istart; i <= npr; ++i) {
            for (j = 1; j <= icd; ++j) {
                switch (iam[ipr + j - 1][2]) {
                    case 2:
                        rnai[i][j] = (ewcy[i][j] - cen[ipr + j - 1][2]) * expn[i][j];
                        break;
                    case 3:
                        rnai[i][j] = (ewcz[i][j] - cen[ipr + j - 1][3]) * expn[i][j];
                        break;
                    default:
                        rnai[i][j] = (ewcx[i][j] - cen[ipr + j - 1][1]) * expn[i][j];
                        break;
                }
            }
        }
        for (i = istart; i <= npr; ++i) {
            for (j = 1; j <= icd; ++j) {
                switch (iam[i][2]) {
                    case 2:
                        rnai[i][j] = (ewcy[i][j] - cen[i][2]) * rnai[i][j];
                        if (iam[ipr + j - 1][2] == iam[i][2])
                            rnai[i][j] += exs[i][j] * 0.5 * expn[i][j];
                        break;
                    case 3:
                        rnai[i][j] = (ewcz[i][j] - cen[i][3]) * rnai[i][j];
                        if (iam[ipr + j - 1][2] == iam[i][2])
                            rnai[i][j] += exs[i][j] * 0.5 * expn[i][j];
                        break;
                    default:
                        rnai[i][j] = (ewcx[i][j] - cen[i][1]) * rnai[i][j];
                        if (iam[ipr + j - 1][2] == iam[i][2])
                            rnai[i][j] += exs[i][j] * 0.5 * expn[i][j];
                        break;
                }
            }
        }
    }
    ips = ic * icd - icd + 1;
    for (i = ic; i <= ncc; ++i) {
        jps = i * icd - icd + 1;
        ovl[ic][i] = 0.0;
        for (j = jps; j <= jps + icd - 1; ++j) {
            double ssum = 0.0;
            for (k = 1; k <= icd; ++k)
                ssum += cc[j] * cc[ips + k - 1] * rnai[j][k];
            ovl[ic][i] += ssum;
        }
        ovl[i][ic] = ovl[ic][i];
    }
}

// ---------------------------------------------------------------------
// setup3() — Stewart's STO-3G expansions (J. Chem. Phys. 52, 431).
// ---------------------------------------------------------------------
void setup3() {
    allz[1][1][1] = 2.227660584;
    allz[2][1][1] = 4.057711562e-1;
    allz[3][1][1] = 1.098175104e-1;
    allc[1][1][1] = 1.543289673e-1;
    allc[2][1][1] = 5.353281423e-1;
    allc[3][1][1] = 4.446345422e-1;
    // 2S
    allz[1][2][1] = 2.581578398;
    allz[2][2][1] = 1.567622104e-1;
    allz[3][2][1] = 6.018332272e-2;
    allc[1][2][1] = -5.994474934e-2;
    allc[2][2][1] = 5.960385398e-1;
    allc[3][2][1] = 4.581786291e-1;
    // 2P
    allz[1][2][2] = 9.192379002e-1;
    allz[2][2][2] = 2.359194503e-1;
    allz[3][2][2] = 8.009805746e-2;
    allc[1][2][2] = 1.623948553e-1;
    allc[2][2][2] = 5.661708862e-1;
    allc[3][2][2] = 4.223071752e-1;
    // 3S
    allz[1][3][1] = 5.641487709e-1;
    allz[2][3][1] = 6.924421391e-2;
    allz[3][3][1] = 3.269529097e-2;
    allc[1][3][1] = -1.782577972e-1;
    allc[2][3][1] = 8.612761663e-1;
    allc[3][3][1] = 2.261841969e-1;
    // 3P
    allz[1][3][2] = 2.692880368;
    allz[2][3][2] = 1.489359592e-1;
    allz[3][3][2] = 5.739585040e-2;
    allc[1][3][2] = -1.061945788e-2;
    allc[2][3][2] = 5.218564264e-1;
    allc[3][3][2] = 5.450015143e-1;
    // 4S
    allz[1][4][1] = 2.267938753e-1;
    allz[2][4][1] = 4.448178019e-2;
    allz[3][4][1] = 2.195294664e-2;
    allc[1][4][1] = -3.349048323e-1;
    allc[2][4][1] = 1.056744667;
    allc[3][4][1] = 1.256661680e-1;
    // 4P
    allz[1][4][2] = 4.859692220e-1;
    allz[2][4][2] = 7.430216918e-2;
    allz[3][4][2] = 3.653340923e-2;
    allc[1][4][2] = -6.147823411e-2;
    allc[2][4][2] = 6.604172234e-1;
    allc[3][4][2] = 3.932639495e-1;
    // 5S
    allz[1][5][1] = 1.080198458e-1;
    allz[2][5][1] = 4.408119382e-2;
    allz[3][5][1] = 2.610811810e-2;
    allc[1][5][1] = -6.617401158e-1;
    allc[2][5][1] = 7.467595004e-1;
    allc[3][5][1] = 7.146490945e-1;
    // 5P
    allz[1][5][2] = 2.127482317e-1;
    allz[2][5][2] = 4.729648620e-2;
    allz[3][5][2] = 2.604865324e-2;
    allc[1][5][2] = -1.389529695e-1;
    allc[2][5][2] = 8.076691064e-1;
    allc[3][5][2] = 2.726029342e-1;
}

// ---------------------------------------------------------------------
// naicap() — (P|P) ESP integrals.
// ---------------------------------------------------------------------
void naicap() {
    int idn, idc, ipx2, np, i, l, j, jstart, iesps, iesp, il, iref, k, ic,
        ipr, istart, in, ir, ir2, ips, jps, mptd;
    double pi, potp1, potp2, potp3, ref, res, term, f, ts, term1, fi, ts1,
        term2, fii, ts2;
    idn = 10;
    mptd = ((ipe - is_esp) * (ipe - is_esp + 1)) / 2;
    std::vector<int> id(ipe - is_esp + 1, 0);
    pf0.assign(mptd + 1, 0.0);
    pf1.assign(mptd + 1, 0.0);
    pf2.assign(mptd + 1, 0.0);
    pexs.assign(mptd + 1, 0.0);
    pce.assign(mptd + 1, 0.0);
    pexpn.assign(mptd + 1, 0.0);
    ptd.assign(mptd + 1, 0.0);
    pewcx.assign(mptd + 1, 0.0);
    pewcy.assign(mptd + 1, 0.0);
    pewcz.assign(mptd + 1, 0.0);
    idc = 0;
    ipx2 = 2 * ipx;
    pi = 4.0 * std::atan(1.0);
    np = is_esp + 1;
    // SETUP INDEX ARRAY
    for (i = np; i <= ipe; ++i) {
        ird[i] = i - is_esp;
        ird[i + ipx] = i - is_esp;
        ird[i + ipx2] = i - is_esp;
    }
    // CALCULATE QUANTITIES INVARIANT WITH ESP POINT FOR (P|P) ESP INTEGRALS
    l = 0;
    for (i = np; i <= ipe; ++i) {
        for (j = i; j <= ipe; ++j) {
            ++l;
            ptd[l] = (cen[i][1] - cen[j][1]) * (cen[i][1] - cen[j][1]) +
                     (cen[i][2] - cen[j][2]) * (cen[i][2] - cen[j][2]) +
                     (cen[i][3] - cen[j][3]) * (cen[i][3] - cen[j][3]);
            pexs[l] = 1.0 / (ex[i] + ex[j]);
            pce[l] = ex[i] * ex[j] * pexs[l];
            pexpn[l] = std::exp(-pce[l] * ptd[l]);
            pewcx[l] = (ex[i] * cen[i][1] + ex[j] * cen[j][1]) * pexs[l];
            pewcy[l] = (ex[i] * cen[i][2] + ex[j] * cen[j][2]) * pexs[l];
            pewcz[l] = (ex[i] * cen[i][3] + ex[j] * cen[j][3]) * pexs[l];
        }
        // SET UP OTHER INDEX ARRAY FOR PACKED SYMMETRIC ARRAY STORAGE
        id[i - is_esp] = l - ipx;
    }
    // READ IN RESTART INFORMATION IF THIS IS A RESTART
    if (keywrd.find("ESPRST") != std::string::npos) {
        // open(unit=iesr, file=esr_fn, status='UNKNOWN', form='UNFORMATTED')
        // read(iesr) jstart, iesps
        // if (jstart /= isc*2) { iesps = 0; close(iesr); goto 50; }
        // do i=1,nesp: read(iesr) es(i)
        // close(iesr); idc = int((iesps/float(nesp))*10)
        iesps = 0;
        jstart = isc * 2;
        goto g50;
    } else {
        iesps = 0;
    }
g50:
    for (iesp = iesps + 1; iesp <= nesp; ++iesp) {
        potp1 = potpt[1][iesp] / a0;
        potp2 = potpt[2][iesp] / a0;
        potp3 = potpt[3][iesp] / a0;
        // CALCULATE QUANTITY U
        l = 0;
        for (i = np; i <= ipe; ++i) {
            for (j = i; j <= ipe; ++j) {
                ++l;
                ptd[l] = ((pewcx[l] - potp1) * (pewcx[l] - potp1) +
                          (pewcy[l] - potp2) * (pewcy[l] - potp2) +
                          (pewcz[l] - potp3) * (pewcz[l] - potp3)) /
                         pexs[l];
                pce[l] = std::sqrt(pi / ptd[l]);
            }
        }
        // CALCULATE F0, F1, AND F2(U) USING TAYLOR SERIES OR ASYMPTOTIC
        // EXPANSION
        il = l;
        l = 0;
        for (i = 1; i <= il; ++i) {
            if (ptd[i] <= tf[0]) {
                iref = (int)std::lround(ptd[i] * 20.0);
                ref = 0.05 * iref;
                res = ptd[i] - ref;
                term = 1.0;
                pf0[i] = 0.0;
                for (k = 0; k <= 6; ++k) {
                    f = fv[k][iref + 1];
                    ts = f * term * fac[k];
                    term = -term * res;
                    pf0[i] += ts;
                }
            } else {
                pf0[i] = pce[i] * 0.5;
            }
            if (ptd[i] <= tf[1]) {
                iref = (int)std::lround(ptd[i] * 20.0);
                ref = 0.05 * iref;
                res = ptd[i] - ref;
                term1 = 1.0;
                pf1[i] = 0.0;
                for (k = 0; k <= 6; ++k) {
                    fi = fv[k + 1][iref + 1];
                    ts1 = fi * term1 * fac[k];
                    term1 = -term1 * res;
                    pf1[i] += ts1;
                }
            } else {
                pf1[i] = pce[i] * 0.25 / ptd[i];
            }
            if (ptd[i] <= tf[2]) {
                iref = (int)std::lround(ptd[i] * 20.0);
                ref = 0.05 * iref;
                res = ptd[i] - ref;
                term2 = 1.0;
                pf2[i] = 0.0;
                for (k = 0; k <= 6; ++k) {
                    fii = fv[k + 2][iref + 1];
                    ts2 = fii * term2 * fac[k];
                    term2 = -term2 * res;
                    pf2[i] += ts2;
                }
            } else {
                pf2[i] = pce[i] * 0.375 / (ptd[i] * ptd[i]);
            }
        }
        // CALCULATE (S||S) TYPE INTEGRALS
        for (i = 1; i <= il; ++i) pf0[i] = 2.0 * pi * pexs[i] * pexpn[i] * pf0[i];
        for (i = 1; i <= il; ++i) ptd[i] = pf0[i];
        for (i = 1; i <= il; ++i)
            pf1[i] = 2.0 * pi * pexs[i] * pexpn[i] * pf1[i];
        for (i = 1; i <= il; ++i)
            pf2[i] = 2.0 * pi * pexs[i] * pexpn[i] * pf2[i];
        for (ic = isc + 1; ic <= ncc; ++ic) {
            ipr = ic * icd - icd + 1;
            istart = ipr;
            for (j = 1; j <= icd; ++j) {
                // CALCULATE (P||S) ESP INTEGRALS
                if (iam[ipr][1] == 1 && is_esp != ip) {
                    for (i = istart; i <= npr; ++i) {
                        in = ipr + j - 1;
                        ir = ird[i] + id[ird[in]];
                        ir2 = id[ird[i]] + ird[in];
                        ir = std::min(ir2, ir);
                        switch (iam[in][2]) {
                            case 2:
                                rnai2[i][j] =
                                    (pewcy[ir] - cen[in][2]) * pf1[ir] -
                                    pf2[ir] * (pewcy[ir] - potp2);
                                rnai[i][j] =
                                    (pewcy[ir] - cen[in][2]) * pf0[ir] -
                                    pf1[ir] * (pewcy[ir] - potp2);
                                break;
                            case 3:
                                rnai2[i][j] =
                                    (pewcz[ir] - cen[in][3]) * pf1[ir] -
                                    pf2[ir] * (pewcz[ir] - potp3);
                                rnai[i][j] =
                                    (pewcz[ir] - cen[in][3]) * pf0[ir] -
                                    pf1[ir] * (pewcz[ir] - potp3);
                                break;
                            default:
                                rnai2[i][j] =
                                    (pewcx[ir] - cen[in][1]) * pf1[ir] -
                                    pf2[ir] * (pewcx[ir] - potp1);
                                rnai[i][j] =
                                    (pewcx[ir] - cen[in][1]) * pf0[ir] -
                                    pf1[ir] * (pewcx[ir] - potp1);
                                break;
                        }
                    }
                }
                // CALCULATE (P||P) ESP INTEGRALS
                if (iam[ipr][1] != 1 || is_esp == ip) continue;
                for (i = istart; i <= npr; ++i) {
                    in = ipr + j - 1;
                    ir = ird[i] + id[ird[in]];
                    ir2 = id[ird[i]] + ird[in];
                    ir = std::min(ir2, ir);
                    switch (iam[i][2]) {
                        case 2:
                            rnai[i][j] =
                                (pewcy[ir] - cen[i][2]) * rnai[i][j] -
                                (pewcy[ir] - potp2) * rnai2[i][j];
                            if (iam[in][2] == iam[i][2])
                                rnai[i][j] +=
                                    pexs[ir] * 0.5 * (ptd[ir] - pf1[ir]);
                            break;
                        case 3:
                            rnai[i][j] =
                                (pewcz[ir] - cen[i][3]) * rnai[i][j] -
                                (pewcz[ir] - potp3) * rnai2[i][j];
                            if (iam[in][2] == iam[i][2])
                                rnai[i][j] +=
                                    pexs[ir] * 0.5 * (ptd[ir] - pf1[ir]);
                            break;
                        default:
                            rnai[i][j] =
                                (pewcx[ir] - cen[i][1]) * rnai[i][j] -
                                (pewcx[ir] - potp1) * rnai2[i][j];
                            if (iam[in][2] == iam[i][2])
                                rnai[i][j] +=
                                    pexs[ir] * 0.5 * (ptd[ir] - pf1[ir]);
                            break;
                    }
                }
            }
            // FORM INTEGRALS OVER CONTRACTED FUNCTIONS
            ips = ic * icd - icd + 1;
            for (i = ic; i <= ncc; ++i) {
                jps = i * icd - icd + 1;
                espi[i][ic] = 0.0;
                for (j = jps; j <= jps + icd - 1; ++j) {
                    double ssum = 0.0;
                    for (k = 1; k <= icd; ++k)
                        ssum += cc[j] * cc[ips + k - 1] * rnai[j][k];
                    espi[i][ic] += ssum;
                }
                es[iesp] += 2.0 * cespm[indc[i]][indc[ic]] * espi[i][ic];
            }
            es[iesp] -= cespm[indc[ic]][indc[ic]] * espi[ic][ic];

        }
        // WRITE OUT RESTART INFORMATION EVERY NESP/10 POINTS
        // if (mod(iesp,nesp/idn) /= 0) continue;
        // open(unit=iesr, file=esr_fn, status='UNKNOWN', form='UNFORMATTED')
        // jstart = isc*2
        // write(iesr) jstart, iesp
        // do i=1,nesp: write(iesr) es(i)
        // close(iesr)
        // idc = idc + 1
    }
}

// ---------------------------------------------------------------------
// setup_esp(mode) — allocate/fill the esp_C arrays.
// ---------------------------------------------------------------------
void setup_esp(int mode) {
    int i;
    const int mesp = 5000;
    if (mode == 1) {
        i = 6 * norbs;
        a2.assign(numat + 5, std::vector<double>(numat + 5, 0.0));
        ovl.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
        cespm2.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
        cespml.assign((size_t)norbs * norbs + 1, 0.0);
        cesp.assign((size_t)norbs * norbs + 1, 0.0);
        potpt.assign(4, std::vector<double>(mesp + 1, 0.0));
        es.assign(mesp + 1, 0.0);
        co.assign(4, std::vector<double>(numat + 1, 0.0));
        al.assign((size_t)(numat + 4) * (numat + 4) + 1, 0.0);
        rad.assign(numat + 1, 0.0);
        b_esp.assign(numat + 5, 0.0);
        qesp.assign(numat + 5, 0.0);
        cespm.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
        cc.assign(i + 1, 0.0);
        cen.assign(i + 1, std::vector<double>(4, 0.0));
        iam.assign(i + 1, std::vector<int>(3, 0));
        ind.assign(i + 1, 0);
        ex.assign(i + 1, 0.0);
        temp.assign(i + 1, 0.0);
        f0.assign(i + 1, std::vector<double>(7, 0.0));
        f1.assign(i + 1, std::vector<double>(7, 0.0));
        ce.assign(i + 1, std::vector<double>(7, 0.0));
        u_esp.assign(i + 1, std::vector<double>(7, 0.0));
        exs.assign(i + 1, std::vector<double>(7, 0.0));
        expn.assign(i + 1, std::vector<double>(7, 0.0));
        rnai.assign(i + 1, std::vector<double>(7, 0.0));
        ewcx.assign(i + 1, std::vector<double>(7, 0.0));
        ewcy.assign(i + 1, std::vector<double>(7, 0.0));
        ewcz.assign(i + 1, std::vector<double>(7, 0.0));
        rnai1.assign(i + 1, std::vector<double>(7, 0.0));
        rnai2.assign(i + 1, std::vector<double>(7, 0.0));
        fc.assign(36 * norbs + 1, 0.0);
        exsr.assign(i + 1, std::vector<double>(7, 0.0));
        indc.assign(norbs + 1, 0);
        itemp.assign(i + 1, 0);
        dx_array.assign(i + 1, 0.0);
        dy_array.assign(i + 1, 0.0);
        dz_array.assign(i + 1, 0.0);
        td.assign(i + 1, 0.0);
        esp_array.assign(mesp + 1, 0.0);
        qsc.assign(numat + 1, 0.0);
        espi.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
        ird.assign(i + 1, 0);
        for (int jj = 1; jj <= mesp; ++jj) es[jj] = 1.0e9;
        for (auto& row : a2)
            for (auto& v : row) v = 1.0e9;
        for (auto& row : ovl)
            for (auto& v : row) v = 1.0e9;
        for (auto& row : cespm2)
            for (auto& v : row) v = 1.0e9;
        for (int jj = 1; jj <= norbs * norbs; ++jj) cesp[jj] = 1.0e9;
        for (auto& row : co)
            for (auto& v : row) v = 1.0e9;
        for (int jj = 1; jj <= (numat + 4) * (numat + 4); ++jj)
            al[jj] = 1.0e9;
        for (int jj = 1; jj <= numat; ++jj) rad[jj] = 1.0e9;
        for (int jj = 1; jj <= numat + 4; ++jj) b_esp[jj] = 1.0e9;
        for (int jj = 1; jj <= numat + 4; ++jj) qesp[jj] = 1.0e9;
        for (auto& row : cespm)
            for (auto& v : row) v = 1.0e9;
        for (int jj = 1; jj <= i; ++jj) cc[jj] = 1.0e9;
    }
}

// ---------------------------------------------------------------------
// getattrib(xmin, xmax) — Williams box limits (shared with pdgrid).
// ---------------------------------------------------------------------
void getattrib(double xmin[4], double xmax[4]) {
    int icntr, i, j, ia;
    double vderw[54] = {0.0};
    double shell, vdmax;
    vderw[1] = 2.4;
    vderw[5] = 3.0;
    vderw[6] = 2.9;
    vderw[7] = 2.7;
    vderw[8] = 2.6;
    vderw[9] = 2.55;
    vderw[15] = 3.1;
    vderw[16] = 3.05;
    vderw[17] = 3.0;
    vderw[35] = 3.15;
    vderw[53] = 3.35;
    shell = 1.2;
    nesp = 0;
    // CONVERT INTERNAL TO CARTESIAN COORDINATES
    gmetry(geo, coord);
    // STRIP COORDINATES AND ATOM LABEL FOR DUMMIES (I.E. 99)
    icntr = 0;
    if (co.size() < (size_t)natoms + 1)
        co.assign(4, std::vector<double>(natoms + 1, 0.0));
    for (i = 1; i <= natoms; ++i) {
        co[1][i] = coord[0][i];
        co[2][i] = coord[1][i];
        co[3][i] = coord[2][i];
        if (labels[i] == 99) continue;
        ++icntr;
        nat[icntr] = labels[i];
    }
    numat = icntr;
    for (i = 1; i <= numat; ++i) {
        j = nat[i];
        if (vderw[j] == 0.0) {
            mopend("VAN DER WAALS' RADIUS NOT DEFINED FOR ATOM IN WILLIAMS "
                   "SURFACE ROUTINE PDGRID!");
            return;
        }
    }
    // NOW CREATE LIMITS FOR A BOX
    for (int k = 1; k <= 3; ++k) {
        xmin[k] = 100000.0;
        xmax[k] = -100000.0;
    }
    for (ia = 1; ia <= numat; ++ia) {
        for (int k = 1; k <= 3; ++k) {
            if (co[k][ia] - xmin[k] < 0.0) xmin[k] = co[k][ia];
            if (co[k][ia] - xmax[k] > 0.0) xmax[k] = co[k][ia];
        }
    }
    // ADD (OR SUBTRACT) THE MAXIMUM VDERW PLUS SHELL
    vdmax = 0.0;
    for (i = 1; i <= 53; ++i)
        if (vderw[i] > vdmax) vdmax = vderw[i];
    for (int k = 1; k <= 3; ++k) {
        xmin[k] -= vdmax + shell;
        xmax[k] += vdmax + shell;
    }
}

// ---------------------------------------------------------------------
// espplane(iplane, xmin, step, ngridpts2, ngridpts3) — ESP plane grid.
// ---------------------------------------------------------------------
void espplane(int iplane, const double xmin[4], const double step[4],
              int ngridpts2, int ngridpts3) {
    int ix, iy;
    double ygrid, zgrid;
    if (nesp == 0) {
        for (iy = 1; iy <= ngridpts2; ++iy) {
            ygrid = xmin[2] + step[2] * (iy - 1);
            for (ix = 1; ix <= ngridpts3; ++ix) {
                ++nesp;
                potpt[1][nesp] = xmin[1] + step[1] * (ix - 1);
                potpt[2][nesp] = ygrid;
            }
        }
    }
    zgrid = xmin[3] + step[3] * (iplane - 1);
    for (ix = 1; ix <= nesp; ++ix) potpt[3][ix] = zgrid;
}

// ---------------------------------------------------------------------
// potcal() — nuclear + electronic ESP, then fit; print charges.
// ---------------------------------------------------------------------
void potcal() {
    int i, j, ieq;
    double dipx = 0.0, dipy = 0.0, dipz = 0.0, slope, dip;
    elesn();
    rms = 0.0;
    rrms = 0.0;
    iz = 0;
    std::string::size_type pos;
    if ((pos = keywrd.find("CHARGE=")) != std::string::npos)
        iz = (int)std::lround(reada(keywrd, (int)pos));
    // DIPOLAR CONSTRAINTS IF DESIRED
    if (keywrd.find("DIPOLE") != std::string::npos) {
        idip = 1;
        if (iz != 0) idip = 0;
    } else {
        idip = 0;
    }
    // GET X,Y,Z DIPOLE COMPONENTS IF DESIRED
    if ((pos = keywrd.find("DIPX=")) != std::string::npos)
        dx = reada(keywrd, (int)pos);
    else
        dx = ux;
    if ((pos = keywrd.find("DIPY=")) != std::string::npos)
        dy = reada(keywrd, (int)pos);
    else
        dy = uy;
    if ((pos = keywrd.find("DIPZ=")) != std::string::npos)
        dz = reada(keywrd, (int)pos);
    else
        dz = uz;
    espfit();
    // WRITE OUT OUR RESULTS
    if (!method_mndo) {
        for (i = 1; i <= numat; ++i) {
            // write(iw,'(17X,I4,9X,A2,1X,F10.4)') i, elemnt(nat(i)), q(i)
        }
        to_screen("To_file: Esp");
    } else {
        // MNDO CALCULATION-SCALE THE CHARGES. TEST FOR SLOPE KEYWORD
        if ((pos = keywrd.find("SLOPE=")) != std::string::npos)
            slope = reada(keywrd, (int)pos);
        else
            slope = 1.422;
        for (i = 1; i <= numat; ++i) qsc[i] = slope * q[i];
        // write(iw,'(7X,''ATOM NO.    TYPE    CHARGE   SCALED CHARGE'')')
        // do i=1,numat: write(iw,'(9X,I4,9X,A2,1X,F10.4,2X,F10.4)') i,
        //               elemnt(nat(i)), q(i), qsc(i)
    }
    // write(iw,'(/12X,A,4X,I6)') 'THE NUMBER OF POINTS IS:', nesp
    // write(iw,'(12X,A,4X,F9.4)') 'THE RMS DEVIATION IS:', rms
    // write(iw,'(12X,A,3X,F9.4)') 'THE RRMS DEVIATION IS:', rrms
    // CALCULATE DIPOLE MOMENT IF NEUTRAL MOLECULE
    if (iz != 0) {
        // goto 60
    } else {
        // write(iw,40)
        dipx = 0.0;
        dipy = 0.0;
        dipz = 0.0;
        for (i = 1; i <= numat; ++i) {
            dipx += co[1][i] * q[i] / a0;
            dipy += co[2][i] * q[i] / a0;
            dipz += co[3][i] * q[i] / a0;
        }
        dip = std::sqrt(dipx * dipx + dipy * dipy + dipz * dipz);
        // write(iw,'(12X,'' X        Y        Z       TOTAL'')')
        // write(iw,'(8X,4F9.4)') dipx*cf, dipy*cf, dipz*cf, dip*cf
    }
    // 60 continue
    if (keywrd.find("SYMAVG") != std::string::npos) {
        cequiv.assign(numat + 1, std::vector<bool>(numat + 1, false));
        for (i = 1; i <= numat; ++i) {
            for (j = 1; j <= numat; ++j) {
                cequiv[i][j] = false;
                if (std::fabs(std::fabs(q[i]) - std::fabs(q[j])) >= 1.0e-5)
                    continue;
                cequiv[i][j] = true;
            }
        }
        for (i = 1; i <= numat; ++i) {
            ieq = 0;
            qsc[i] = 0.0;
            for (j = 1; j <= numat; ++j) {
                if (!cequiv[i][j]) continue;
                qsc[i] += std::fabs(q[j]);
                ++ieq;
            }
            q[i] = q[i] / std::fabs(q[i]) * qsc[i] / ieq;
        }
        // write(iw,*) ' '
        // write(iw,*) '   ELECTROSTATIC POTENTIAL CHARGES AVERAGED FOR'
        // write(iw,*) '   SYMMETRY EQUIVALENT ATOMS'
        // write(iw,*) ' '
        if (keywrd.find("AM1") != std::string::npos ||
            keywrd.find("PM3") != std::string::npos) {
            // write(iw,'(7X,''ATOM NO.    TYPE    CHARGE'')')
            // do i=1,numat: write(iw,'(9X,I4,9X,A2,1X,F10.4)') i,
            //               elemnt(nat(i)), q(i)
        } else {
            // write(iw,'(7X,''ATOM NO.    TYPE     CHARGE   SCALED CHARGE'')')
            // do i=1,numat: write(iw,'(9X,I4,9X,A2,1X,F10.4,2X,F10.4)') i,
            //               elemnt(nat(i)), q(i), q(i)*slope
        }
    }
}

double dex2(int m) {
    if (m < 2) return 1.0;
    double r = 1.0;
    for (int i = 1; i <= m; i += 2) r *= i;
    return r;
}

void fsub(int n, double x, double& fval) {
    const double a1s2 = 0.5, pie4 = 0.7853981633974483, a1 = 1.0, xsw = 24.0;
    double e = a1s2 * std::exp(-x);
    double fac0 = n + a1s2;
    std::vector<double> term(200, 0.0);
    double sum, term0, fact, suma, sum1;
    int ku, i;
    if (x > xsw) {
        sum = std::sqrt(pie4 / x);
        if (n != 0) {
            fact = -a1s2;
            for (int k = 1; k <= n; ++k) {
                fact += a1;
                sum *= fact / x;
            }
        }
        term[1] = -e / x;
        suma = sum + term[1];
        if (std::fabs(sum - suma) < 1e-20) {
            fval = suma;
            return;
        }
        fact = fac0;
        ku = (int)(x + fac0 - a1);
        for (i = 2; i <= ku; ++i) {
            fact -= a1;
            term[i] = term[i - 1] * fact / x;
            sum1 = suma;
            suma += term[i];
            if (std::fabs(sum - suma) < 1e-20) break;
        }
    } else {
        fact = fac0;
        term0 = e / fact;
        sum = term0;
        ku = (int)(x - fac0);
        if (ku >= 1)
            for (int k = 1; k <= ku; ++k) {
                fact += a1;
                term0 = term0 * x / fact;
                sum += term0;
            }
        i = 1;
        fact += a1;
        term[1] = term0 * x / fact;
        suma = sum + term[1];
        if (std::fabs(sum - suma) < 1e-20) {
            fval = suma;
            return;
        }
        i = 2;
        fact += a1;
        term[2] = term[1] * x / fact;
        sum1 = suma;
        suma += term[2];
        while (sum1 - suma != 0.0) {
            ++i;
            fact += a1;
            term[i] = term[i - 1] * x / fact;
            sum1 = suma;
            suma += term[i];
        }
    }
    fval = suma;
}
