// gdisp.cpp — C++ translation of MOPAC 2016 "gdisp.F90".

#include "gdisp.h"

#include <cmath>
#include <vector>
#include "molkst_C.h"
#include "dftd3_bits.h"
#include "common_arrays_C.h"

void ncoord(int natoms, const std::vector<double>& rcov,
            const std::vector<int>& nat,
            const std::vector<std::vector<double>>& xyz,
            std::vector<double>& cn) {
    for (int i = 1; i <= natoms; ++i) {
        double xn = 0.0;
        for (int j = 1; j <= natoms; ++j) {
            if (j == i) continue;
            double dx = xyz[1][j] - xyz[1][i];
            double dy = xyz[2][j] - xyz[2][i];
            double dz = xyz[3][j] - xyz[3][i];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            double rco = rcov[nat[i]] + rcov[nat[j]];
            double rr = rco / r;
            double damp = 1.0 / (1.0 + std::exp(-16.0 * (rr - 1.0)));
            xn += damp;
        }
        cn[i] = xn;
    }
}

void get_dC6_dCNij(int, int,
                   const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& c6ab,
                   int mxci, int mxcj, double cni, double cnj,
                   int izi, int izj, double& c6check, double& dc6i, double& dc6j) {
    double c6mem = -1e99, r_save = 9999.0;
    double zaehler = 0, nenner = 0, dzi = 0, dni = 0, dzj = 0, dnj = 0;
    for (int a = 1; a <= mxci; ++a)
        for (int b = 1; b <= mxcj; ++b) {
            double c6ref = c6ab[izi][izj][a][b][1];
            if (c6ref > 0.0) {
                double cri = c6ab[izi][izj][a][b][2];
                double crj = c6ab[izi][izj][a][b][3];
                double r = (cri - cni) * (cri - cni) + (crj - cnj) * (crj - cnj);
                if (r < r_save) { r_save = r; c6mem = c6ref; }
                double et = std::exp(-4.0 * r);
                zaehler += c6ref * et;
                nenner += et;
                et = -4.0 * et * 2.0;
                double term = et * (cni - cri);
                dzi += c6ref * term; dni += term;
                term = et * (cnj - crj);
                dzj += c6ref * term; dnj += term;
            }
        }
    if (nenner > 1e-99) {
        c6check = zaehler / nenner;
        dc6i = (dzi * nenner - dni * zaehler) / (nenner * nenner);
        dc6j = (dzj * nenner - dnj * zaehler) / (nenner * nenner);
    } else {
        c6check = c6mem; dc6i = 0; dc6j = 0;
    }
}

void gdisp(const std::vector<std::vector<double>>& xyz,
           const std::vector<std::vector<double>>& r0ab,
           double rs6, double alp6, const c6ab_t& c6ab, double s6,
           const std::vector<int>& mxc,
           const std::vector<double>& rcov,
           std::vector<std::vector<double>>& dxyz_temp) {
    const int max_elem = 94, maxc = 5;
    const double k1 = 16.0;
    using namespace molkst_C;
    using namespace common_arrays_C;
    int npair = (numat * (numat + 1)) / 2;
    std::vector<double> drij(npair + 1, 0.0);
    std::vector<double> dc6i(numat + 1, 0.0);
    std::vector<double> cn(numat + 1, 0.0);
    ncoord(numat, rcov, nat, xyz, cn);
    double dc6_rest = 0.0;
    auto lin = [](int i, int j) { return (i * (i - 1)) / 2 + j; };
    for (int i = 2; i <= numat; ++i) {
        for (int j = 1; j <= i - 1; ++j) {
            double rij[3] = {
                xyz[1][j] - xyz[1][i],
                xyz[2][j] - xyz[2][i],
                xyz[3][j] - xyz[3][i]};
            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
            if (r2 > 10000.0) continue;
            int linij = lin(i, j);
            double R0 = r0ab[nat[j]][nat[i]];
            double c6, dc6iji, dc6ijj;
            get_dC6_dCNij(maxc, max_elem, c6ab, mxc[nat[i]], mxc[nat[j]],
                          cn[i], cn[j], nat[i], nat[j], c6, dc6iji, dc6ijj);
            double r = std::sqrt(r2);
            double r6 = r2 * r2 * r2;
            double r7 = r6 * r;
            double t6 = std::pow(r / (rs6 * R0), -alp6);
            double damp6 = 1.0 / (1.0 + 6.0 * t6);
            double tmp1 = s6 * 6.0 * damp6 * c6 / r7;
            drij[linij] -= tmp1;
            drij[linij] += tmp1 * alp6 * t6 * damp6;
            dc6_rest = s6 / r6 * damp6;
            dc6i[i] += dc6_rest * dc6iji;
            dc6i[j] += dc6_rest * dc6ijj;
        }
    }
    for (int i = 2; i <= numat; ++i) {
        for (int j = 1; j <= i - 1; ++j) {
            int linij = lin(i, j);
            double rij[3] = {
                xyz[1][j] - xyz[1][i],
                xyz[2][j] - xyz[2][i],
                xyz[3][j] - xyz[3][i]};
            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
            double r = std::sqrt(r2);
            double dcn;
            if (r2 < 100.0) {
                double rcovij = rcov[nat[i]] + rcov[nat[j]];
                double expterm = std::exp(-k1 * (rcovij / r - 1.0));
                dcn = -k1 * rcovij * expterm / (r * r * (expterm + 1.0) * (expterm + 1.0));
            } else {
                dcn = 0.0;
            }
            double x1 = drij[linij] + dcn * (dc6i[i] + dc6i[j]);
            dxyz_temp[1][i] += x1 * rij[0] / r;
            dxyz_temp[2][i] += x1 * rij[1] / r;
            dxyz_temp[3][i] += x1 * rij[2] / r;
            dxyz_temp[1][j] -= x1 * rij[0] / r;
            dxyz_temp[2][j] -= x1 * rij[1] / r;
            dxyz_temp[3][j] -= x1 * rij[2] / r;
        }
    }
}
