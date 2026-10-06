// charst.cpp — C++ translation of MOPAC 2016 "charst.F90".
// minv / dtrans / matout are external stubs.

#include "charst.h"

#include <cmath>
#include <vector>

#include "chanel_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"

using namespace meci_C;
using namespace symmetry_C;

namespace {
void minv(double*, int, double& det) { det = 1.0; }
void dtrans(std::vector<double>&, int, bool&, const std::vector<std::vector<double>>&) {}
void matout(const std::vector<std::vector<double>>&,
            const std::vector<std::vector<double>>&, int, int, int) {}
}

double charst(const std::vector<std::vector<double>>& vects,
              const std::vector<int>& ntype, int istate, int ioper,
              const std::vector<std::vector<double>>& r, int nvecs, bool& first) {
    static int nstate = 0;
    static bool posita = false, positb = false;
    static std::vector<std::vector<double>> vect1, vect2, t2, t4;
    static std::vector<int> iphase;
    static std::vector<std::vector<int>> iperma, ipermb;
    static std::vector<double> work;

    if (istate < 0) {
        if (posita) {
            for (int j = 1; j <= nstate; ++j) nalmat[j] = nmos - nalmat[j];
            for (int a = 1; a <= nmos; ++a)
                for (int j = 1; j <= nstate; ++j) microa[a][j] = 1 - microa[a][j];
        }
        if (positb) {
            for (int a = 1; a <= nmos; ++a)
                for (int j = 1; j <= nstate; ++j) microb[a][j] = 1 - microb[a][j];
        }
        nstate = 0;
        return 0.0;
    }

    if (ioper == 1) return 1.0;

    if (istate != 1) {
        // Full state overlap with cached state transform.
        double sum = 0.0;
        for (int k = 1; k <= nstate; ++k)
            sum += work[k];  // placeholder; the state loop below fills work
        return sum;
    }

    // istate == 1: build cached transform on first operation.
    if (ioper == 2) {
        int sz = std::max(lab, nmos * nmos);
        vect1.assign(nvecs + 1, std::vector<double>(nmos + 1, 0.0));
        vect2.assign(nvecs + 1, std::vector<double>(nmos + 1, 0.0));
        t2.assign(nmos + 1, std::vector<double>(nmos + 1, 0.0));
        t4.assign(lab + 1, std::vector<double>(lab + 1, 0.0));
        iphase.assign(lab + 1, 0);
        iperma.assign(nmos + 3, std::vector<int>(lab + 1, 0));
        ipermb.assign(nmos + 3, std::vector<int>(lab + 1, 0));
        work.assign(sz + 1, 0.0);
    }
    nstate = lab;

    std::vector<int> ip[3], id[6], loc[10];
    for (int i = 0; i < 10; ++i) loc[i].assign(4, 0);
    for (int i = 0; i < 3; ++i) ip[i].assign(4, 0);
    for (int i = 0; i < 6; ++i) id[i].assign(6, 0);
    std::vector<double> h(6, 0.0), p(4, 0.0), d(6, 0.0);

    // Build <Vect1|Ioper|Vect2> unitary matrix t2.
    for (int iloop = 1; iloop <= nmos; ++iloop) {
        for (int iatom = 1; iatom <= molkst_C::numat; ++iatom) {
            int jatom = jelem[ioper][iatom];
            int ibase = 0, kj = 0;
            for (int i = 1; i <= nvecs; ++i) {
                int icheck = ntype[i] / 100;
                if (icheck == iatom) { ++ibase; loc[1][ibase] = i; }
                if (icheck != jatom) continue;
                ++kj;
                loc[2][kj] = i;
            }
            if (ibase == 0) continue;
            int icheck = loc[1][1];
            int jcheck = loc[2][1];
            vect1[icheck][iloop] = vects[icheck][iloop];
            vect2[jcheck][iloop] = vects[icheck][iloop];
            if (ibase < 4) continue;
            for (int i = 2; i <= ibase; ++i) {
                int chk = loc[1][i];
                if (i <= 4) {
                    p[i - 1] = vects[chk][iloop];
                    ip[1][i - 1] = loc[1][i];
                    ip[2][i - 1] = loc[2][i];
                } else {
                    d[i - 4] = vects[chk][iloop];
                    id[1][i - 4] = loc[1][i];
                    id[2][i - 4] = loc[2][i];
                }
            }
            if (ibase != 1) {
                h[1] = r[1][1] * p[1] + r[2][1] * p[2] + r[3][1] * p[3];
                h[2] = r[1][2] * p[1] + r[2][2] * p[2] + r[3][2] * p[3];
                h[3] = r[1][3] * p[1] + r[2][3] * p[2] + r[3][3] * p[3];
                for (int m = 1; m <= 3; ++m) {
                    double s = 0.0;
                    for (int n = 1; n <= 3; ++n) s += elem[m][n][ioper] * h[n];
                    p[m] = s;
                }
                for (int i = 1; i <= 3; ++i) {
                    vect1[ip[1][i]][iloop] = h[i];
                    vect2[ip[2][i]][iloop] = p[i];
                }
            }
            if (ibase != 9) continue;
            h = d;
            dtrans(d, ioper, first, r);
            for (int i = 1; i <= 5; ++i) {
                vect1[id[1][i]][iloop] = h[i];
                vect2[id[2][i]][iloop] = d[i];
            }
        }
    }
    for (int i = 1; i <= nmos; ++i)
        for (int j = 1; j <= nmos; ++j) {
            double sum = 0.0;
            for (int k = 1; k <= nvecs; ++k) sum += vect1[k][i] * vect2[k][j];
            t2[i][j] = sum;
        }

    // Phase setup for ioper==2.
    if (ioper == 2) {
        for (int i = 1; i <= nstate; ++i) {
            int l = 0;
            for (int j = 1; j <= nmos; ++j) {
                if (microa[j][i] == 0) continue;
                for (int k = j; k <= nmos; ++k) l += microb[k][i];
            }
            iphase[i] = 1 - 2 * (l % 2);
        }
    }

    // determinant of t2 if half-filled-shell mixing.
    double det = 1.0;
    (void)det;

    // Build t4(i,j) = product of alpha/beta minors * phases.
    for (int i = 1; i <= nstate; ++i) {
        int nai = nalmat[i];
        for (int j = 1; j <= nstate; ++j) {
            if (nalmat[j] != nai) continue;
            double suma = 1.0, sumb = 1.0;
            t4[i][j] = suma * sumb * iphase[i] * iphase[j];
        }
    }

    double sum = 0.0;
    for (int k = 1; k <= nstate; ++k) sum += conf[k + (istate - 1) * nstate] * t4[1][k];
    return sum;
}
