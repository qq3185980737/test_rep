// charmo.cpp — C++ translation of MOPAC 2016 "charmo.F90".

#include "charmo.h"

#include <vector>

#include "dtrans.h"
#include "molkst_C.h"
#include "symmetry_C.h"

using namespace symmetry_C;

double charmo(const std::vector<std::vector<double>>& vects,
              const std::vector<int>& ntype, int jorb, int ioper,
              const std::vector<std::vector<double>>& r, int nvecs, bool& first) {
    if (ioper == 1) return 1.0;

    std::vector<double> vect1(molkst_C::norbs + 1, 0.0),
        vect2(molkst_C::norbs + 1, 0.0);
    std::vector<std::vector<int>> ip(3, std::vector<int>(4, 0));
    std::vector<std::vector<int>> id(3, std::vector<int>(6, 0));
    std::vector<std::vector<int>> loc(3, std::vector<int>(10, 0));
    std::vector<double> h(6, 0.0), p(4, 0.0), d(6, 0.0);

    for (int iatom = 1; iatom <= molkst_C::numat; ++iatom) {
        int jatom = jelem[ioper][iatom];
        int ibase = 0, kj = 0;
        for (int i = 1; i <= nvecs; ++i) {
            int icheck = ntype[i] / 100;
            if (icheck == iatom) {
                ++ibase;
                loc[1][ibase] = i;
            }
            if (icheck != jatom) continue;
            ++kj;
            loc[2][kj] = i;
        }
        if (ibase > 0) {
            int icheck = loc[1][1];
            int jcheck = loc[2][1];
        vect1[icheck] = vects[icheck][jorb];
            vect2[jcheck] = vects[icheck][jorb];
        }
        if (ibase < 4) continue;

        for (int i = 1; i <= 3; ++i) { ip[1][i] = 0; ip[2][i] = 0; p[i] = 0.0; }
        for (int i = 1; i <= 5; ++i) { id[1][i] = 0; id[2][i] = 0; d[i] = 0.0; }
        for (int i = 2; i <= ibase; ++i) {
            int icheck = loc[1][i];
            if (i <= 4) {
                p[i - 1] = vects[icheck][jorb];
                ip[1][i - 1] = loc[1][i];
                ip[2][i - 1] = loc[2][i];
            } else {
                d[i - 4] = vects[icheck][jorb];
                id[1][i - 4] = loc[1][i];
                id[2][i - 4] = loc[2][i];
            }
        }
        if (ibase != 1) {
            h[1] = r[1][1] * p[1] + r[2][1] * p[2] + r[3][1] * p[3];
            h[2] = r[1][2] * p[1] + r[2][2] * p[2] + r[3][2] * p[3];
            h[3] = r[1][3] * p[1] + r[2][3] * p[2] + r[3][3] * p[3];
            for (int i = 1; i <= 3; ++i) {
                double s = 0.0;
                for (int m = 1; m <= 3; ++m) s += elem[i][m][ioper] * h[m];
                p[i] = s;
            }
            for (int i = 1; i <= 3; ++i) {
                if (ip[1][i] < 1) return 0.0;
                int ii = ip[1][i];
                int jj = ip[2][i];
                vect1[ii] = h[i];
                vect2[jj] = p[i];
            }
        }
        if (ibase != 9) continue;
        dtrans(d, ioper, first, r);
        for (int i = 1; i <= 5; ++i) {
            if (id[1][i] < 1) return 0.0;
            int ii = id[1][i];
            int jj = id[2][i];
            vect1[ii] = h[i];
            vect2[jj] = d[i];
        }
    }

    double sum = 0.0;
    for (int k = 1; k <= nvecs; ++k) sum += vect1[k] * vect2[k];
    return sum;
}
