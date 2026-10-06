// charvi.cpp — C++ translation of MOPAC 2016 "charvi.F90".

#include "charvi.h"

#include "molkst_C.h"
#include "symmetry_C.h"

using namespace symmetry_C;

double charvi(const std::vector<std::vector<double>>& vects, int jorb,
              int ioper, const std::vector<std::vector<double>>& r, int nvecs) {
    if (ioper == 1) return 1.0;

    std::vector<double> vect1(nvecs + 1, 0.0), vect2(nvecs + 1, 0.0);
    for (int iatom = 1; iatom <= molkst_C::numat; ++iatom) {
        int jatom = jelem[ioper][iatom];
        double p1 = vects[iatom * 3 - 2][jorb];
        double p2 = vects[iatom * 3 - 1][jorb];
        double p3 = vects[iatom * 3][jorb];
        double h1 = r[1][1] * p1 + r[2][1] * p2 + r[3][1] * p3;
        double h2 = r[1][2] * p1 + r[2][2] * p2 + r[3][2] * p3;
        double h3 = r[1][3] * p1 + r[2][3] * p2 + r[3][3] * p3;
        double q1 = elem[1][1][ioper] * h1 + elem[1][2][ioper] * h2 + elem[1][3][ioper] * h3;
        double q2 = elem[2][1][ioper] * h1 + elem[2][2][ioper] * h2 + elem[2][3][ioper] * h3;
        double q3 = elem[3][1][ioper] * h1 + elem[3][2][ioper] * h2 + elem[3][3][ioper] * h3;
        vect1[iatom * 3 - 2] = h1;
        vect1[iatom * 3 - 1] = h2;
        vect1[iatom * 3] = h3;
        vect2[jatom * 3 - 2] = q1;
        vect2[jatom * 3 - 1] = q2;
        vect2[jatom * 3] = q3;
    }
    double sum = 0.0;
    for (int k = 1; k <= nvecs; ++k) sum += vect1[k] * vect2[k];
    return sum;
}
