// local.cpp
#include "local.h"
#include <cmath>
#include <vector>
#include <string>
#include "molkst_C.h"
#include "common_arrays_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

extern "C" void resolv_(double*, double*, int, double*, int) {}
extern "C" void phase_lock_(double*, int) {}
extern "C" void molval_(double*, double*, double) {}
extern "C" void matout_(double*, double*, int, int, int) {}
extern "C" void to_screen_(const char*) {}

void local(double* c, int nocc, double* eig, int iprint, const char* txt) {
    int niter = 100;
    double eps = 1e-7;
    std::vector<double> refeig(eig, eig + norbs);
    std::vector<double> cold(c, c + norbs*norbs);
    std::vector<double> eig1(norbs+1, 0.0);
    int iter = 0;
    while (true) {
        double sum = 0.0;
        ++iter;
        for (int i = 1; i <= nocc; ++i) {
            for (int j = 1; j <= nocc; ++j) {
                if (j == i) continue;
                double xijjj=0, xjiii=0, xiiii=0, xjjjj=0, xijij=0, xiijj=0;
                std::vector<double> psi1(norbs+1), psi2(norbs+1);
                for (int k = 1; k <= norbs; ++k) {
                    psi1[k] = c[(i-1)*norbs + k-1];
                    psi2[k] = c[(j-1)*norbs + k-1];
                }
                for (int k1 = 1; k1 <= numat; ++k1) {
                    int kl = nfirst[k1], ku = nlast[k1];
                    double dij=0, dii=0, djj=0;
                    for (int k = kl; k <= ku; ++k) {
                        dij += psi1[k]*psi2[k];
                        dii += psi1[k]*psi1[k];
                        djj += psi2[k]*psi2[k];
                    }
                    xijjj += dij*djj;
                    xjiii += dij*dii;
                    xiiii += dii*dii;
                    xjjjj += djj*djj;
                    xijij += dij*dij;
                    xiijj += dii*djj;
                }
                double aij = xijij - (xiiii + xjjjj - 2.0*xiijj)/4.0;
                double bij = xjiii - xijjj;
                double ca = std::sqrt(aij*aij + bij*bij);
                double sa = aij + ca;
                if (sa < 1e-14) continue;
                sum += sa;
                ca = -aij/ca;
                ca = (1.0 + std::sqrt((1.0 + ca)/2.0))/2.0;
                if ((2.0*ca - 1.0)*bij < 0.0) ca = 1.0 - ca;
                sa = std::sqrt(1.0 - ca);
                ca = std::sqrt(ca);
                for (int k = 1; k <= norbs; ++k) {
                    double newi = ca*psi1[k] + sa*psi2[k];
                    double newj = -sa*psi1[k] + ca*psi2[k];
                    c[(i-1)*norbs + k-1] = newi;
                    c[(j-1)*norbs + k-1] = newj;
                }
            }
        }
        double sum1 = 0.0;
        for (int i = 1; i <= nocc; ++i)
            for (int j = 1; j <= numat; ++j) {
                int il = nfirst[j], iu = nlast[j];
                double x = 0.0;
                for (int k = il; k <= iu; ++k) x += c[(i-1)*norbs + k-1]*c[(i-1)*norbs + k-1];
                sum1 += x*x;
            }
        if (sum <= eps || iter >= niter) break;
    }
    resolv_(c, cold.data(), norbs, eig, nocc);
    for (int i = 1; i <= nocc; ++i) {
        double sum = 0.0;
        for (int j = 1; j <= nocc; ++j) {
            double co = 0.0;
            for (int k = 1; k <= norbs; ++k)
                co += cold[(j-1)*norbs + k-1]*c[(i-1)*norbs + k-1];
            sum += co*co*eig[j-1];
        }
        eig1[i] = sum;
    }
    for (int i = 1; i <= nocc; ++i) {
        double x = 1e10;
        int i1 = i;
        for (int j = i; j <= nocc; ++j)
            if (eig1[j] < x) { x = eig1[j]; i1 = j; }
        eig[i-1] = eig1[i1];
        std::swap(eig1[i1], eig1[i]);
        for (int j = 1; j <= norbs; ++j) {
            std::swap(c[(i1-1)*norbs + j-1], c[(i-1)*norbs + j-1]);
        }
    }
    if (iprint == 1) {
        phase_lock_(c, norbs);
        matout_(c, eig, nocc, norbs, norbs);
    }
    std::string t(txt, txt+2);
    if (t == "c ") {
        if (nbeta == 0) molval_(c, p.data(), 2.0);
        else molval_(c, pa.data(), 2.0);
    } else {
        molval_(c, pb.data(), 2.0);
    }
    // F90 tail: unless GRAPH is present, restore the input
    // eigenvectors/energies (localization is only used as a
    // diagnostic in that case).
    if (keywrd.find(" GRAPH") == std::string::npos) {
        for (int i = 1; i <= nocc; ++i) eig[i - 1] = refeig[i - 1];
        for (int j = 0; j < norbs * norbs; ++j)
            c[j] = cold[j];
    }
}
