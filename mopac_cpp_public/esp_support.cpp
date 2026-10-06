// esp_support.cpp — helpers for new_esp (MOPAC 2016).
// get_minus_point_five_overlap: direct translation of the routine in
// esp_utilities.F90 — rebuilds the (normalised) overlap from the packed h,
// diagonalises it and forms S^{-1/2} in matrix form.
#include "esp_support.h"

#include <cmath>
#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "rsp.h"

void get_minus_point_five_overlap(std::vector<std::vector<double>>& s) {
    using common_arrays_C::nfirst;
    using common_arrays_C::nlast;
    using common_arrays_C::nat;
    using common_arrays_C::h;
    using common_arrays_C::i1fact;
    using molkst_C::numat;
    using molkst_C::norbs;
    using parameters_C::betas;
    using parameters_C::betap;
    using parameters_C::betad;

    // Lower-half-triangle indices (Pascal's triangle).
    if (static_cast<int>(i1fact.size()) < 3 + norbs) i1fact.assign(3 + norbs, 0);
    for (int i = 1; i <= norbs; ++i) i1fact[i] = (i * (i + 1)) / 2;

    std::vector<double> eigs(norbs + 1, 0.0), vecs(norbs * norbs, 0.0);
    for (int i = 1; i <= numat; ++i) {
        int if_ = nfirst[i];
        int il = nlast[i];
        if (il >= if_) {
            eigs[if_] = betas[nat[i]];
            if (il > if_) {
                eigs[if_ + 1] = betap[nat[i]];
                eigs[if_ + 2] = eigs[if_ + 1];
                eigs[if_ + 3] = eigs[if_ + 1];
                if (il > if_ + 3) {
                    eigs[if_ + 4] = betad[nat[i]];
                    eigs[if_ + 5] = eigs[if_ + 4];
                    eigs[if_ + 6] = eigs[if_ + 4];
                    eigs[if_ + 7] = eigs[if_ + 4];
                    eigs[if_ + 8] = eigs[if_ + 4];
                }
            }
        }
        for (int k = if_; k <= il; ++k) {
            double bi = eigs[k];
            int ii = (k * (k - 1)) / 2;
            for (int j = 1; j <= i - 1; ++j) {
                int jf = nfirst[j];
                int jl = nlast[j];
                for (int jj = jf; jj <= jl; ++jj) {
                    double bj = eigs[jj];
                    int ij = ii + jj;
                    h[ij] = 2.0 * h[ij] / (bi + bj) + 1.0e-14;
                }
            }
            for (int jj = if_; jj <= k; ++jj) {
                int ij = ii + jj;
                h[ij] = 0.0;
            }
        }
    }
    for (int i = 1; i <= norbs; ++i) h[i1fact[i]] = 1.0;
    // rsp() consumes 0-based arrays; h/eigs are 1-based (index 0 is pad).
    std::vector<double> hdiag = h;
    rsp(hdiag.data() + 1, norbs, eigs.data() + 1, vecs.data());
    for (int i = 1; i <= norbs; ++i) eigs[i] = 1.0 / std::sqrt(std::fabs(eigs[i]));
    // s(i,j) = sum_k vecs(i,k) eigs(k) vecs(j,k)   (S^{-1/2}, symmetric)
    if (static_cast<int>(s.size()) <= norbs) s.assign(norbs + 1, std::vector<double>(norbs + 1, 0.0));
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 1; k <= norbs; ++k)
                sum += vecs[(k - 1) * norbs + (i - 1)] * eigs[k] *
                       vecs[(k - 1) * norbs + (j - 1)];
            s[i][j] = sum;
            s[j][i] = sum;
        }
}

void density_for_GPU(const double* vecs, double fract, int nclose, int nopen,
                     double two, int ij, int norbs, int mode,
                     std::vector<double>& p, int mode2) {
  (void)mode;
  (void)mode2;
  if ((int)p.size() < ij + 1) p.resize(ij + 1, 0.0);
  // vecs flat layout from mult(): vecs[(i-1)*norbs + (k-1)] = C(i,k),
  // AO i, MO k (row-major interpretation).
  for (int i = 1; i <= norbs; ++i) {
    for (int j = 1; j <= i; ++j) {
      double sum = 0.0;
      for (int k = 1; k <= nclose; ++k)
        sum += two * vecs[(i - 1) * norbs + (k - 1)] *
                    vecs[(j - 1) * norbs + (k - 1)];
      for (int k = nclose + 1; k <= nopen; ++k)
        sum += fract * vecs[(i - 1) * norbs + (k - 1)] *
                       vecs[(j - 1) * norbs + (k - 1)];
      p[(i * (i - 1)) / 2 + j] = sum;
    }
  }
}

