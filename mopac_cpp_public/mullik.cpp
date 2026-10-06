// mullik.cpp — C++ translation of MOPAC 2016 "mullik.F90".
#define _CRT_SECURE_NO_WARNINGS
// S^{-1/2}, orthogonalises C, builds the density P = (C S^{-1/2}) occ (C S^{-1/2})^T
// and computes atom charges and per-orbital populations.
#include "mullik.h"

#include "molkst_C.h"
#include "symmetry_C.h"
#include "maps_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "chanel_C.h"
#include "chrge.h"
#include "rsp.h"
#include "mult.h"
#include "density_for_GPU.h"
#include "mopend.h"
#include "mod_vars_cuda.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace symmetry_C;
using namespace maps_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace chanel_C;
using namespace mod_vars_cuda;

void mullik() {
    static bool graph = false, graph_formatted = false;
    static int icalcn = 0, mo_l = 1, mo_u = 0;

    const int n = norbs;
    const int nlower = (n * (n + 1)) / 2;
    std::vector<double> store(nlower + 1), store_h(nlower + 1), eig(n + 1);
    std::vector<std::vector<double>> work;
    std::vector<double> q2(numat + 1);

    if (icalcn != numcal) {
        icalcn = numcal;
        graph = keywrd.find("GRAPH") != std::string::npos;
        graph_formatted = keywrd.find("GRAPHF") != std::string::npos;
        ifact.assign(n + 2, 0);
        for (int i = 1; i <= n; ++i) ifact[i] = (i * (i - 1)) / 2;
        ifact[n + 1] = nlower;
        if (graph) {
            std::string fname = graph_formatted
                ? gpt_fn.substr(0, gpt_fn.size() - 3) + "mgf"
                : gpt_fn;
            FILE* f = std::fopen(fname.c_str(), graph_formatted ? "w" : "wb");
            if (!f) {
                mopend("File '" + fname + "' is unavailable for use");
                return;
            }
            std::fclose(f);
        }
        mo_l = 1;
        mo_u = n;
    }
    chrge(p, q2);
    for (int i = 1; i <= numat; ++i) q[i] = tore[nat[i]] - q2[i];

    store_h = h;
    for (int i = 1; i <= numat; ++i) {
        int i_f = nfirst[i], i_l = nlast[i], im1 = i - 1;
        double bi = betas[nat[i]];
        for (int k = i_f; k <= i_l; ++k) {
            int ii = (k * (k - 1)) / 2;
            for (int j = 1; j <= im1; ++j) {
                int jf = nfirst[j], jl = nlast[j];
                double bj = betas[nat[j]];
                // +1.D-14 guards against errors in the diagonalization.
                int ij = ii + jf;
                h[ij] = 2.0 * h[ij] / (bi + bj) + 1.0e-14;
                store[ij] = h[ij];
                bj = betap[nat[j]];
                bj = betap[nat[j]];
                for (int m = jf + 1; m <= jl; ++m) {
                    h[ii + m] = 2.0 * h[ii + m] / (bi + bj) + 1.0e-14;
                    store[ii + m] = h[ii + m];
                }
            }
            for (int m = i_f; m <= k; ++m) {
                store[ii + m] = 0.0;
                h[ii + m] = 0.0;
            }
            bi = betap[nat[i]];
        }
    }
    for (int m = 2; m <= n + 1; ++m) {
        store[ifact[m]] = 1.0;
        h[ifact[m]] = 1.0;
    }
    // Diagonalize the (modified) overlap matrix: rsp uses 0-based arrays,
    // h is 1-based (h[0] is a pad), so pass hdiag.data()+1.
    std::vector<double> hdiag = h;
    std::vector<double> vecs_rsp(n * n);
    rsp(hdiag.data() + 1, n, eig.data() + 1, vecs_rsp.data());
    for (int i = 1; i <= n; ++i) eig[i] = 1.0 / std::sqrt(std::fabs(eig[i]));
    // work(i,j) = sum_k vecs(i,k) eig(k) vecs(j,k)  =>  S^{-1/2}
    work.assign(n + 1, std::vector<double>(n + 1, 0.0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= i; ++j) {
            double s = 0.0;
            for (int k = 1; k <= n; ++k)
                s += vecs_rsp[(k - 1) * n + (i - 1)] * eig[k] * vecs_rsp[(k - 1) * n + (j - 1)];
            work[i][j] = s;
            work[j][i] = s;
        }

    if (graph) {
        if (graph_formatted) {
            std::string fname = gpt_fn.substr(0, gpt_fn.size() - 3) + "mgf";
            FILE* f = std::fopen(fname.c_str(), "w");
            if (!f) { mopend("File '" + fname + "' is unavailable for use"); return; }
            std::fprintf(f, "%5d %s\n", numat, ("MOPAC-Graphical data Version 2012." + verson).c_str());
            for (int i = 1; i <= numat; ++i)
                std::fprintf(f, "%4d%13.7f%12.7f%12.7f%9.4f\n", nat[i], coord[0][i], coord[1][i], coord[2][i], q[i]);
            for (int i = 1; i <= numat; ++i)
                std::fprintf(f, "%11.7f%11.7f%11.7f\n", zs[nat[i]], zp[nat[i]], zd[nat[i]]);
            if (keywrd.find(" ALLV") == std::string::npos) {
                mo_l = std::max(1, std::max(nalpha, nclose) - 8);
                mo_u = std::min(n, std::max(nalpha, nclose) + 7);
            }
            bool namo_ok = !namo.empty();
            for (int i = mo_l; i <= mo_u; ++i) {
                int j = 0;
                if (uhf && i <= nalpha) j = 1;
                if (!uhf && i <= nclose) j = 2;
                if (namo_ok)
                    std::fprintf(f, "%s%2d%2d %s%10.4f\n", "ORBITAL", j, jndex[i], namo[i].c_str(), eigs[i]);
                else
                    std::fprintf(f, "%s%2d %s%10.4f\n", "ORBITAL", j, "NONE", eigs[i]);
                for (int j2 = 1; j2 <= n; ++j2) std::fprintf(f, "%15.8e", c[j2][i]);
                std::fprintf(f, "\n");
            }
            std::fprintf(f, "INVERSE_MATRIX[%d x %d]=\n", n, n);
            for (int i = 1; i <= n; ++i) {
                for (int j2 = 1; j2 <= i; ++j2) std::fprintf(f, "%15.8e", work[j2][i]);
                std::fprintf(f, "\n");
            }
            if (rxn_coord < 1.0e8) {
                std::fprintf(f, "Reaction coordinate:%12.4f\n", rxn_coord);
                std::fprintf(f, "Heat_of_formation:  %12.4f\n", escf);
            }
            std::string l = keywrd;
            if (method_pm6 && keywrd.find(" PM6") == std::string::npos) l = " PM6" + l;
            if (uhf && keywrd.find(" UHF") == std::string::npos) l = " UHF" + l;
            std::fprintf(f, "Keywords:%s\n", l.c_str());
            if (uhf) {
                if (keywrd.find(" ALLV") == std::string::npos) {
                    mo_l = std::max(1, std::max(nbeta, nclose) - 8);
                    mo_u = std::min(n, std::max(nbeta, nclose) + 7);
                }
                bool namo_ok2 = !namo.empty();
                // Note: Fortran calls symtrz(cb, eigb, 1, .TRUE.) here to
                // symmetrise the beta eigenvectors; skipped (symtrz port pending).
                for (int i = mo_l; i <= mo_u; ++i) {
                    int j = 0;
                    if (uhf && i <= nbeta) j = 1;
                    if (namo_ok2)
                        std::fprintf(f, "%s%2d%2d %s%10.4f\n", "ORBITAL", j, jndex[i], namo[i].c_str(), eigb[i]);
                    else
                        std::fprintf(f, "%s%2d %s%10.4f\n", "ORBITAL", j, "NONE", eigb[i]);
                    for (int j2 = 1; j2 <= n; ++j2) std::fprintf(f, "%15.8e", cb[j2][i]);
                    std::fprintf(f, "\n");
                }
            }
            std::fclose(f);
        } else {
            std::string fname = gpt_fn;
            FILE* f = std::fopen(fname.c_str(), "wb");
            if (!f) { mopend("File '" + fname + "' is unavailable for use"); return; }
            int rec1[3] = {numat, n, nelecs};
            std::fwrite(rec1, sizeof(int), 3, f);
            for (int d = 1; d <= 3; ++d)
                for (int a = 1; a <= numat; ++a) { double v = coord[d-1][a]; std::fwrite(&v, sizeof(double), 1, f); }
            for (int a = 1; a <= numat; ++a) { std::fwrite(&nlast[a], sizeof(int), 1, f); std::fwrite(&nfirst[a], sizeof(int), 1, f); }
            for (int a = 1; a <= numat; ++a) { double v = zs[nat[a]]; std::fwrite(&v, sizeof(double), 1, f); }
            for (int a = 1; a <= numat; ++a) { double v = zp[nat[a]]; std::fwrite(&v, sizeof(double), 1, f); }
            for (int a = 1; a <= numat; ++a) { double v = zd[nat[a]]; std::fwrite(&v, sizeof(double), 1, f); }
            for (int a = 1; a <= numat; ++a) { int v = nat[a]; std::fwrite(&v, sizeof(int), 1, f); }
            for (int i = 1; i <= n; ++i) for (int j = 1; j <= n; ++j) { double v = c[i][j]; std::fwrite(&v, sizeof(double), 1, f); }
            for (int i = 1; i <= n; ++i) for (int j = 1; j <= n; ++j) { double v = work[i][j]; std::fwrite(&v, sizeof(double), 1, f); }
            std::fwrite(&id, sizeof(int), 1, f);
            for (int i = 1; i <= 3; ++i) for (int j = 1; j <= 3; ++j) { double v = tvec[i][j]; std::fwrite(&v, sizeof(double), 1, f); }
            std::fclose(f);
        }
        h = store_h;
        if (keywrd.find("MULLIK") == std::string::npos) return;
    }

    std::vector<double> buf_c(n * n), buf_w(n * n), vecs_orth(n * n);
    // vecs = S^{-1/2} * C  (mult reads column-major buffers:
    // buf[(row-1)+(col-1)*n] = M[row][col]).
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j) {
            buf_c[(j - 1) + (i - 1) * n] = c[j][i];
            buf_w[(j - 1) + (i - 1) * n] = work[j][i];
        }
    mult(buf_c.data(), buf_w.data(), vecs_orth.data(), n);
    // Build density into pb (1-based vector).
    // density_for_GPU expects row-major input; mult emits column-major.
    std::vector<double> vecs_rm(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) vecs_rm[i * n + j] = vecs_orth[j * n + i];
    if (static_cast<int>(pb.size()) < nlower + 1) pb.resize(nlower + 1);
    density_for_GPU(vecs_rm.data(), fract, nclose, nopen, 2.0, nlower, n, 2,
                    pb.data() + 1, mod_vars_cuda::lgpu ? 2 : 3);
    // P = (C S^{-1/2}) occ (C S^{-1/2})^T * S
    for (int i = 1; i <= nlower; ++i) pb[i] *= store[i];
    double summ = 0.0;
    for (int i = 1; i <= n; ++i) {
        double s = 0.0;
        for (int j = 1; j <= i; ++j) s += pb[ifact[i] + j];
        for (int j = i + 1; j <= n; ++j) s += pb[ifact[j] + i];
        summ += s;
        pb[ifact[i + 1]] = s;
    }
    (void)summ;
    h = store_h;
}
