// freqcy.cpp — C++ translation of MOPAC 2016 "freqcy.F90".
#include "freqcy.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "to_screen_C.h"
#include "symt.h"
#include "frame.h"
#include "rsp.h"
#include "symtrz.h"
#include "force.h"

using namespace molkst_C;
using namespace common_arrays_C;

void freqcy_mass_weight(std::vector<double>& fmatrx,
                        std::vector<double>& oldf,
                        const std::vector<double>& wtmass, int nvar) {
    int linear = 0;
    for (int i = 1; i <= nvar; ++i) {
        for (int k = 1; k <= i; ++k) oldf[linear + k] = fmatrx[linear + k] * 1e5;
        for (int k = 1; k <= i; ++k)
            fmatrx[linear + k] = fmatrx[linear + k] * wtmass[i] * wtmass[k];
        linear += i;
    }
}

void freqcy(std::vector<double>& fmatrx,
            std::vector<double>& freq,
            std::vector<double>& travel,
            bool eorc,
            std::vector<std::vector<double>>& deldip,
            std::vector<double>& ff,
            std::vector<double>& oldf,
            bool ts) {
    double fact = funcon_C::fpc_10;
    double c2pi = 1.0 / (funcon_C::fpc_8 * funcon_C::pi * 2.0);
    if (!ts && keywrd.find(" NOSYM") == std::string::npos) {
        for (int loop = 1; loop <= 20; ++loop) {
            double sumerr = 0.0;
            for (int i = 1; i <= nvar; ++i) {
                double sum = 0.0, err = 0.0;
                for (int j = 1; j <= i - 1; ++j)
                    sum += fmatrx[(i * (i - 1)) / 2 + j];
                for (int j = i + 1; j <= nvar; ++j)
                    sum += fmatrx[(j * (j - 1)) / 2 + i];
                err += fmatrx[(i * (i + 1)) / 2] + sum;
                sumerr += std::abs(err);
                fmatrx[(i * (i + 1)) / 2] = (-sum) - err * 0.5;
            }
            // F90: call symt(fmatrx, deldip, ff); deldip(3,3*numat) column-major.
            std::vector<double> deldip_flat(9 * numat + 1, 0.0);
            for (int i = 1; i <= 3; ++i)
                for (int j = 1; j <= 3 * numat; ++j)
                    deldip_flat[(j - 1) * 3 + (i - 1)] = deldip[i][j];
            symt(fmatrx.data(), deldip_flat.data(), ff.data());
            if (sumerr < 1e-6) break;
        }
    }
    std::vector<double> wtmass(3 * numat + 1, 0.0);
    int l = 0;
    for (int i = 1; i <= numat; ++i) {
        double weight = 1.0 / std::sqrt(atmass[i]);
        wtmass[l + 1] = weight;
        wtmass[l + 2] = weight;
        wtmass[l + 3] = weight;
        l += 3;
    }
    int linear = 0;
    for (int i = 1; i <= nvar; ++i) {
        for (int k = 1; k <= i; ++k) oldf[linear + k] = fmatrx[linear + k] * 1e5;
        for (int k = 1; k <= i; ++k)
            fmatrx[linear + k] = fmatrx[linear + k] * wtmass[i] * wtmass[k];
        linear += i;
    }
    if (!ts) frame(fmatrx, numat, 1);
    std::vector<double> cnorml_flat(nvar * nvar + 1, 0.0);
    rsp(fmatrx.data() + 1, nvar, freq.data() + 1, cnorml_flat.data() + 1);
    phase_lock(cnorml_flat, nvar);
    if (eorc && nvar == 3 * numat) symtrz(cnorml_flat.data(), freq.data(), 2, 1);
    for (int i = 1; i <= nvar; ++i) {
        int j = (int)((freq[i] + 50.0) * 0.01);
        freq[i] = freq[i] - (double)j * 100.0;
    }
    for (int i = 1; i <= nvar; ++i) freq[i] *= 1e5;
    double const_ = std::sqrt(2.0 * funcon_C::fpc_6 * funcon_C::fpc_8 * 1e11);
    for (int i = 1; i <= nvar; ++i) {
        int ii = (i - 1) * nvar;
        double summ = 0.0;
        for (int j = 1; j <= nvar / 3; ++j) {
            double a = cnorml_flat[ii + j * 3 - 2];
            double b = cnorml_flat[ii + j * 3 - 1];
            double c = cnorml_flat[ii + j * 3];
            summ += std::pow(a * a + b * b + c * c, 2) * atmass[j];
        }
        double sum = 0.0;
        for (int j = 1; j <= nvar; ++j) {
            int jii = j + ii;
            int jj = (j * (j - 1)) / 2;
            for (int k = 1; k <= j; ++k)
                sum += cnorml_flat[jii] * oldf[jj + k] * cnorml_flat[k + ii];
            for (int k = j + 1; k <= nvar; ++k)
                sum += cnorml_flat[jii] * oldf[(k * (k - 1)) / 2 + j] * cnorml_flat[k + ii];
        }
        sum *= 0.5;
        double sum1 = sum * 2.0;
        if (std::abs(freq[i]) > std::abs(sum) * 1e-20) sum = 1.0 * sum / freq[i];
        else sum = 0.0;
        to_screen_C::redmas[i][0] = summ;
        if (std::abs(freq[i]) > std::abs(sum) * 1e-20) {
            to_screen_C::redmas[i][1] = std::abs(sum1 / freq[i]);
            sum = sum / freq[i];
        } else {
            to_screen_C::redmas[i][1] = 0.0;
            sum = 0.0;
        }
        freq[i] = (freq[i] >= 0 ? 1 : -1) * std::sqrt(fact * std::abs(freq[i])) * c2pi;
        if (std::abs(freq[i]) < std::abs(sum1) * 1e20)
            sum1 = std::sqrt(std::abs(freq[i] / (sum1 * 1e-5)));
        else
            sum1 = 0.0;
        if (sum < 0 || sum > 100) sum = 0.0;
        travel[i] = sum1 * const_;
        if (travel[i] > 1.0) travel[i] = 0.0;
    }
    if (eorc) {
        int ij = 0;
        for (int i = 1; i <= nvar; ++i) {
            double sum = 0.0;
            int j = 0;
            for (int jj = 1; jj <= nvar / 3; ++jj) {
                double sum1 = 0.0;
                cnorml_flat[ij + 1] *= wtmass[j + 1];
                sum1 += cnorml_flat[ij + 1] * cnorml_flat[ij + 1];
                cnorml_flat[ij + 2] *= wtmass[j + 2];
                sum1 += cnorml_flat[ij + 2] * cnorml_flat[ij + 2];
                cnorml_flat[ij + 3] *= wtmass[j + 3];
                sum1 += cnorml_flat[ij + 3] * cnorml_flat[ij + 3];
                j += 3;
                ij += 3;
                sum += std::sqrt(sum1);
            }
            sum = 1.0 / sum;
            ij -= nvar;
            for (int k = 1; k <= nvar; ++k) cnorml_flat[ij + k] *= sum;
            ij += nvar;
        }
        for (int k = 1; k <= linear; ++k) fmatrx[k] = oldf[k] * 1e-5;
    } else {
        linear = 0;
        for (int i = 1; i <= nvar; ++i) {
            for (int k = 1; k <= i; ++k)
                fmatrx[linear + k] = oldf[linear + k] * 1e-5 * wtmass[i] * wtmass[k];
            linear += i;
        }
    }
}

