// pinout.cpp — C++ translation of MOPAC 2016 "pinout.F90".
// Writes/reads all localised-MO information to/from the density file.
#include "pinout.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
using namespace MOZYME_C;

namespace {
using molkst_C::numat;
using molkst_C::keywrd;
using molkst_C::nelecs;
using molkst_C::norbs;
using molkst_C::line;
using chanel_C::iw;
using chanel_C::iden;
using chanel_C::density_fn;
using common_arrays_C::nbonds;
using common_arrays_C::ibonds;

// read lo..hi (1-based) elements of an int vector from a binary stream.
bool read_block(std::ifstream& in, std::vector<int>& v, int hi, int lo = 1) {
    for (int i = lo; i <= hi; ++i) in.read((char*)&v[i], sizeof(int));
    return !in.fail();
}
// read lo..hi (1-based) elements of a double vector from a binary stream.
bool read_block(std::ifstream& in, std::vector<double>& v, int hi, int lo) {
    for (int i = lo; i <= hi; ++i) in.read((char*)&v[i], sizeof(double));
    return !in.fail();
}
}  // namespace

extern void prtlmo();
extern void mopend(const std::string&);
extern void to_screen(const std::string&);

// MOZYME module data (global linkage).

void pinout(int mode) {
    using MOZYME_C::iorbs;
    int nocc = nelecs / 2;
    int nvir = norbs - nocc;
    std::ifstream in;
    std::ofstream out;
    // mode == 0: verify the density file exists.
    if (mode == 0) {
        in.open(density_fn, std::ios::binary);
        if (!in.is_open()) {
            line = "         FILE " + density_fn + " IS MISSING";
            std::fprintf(stdout, "\n\n%s\n", line.c_str());
            to_screen(line);
            mopend("DENSITY file is missing");
            return;
        }
        in.close();
    }
    if (mode == 1) {
        out.open(density_fn, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) goto err1030;
        // OUTPUT ALL DATA FOR LMOs
        for (int i = 1; i <= nocc; ++i) out.write((const char*)&ncf[i], sizeof(int));
        for (int i = 1; i <= nvir; ++i) out.write((const char*)&nce[i], sizeof(int));
        for (int i = 1; i <= nocc; ++i)
            for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j)
                out.write((const char*)&icocc[j], sizeof(int));
        for (int i = 1; i <= nvir; ++i)
            for (int j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j)
                out.write((const char*)&icvir[j], sizeof(int));
        for (int i = 1; i <= numat; ++i) out.write((const char*)&iorbs[i], sizeof(int));
        for (int i = 1; i <= numat; ++i) out.write((const char*)&nbonds[i], sizeof(int));
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= 9; ++j) out.write((const char*)&ibonds[j][i], sizeof(int));
        // Now write the LMOs.
        for (int i = 1; i <= nocc; ++i) {
            int l = 0;
            for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j)
                l = l + iorbs[icocc[j]];
            for (int j = ncocc[i] + 1; j <= ncocc[i] + l; ++j)
                out.write((const char*)&cocc[j], sizeof(double));
        }
        for (int i = 1; i <= nvir; ++i) {
            int l = 0;
            for (int j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j)
                l = l + iorbs[icvir[j]];
            for (int j = ncvir[i] + 1; j <= ncvir[i] + l; ++j)
                out.write((const char*)&cvir[j], sizeof(double));
        }
        out.close();
        goto done;
    } else {
        // READ IN ALL DATA FOR LMOs
        in.open(density_fn, std::ios::binary);
        if (!in.is_open() || !read_block(in, ncf, nocc)) goto err1010;
        if (!read_block(in, nce, nvir)) goto err1020;
        int l = 0;
        for (int i = 1; i <= nocc; ++i) l = l + ncf[i];
        if ((l * 8) / 5 > icocc_dim) {
            icocc_dim = (l * 8) / 5;
            cocc_dim = (int)std::lround(cocc_dim * (l * 1.6 / (l * 8 / 5)));
            if (cocc_dim < icocc_dim) cocc_dim = icocc_dim;
            icocc.resize(icocc_dim + 1, 0);
            cocc.resize(cocc_dim + 1, 0.0);
        }
        l = 0;
        for (int i = 1; i <= nvir; ++i) l = l + nce[i];
        if ((l * 8) / 5 > icvir_dim) {
            icvir_dim = (l * 8) / 5;
            cvir_dim = (int)std::lround(cvir_dim * (l * 1.6 / (l * 8 / 5)));
            if (cvir_dim < icvir_dim) cvir_dim = icvir_dim;
            icvir.resize(icvir_dim + 1, 0);
            cvir.resize(cvir_dim + 1, 0.0);
        }
        // COMPRESS INCOMING DATA.
        int j = 0;
        for (int i = 1; i <= nocc; ++i) { nncf[i] = j; j = j + ncf[i]; }
        j = 0;
        for (int i = 1; i <= nvir; ++i) { nnce[i] = j; j = j + nce[i]; }
        for (int i = 1; i <= nocc; ++i)
            if (!read_block(in, icocc, nncf[i] + ncf[i], nncf[i] + 1)) goto err1020;
        for (int i = 1; i <= nvir; ++i)
            if (!read_block(in, icvir, nnce[i] + nce[i], nnce[i] + 1)) goto err1020;
        for (int i = 1; i <= numat; ++i) in.read((char*)&iorbs[i], sizeof(int));
        if (in.fail()) goto err1020;
        for (int i = 1; i <= numat; ++i) in.read((char*)&nbonds[i], sizeof(int));
        if (in.fail()) goto err1020;
        for (int i = 1; i <= numat; ++i)
            for (int j1 = 1; j1 <= 9; ++j1) in.read((char*)&ibonds[j1][i], sizeof(int));
        if (in.fail()) goto err1020;
        // NOW READ THE LMOs
        int k = 0;
        for (int i = 1; i <= (int)cocc.size() - 1; ++i) cocc[i] = 0.0;
        for (int i = 1; i <= nocc; ++i) {
            ncocc[i] = k;
            int l2 = 0;
            for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j)
                l2 = l2 + iorbs[icocc[j]];
            k = k + l2;
            if (!read_block(in, cocc, ncocc[i] + l2, ncocc[i] + 1)) goto err1020;
        }
        k = 0;
        for (int i = 1; i <= (int)cvir.size() - 1; ++i) cvir[i] = 0.0;
        for (int i = 1; i <= nvir; ++i) {
            ncvir[i] = k;
            int l2 = 0;
            for (int j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j)
                l2 = l2 + iorbs[icvir[j]];
            k = k + l2;
            if (!read_block(in, cvir, ncvir[i] + l2, ncvir[i] + 1)) goto err1020;
        }
        in.close();
        goto done;
    err1010:
        std::fprintf(stdout, "\n\n%20s\n", (" FILE " + density_fn + " EXISTS, BUT IS FAULTY").c_str());
        mopend("DENSITY file is faulty");
        return;
    err1020:
        std::fprintf(stdout, "\n\n%20s\n", (" OLDENS FILE FOR " + density_fn + " IS CORRUPT").c_str());
        mopend("DENSITY file is corrupt");
        return;
    }
err1030:
    line = " Cannot write density matrix to \"" + density_fn + "\"";
    mopend(line);
    return;
done:
    if (keywrd.find(" PINOUT") != std::string::npos) {
        prtlmo();
    }
}
