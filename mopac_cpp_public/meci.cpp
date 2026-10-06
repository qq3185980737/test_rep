// meci.cpp — full C++ translation of the meci.F90 driver (MECI).
//
// CI orchestrator: parse keywords, build active-space occupancy, transform
// two-electron integrals, enumerate microstates, build/diagonalize the CI
// matrix, and select the requested root. Print/analysis blocks (ciosci,
// dmecip, spin-density tables, ESR output) are retained as comments.

#include "meci.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "diagi.h"
#include "ijkl.h"
#include "meci_C.h"
#include "mecih.h"
#include "mndod_C.h"
using mndod_C::fx;
#include "molkst_C.h"
#include "mopend.h"
#include "perm.h"
#include "reada.h"
#include "rsp.h"
#include "symmetry_C.h"

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;
using namespace symmetry_C;
using symmetry_C::namo;
using symmetry_C::jndex;

// F90 symmetry_C members used by MECI; symtrz/matout blocks commented out.
const char* const tspin_tbl[10] = {
    "SINGLET ", "DOUBLET ", "TRIPLET ", "QUARTET ", "QUINTET ",
    "SEXTET  ", "SEPTET  ", "OCTET   ", "NONET   ", "???????"};

// Column-major 4-D xy accessor with F90 1-based semantics.
static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double meci() {
    // Cleanup branch: numat == 0 deallocate all CI state.
    if (numat == 0) {
        spin.clear();
        nalmat.clear();
        eig.clear();
        conf.clear();
        vectci.clear();
        ispin.clear();
        ispqr.clear();
        occa.clear();
        deltap.clear();
        nfa.clear();
        rjkaa.clear();
        rjkab.clear();
        eiga.clear();
        cimat.clear();
        return 0.0;
    }

    static int icalcn = 0;
    static bool first1 = true;
    static bool geook = false, lspin1 = false, peci = false, cis = false,
                cisd = false, cisdt = false, prnt2 = false;
    static bool sing = false, doub = false, trip = false, quar = false,
                quin = false, sext = false, sept = false, octe = false,
                none = false;
    static int ndoubl = 99, lroot = 1, smult = -1, ne = 0;
    static int limci = 0, nupp = 0, ndown = 0, lima = 0, limb = 0;
    static double xx = 0.0;
    static std::string root_ir = "XXXX";
    static std::vector<std::vector<int>> nperma, npermb;

    if (icalcn != numcal) {
        icalcn = numcal;
        first1 = true;
        limci = 0;
        int mdim = maxci;
        geook = keywrd.find("GEO-OK") != std::string::npos;
        lspin1 = keywrd.find("ESR") != std::string::npos;
        bool debug = keywrd.find("DEBUG") != std::string::npos;
        prnt2 = keywrd.find("MECI") != std::string::npos;
        debug = debug && prnt2;
        bool large = keywrd.find("LARGE") != std::string::npos;
        peci = keywrd.find(" PECI") != std::string::npos;
        cis = keywrd.find(" CIS ") != std::string::npos;
        cisd = keywrd.find(" CISD ") != std::string::npos;
        cisdt = keywrd.find(" CISDT ") != std::string::npos;
        root_ir = "XXXX";
        ndoubl = 99;
        size_t i0 = keywrd.find("C.I.=(");
        if (i0 != std::string::npos) {
            size_t j0 = keywrd.find(',', i0);
            // F90: j = index(keywrd(i:i+10),',') + i - 1 (1-based).
            ndoubl = (int)std::floor(reada(keywrd, (int)j0 + 1) + 0.5);
            nmos = (int)std::floor(
                reada(keywrd, (int)i0 + 6) + 0.5);  // F90: index('C.I.=(')+5
            ndoubl = std::min(ndoubl, nclose);
        } else if (keywrd.find("C.I.=") != std::string::npos) {
            nmos = (int)std::floor(
                reada(keywrd, (int)keywrd.find("C.I.=") + 5) + 0.5);  // F90: index('C.I.=')+5
        } else {
            nmos = nopen - nclose;
        }
        nmos = std::min(nmos, norbs);
        occa.assign(nmos + 1, 0.0);
        deltap.assign(nmos + 1, std::vector<double>(nmos + 1, 0.0));
        lroot = 1;
        if (keywrd.find("EXCI") != std::string::npos) lroot = 2;
        i0 = keywrd.find(" ROOT");
        if (i0 != std::string::npos) {
            size_t j0 = i0 + 6;
            for (; j0 < i0 + 10; ++j0) {
                if (keywrd[j0] > '9' || keywrd[j0] < '0') break;
            }
            if (j0 > keywrd.size()) j0 = keywrd.size();
            lroot = (int)std::floor(
                reada(keywrd.substr(0, j0), (int)i0 + 7) + 0.5);
            for (j0 = i0 + 7; j0 < i0 + 15 && j0 < keywrd.size(); ++j0) {
                if (keywrd[j0] == ' ') break;
                if (keywrd[j0] >= 'A' && keywrd[j0] <= 'Z') {
                    root_ir = keywrd.substr(j0);
                    size_t ksp = root_ir.find(' ');
                    if (ksp != std::string::npos) root_ir = root_ir.substr(0, ksp);
                    for (int kk = 2; kk <= 3 && kk < (int)root_ir.size(); ++kk) {
                        if (root_ir[kk] >= 'A' && root_ir[kk] <= 'Z')
                            root_ir[kk] = root_ir[kk] - 'A' + 'a';
                    }
                }
                if (root_ir != " ") break;
            }
        }
        int j;
        if (ndoubl == 99) {
            j = std::max(std::min((nclose + nopen + 1) / 2 - (nmos - 1) / 2,
                                  norbs - nmos + 1),
                         1);
        } else {
            j = nclose - ndoubl + 1;
            if (fract > 1.99) j = j + 1;
        }
        int l = 0;
        if (nclose - j + 1 > 0) {
            for (int p = 1; p <= nclose - j + 1; ++p) occa[p] = 1.0;
            l = nclose - j + 1;
        }
        j = std::max(j, nclose + 1);
        if (nopen - j + 1 > 0) {
            for (int p = 1; p <= nopen - j + 1; ++p) occa[l + p] = fract * 0.5;
            l = nopen - j + 1 + l;
        }
        for (int i = 1; i <= nmeci; ++i) {
            l = l + 1;
            if (l > nmos) break;
            occa[l] = 0.0;
        }
        fprintf(stderr, "[MECI] nclose=%d nopen=%d nmos=%d fract=%.4f ndoubl=%d nmeci=%d occa:", nclose, nopen, nmos, fract, ndoubl, nmeci);
        for (int ii = 1; ii <= nmos; ++ii) fprintf(stderr, " %+.4f", occa[ii]);
        fprintf(stderr, "\n"); fflush(stderr);
        sing = keywrd.find(" SING") != std::string::npos ||
               keywrd.find(" EXCI") != std::string::npos ||
               keywrd.find(" BIRAD") != std::string::npos;
        doub = keywrd.find(" DOUB") != std::string::npos;
        trip = keywrd.find(" TRIP") != std::string::npos;
        quar = keywrd.find(" QUAR") != std::string::npos;
        quin = keywrd.find(" QUIN") != std::string::npos;
        sext = keywrd.find(" SEXT") != std::string::npos;
        sept = keywrd.find(" SEPT") != std::string::npos;
        octe = keywrd.find(" OCTE") != std::string::npos;
        none = keywrd.find(" NONE") != std::string::npos;
        smult = -1;
        if (sing) smult = 1;
        if (doub) smult = 2;
        if (trip) smult = 3;
        if (quar) smult = 4;
        if (quin) smult = 5;
        if (sext) smult = 6;
        if (sept) smult = 7;
        double x = 0.0;
        for (j = 1; j <= nmos; ++j) x += occa[j];
        xx = x + x;
        ne = (int)std::floor(xx + 0.5);
        nelec = (nelecs - ne + 1) / 2;
    }
    root_requested = 0;

    if (!geook && nelec > 0) {
        if (std::fabs(eigs[nelec + 1] - eigs[nelec]) < 2e-2 ||
            (nelec + 1 + nmos < norbs &&
             std::fabs(eigs[nelec + 1 + nmos] - eigs[nelec + nmos]) < 2e-2)) {
            mopend("JOB STOPPED. TO CONTINUE, SPECIFY GEO-OK.");
            return 0.0;
        }
    }

    // eiga = active-space eigenvalues.
    eiga.assign(nmos + 1, 0.0);
    for (int i = 1; i <= nmos; ++i) eiga[i] = eigs[nelec + i];
    {
        fprintf(stderr, "[MECI] eigs1..%d:", norbs);
        for (int ii = 1; ii <= norbs; ++ii) fprintf(stderr, " %+.6f", eigs[ii]);
        fprintf(stderr, " nelec=%d\n", nelec); fflush(stderr);
    }

    // Pack MO coefficients into contiguous column-major buffers.
    std::vector<double> cp(norbs * nmos, 0.0);
    for (int col = 1; col <= nmos; ++col)
        for (int row = 1; row <= norbs; ++row)
            cp[(col - 1) * norbs + (row - 1)] = c[row][nelec + col];
    std::vector<double> cfull(norbs * norbs, 0.0);
    for (int col = 1; col <= norbs; ++col)
        for (int row = 1; row <= norbs; ++row)
            cfull[(col - 1) * norbs + (row - 1)] = c[row][col];

    // Two-electron integral transformation.
    dijkl.assign(norbs * nmos * (nmos * (nmos + 1) / 2), 0.0);
    xy.assign(nmos * nmos * nmos * nmos, 0.0);
    std::vector<double> cij(lm61), ckl(lm61), wcij(lm61);
    {
        fprintf(stderr, "[MECI] w0..10:");
        for (int ii = 0; ii < 11; ++ii) fprintf(stderr, " %+.6f", common_arrays_C::w[ii]);
        fprintf(stderr, "\n"); fflush(stderr);
    }
    ijkl(cp.data(), cfull.data(), nelec, nmos, dijkl.data(), cij.data(),
         ckl.data(), wcij.data(), xy.data());

    rjkaa.assign(nmos + 1, std::vector<double>(nmos + 1, 0.0));
    rjkab.assign(nmos + 1, std::vector<double>(nmos + 1, 0.0));
    for (int i = 1; i <= nmos; ++i)
        for (int j = 1; j <= nmos; ++j) {
            rjkaa[i][j] = xy[xyidx(i, i, j, j, nmos)] -
                          xy[xyidx(i, j, i, j, nmos)];
            rjkab[i][j] = xy[xyidx(i, i, j, j, nmos)];
        }
    rjkab1 = rjkab[1][1];

    // Remove inter-electronic interactions from eiga (before gse/diagi).
    for (int i = 1; i <= nmos; ++i) {
        double x = 0.0;
        for (int j = 1; j <= nmos; ++j)
            x += (rjkaa[i][j] + rjkab[i][j]) * occa[j];
        eiga[i] -= x;
    }
    for (int i = 1; i <= nmos; ++i)
        for (int j = 1; j <= nmos; ++j) rjkaa[i][j] *= 0.5;

    if (first1) {
        bool getmic = keywrd.find("MICROS") != std::string::npos;
        // getmic: microstate input from data file — not implemented (I/O
        // block retained as comment); falls through to enumeration.
        if (getmic) {
            std::printf(" MECI: MICROS keyword input not supported.\n");
            mopend("MICROSTATES SPECIFIED BY KEYWORDS BUT MISSING FROM DATA");
            return 0.0;
        }
        if (molkst_C::msdel == 0 && ne % 2 == 1) molkst_C::msdel = 1;
        nupp = (ne + molkst_C::msdel) / 2;
        ndown = ne - nupp;
        if (nupp * ndown < 0) {
            mopend("IMPOSSIBLE VALUE OF DELTA S");
            return 0.0;
        }
        if (limci == 0) {
            lima = (int)std::floor(
                fx[nmos + 1] / (fx[nupp + 1] * fx[nmos - nupp + 1]) + 0.5);
            limb = (int)std::floor(
                fx[nmos + 1] / (fx[ndown + 1] * fx[nmos - ndown + 1]) + 0.5);
            nperma.assign(nmos + 1, std::vector<int>(lima + 1, 0));
            npermb.assign(nmos + 1, std::vector<int>(limb + 1, 0));
        }
        lab = lima * limb;
        if (peci) limci = 4;
        if (cis) limci = 2;
        if (cisd) limci = 4;
        if (cisdt) limci = 6;
        int lima_out = lima, limb_out = limb;
        // perm expects a flat column-major (nmos, nperms) buffer.
        std::vector<int> pa(nmos * (lima + 1) + 1, 0), pb(nmos * (limb + 1) + 1, 0);
        auto npm = [&](const std::vector<int>& b, int i, int p, int nperm) {
            return b[(p - 1) * nmos + (i - 1)];
        };
        perm(pa.data(), nupp, nmos, lima_out, limci);
        perm(pb.data(), ndown, nmos, limb_out, limci);
        lima = lima_out;
        limb = limb_out;
        lab = lima * limb;
        microa.assign(nmos + 1, std::vector<int>(lab + 1, 0));
        microb.assign(nmos + 1, std::vector<int>(lab + 1, 0));
        // Keep flat permutations for later (limci==0 path uses nperma/npermb).
        nperma.assign(nmos + 1, std::vector<int>(lima + 1, 0));
        npermb.assign(nmos + 1, std::vector<int>(limb + 1, 0));
        for (int p = 1; p <= lima; ++p)
            for (int i = 1; i <= nmos; ++i) nperma[i][p] = npm(pa, i, p, lima);
        for (int p = 1; p <= limb; ++p)
            for (int i = 1; i <= nmos; ++i) npermb[i][p] = npm(pb, i, p, limb);
    } else {
        if (limci != 0) {
            // microa/microb already set on the first pass; keep sizes.
        }
        lab = lima * limb;
    }

    spin.assign(lab + 1, 0.0);
    nalmat.assign(lab + 1, 0);
    eig.assign(lab + 2, 0.0);
    vectci.assign(30 * lab, 0.0);
    ispin.assign(lab + 1, 0);
    ispqr.assign(lab + 1, std::vector<int>(nmeci + 2, -99999));

    // Ground-state energy of the (fractional) reference configuration.
    double gse = 0.0;
    for (int i = 1; i <= nmos; ++i) {
        gse += eiga[i] * occa[i] * 2.0;
        gse += xy[xyidx(i, i, i, i, nmos)] * occa[i] * occa[i];
        for (int j = i + 1; j <= nmos; ++j)
            gse += 2.0 * (2.0 * xy[xyidx(i, i, j, j, nmos)] -
                          xy[xyidx(i, j, i, j, nmos)]) *
                   occa[i] * occa[j];
    }

    // Microstate enumeration / excitation filter.
    if (limci != 0) {
        int loc = 0;
        for (int i = 1; i <= lima; ++i)
            for (int j = 1; j <= limb; ++j) {
                int k = 0;
                for (int l = 1; l <= nmos; ++l)
                    k += std::abs(nperma[l][i] - nperma[l][1]) +
                         std::abs(npermb[l][j] - npermb[l][1]);
                if (peci) {
                    for (int l = 1; l <= nmos; ++l)
                        if (k > 2 && nperma[l][i] != npermb[l][j]) k = 1000;
                }
                if (k > limci) continue;
                ++loc;
                for (int l = 1; l <= nmos; ++l) {
                    microa[l][loc] = nperma[l][i];
                    microb[l][loc] = npermb[l][j];
                }
            }
        lab = loc;
    }

    if (lab > maxci) {
        mopend("Too many configurations requested");
        return 0.0;
    }
    conf.assign(lab * lab, 0.0);
    cimat.assign((lab * (lab + 1)) / 2 + 9, 0.0);
    std::vector<double> diag(lab + 1, 0.0);
    cdiag.assign(lab + 1, 0.0);

    std::vector<int> ma(nmos + 1, 0), mb(nmos + 1, 0);
    auto pack_col = [&](int ic) {
        for (int j = 1; j <= nmos; ++j) {
            ma[j] = microa[j][ic];
            mb[j] = microb[j][ic];
        }
    };
    if (limci != 0) {
        for (int i = 1; i <= lab; ++i) {
            pack_col(i);
            diag[i] = diagi(ma.data(), mb.data(), eiga.data(), xy.data(),
                            nmos) -
                      gse;
            cdiag[i] = cdiagi;
        }
    } else {
        // Full CI: diag from every (alpha,beta) permutation pair.
        std::vector<double> cimat2(lima * limb + 1, 0.0);
        int ii = 0;
        std::vector<int> ma2(nmos + 1, 0), mb2(nmos + 1, 0);
        for (int i1 = 1; i1 <= lima; ++i1)
            for (int i2 = 1; i2 <= limb; ++i2) {
                ++ii;
                for (int j = 1; j <= nmos; ++j) {
                    ma2[j] = nperma[j][i1];
                    mb2[j] = npermb[j][i2];
                }
                cimat2[ii] = diagi(ma2.data(), mb2.data(), eiga.data(),
                                   xy.data(), nmos) -
                             gse;
            }
        int mdim = maxci;
        if (lab <= mdim) {
            int k = 0;
            for (int i = 1; i <= lima; ++i)
                for (int j = 1; j <= limb; ++j) {
                    ++k;
                    diag[k] = cimat2[k];
                    for (int m = 1; m <= nmos; ++m) {
                        microa[m][k] = nperma[m][i];
                        microb[m][k] = npermb[m][j];
                    }
                }
        } else {
            for (int labi = 1; labi <= mdim; ++labi) {
                double x = 1000.0;
                int jj = 0, j1 = 0, j2 = 0, ii2 = 0;
                for (int i1 = 1; i1 <= lima; ++i1)
                    for (int i2 = 1; i2 <= limb; ++i2) {
                        ++ii2;
                        if (cimat2[ii2] >= x) continue;
                        x = cimat2[ii2];
                        jj = ii2;
                        j1 = i1;
                        j2 = i2;
                    }
                for (int m = 1; m <= nmos; ++m) {
                    microa[m][labi] = nperma[m][j1];
                    microb[m][labi] = npermb[m][j2];
                }
                diag[labi] = cimat2[jj];
                cimat2[jj] = 1e8;
            }
            lab = mdim;
        }
    }

    // nalmat (alpha occupancy) and spin quantum number.
    for (int i = 1; i <= lab; ++i) {
        int k = 0;
        double x = 0.0;
        for (int j = 1; j <= nmos; ++j) {
            x += microa[j][i] * microb[j][i];
            k += microa[j][i];
        }
        nalmat[i] = k;
        double dm = xx - 2 * nalmat[i];
        spin[i] = 4.0 * x - dm * dm;
    }

    if (lab < lroot) {
        mopend("C.I. IS OF SIZE LESS THAN ROOT SPECIFIED. MODIFY SIZE OF "
               "C.I. OR ROOT NUMBER.");
        return 0.0;
    }

    // Build the CI matrix and diagonalize.
    mecih(diag.data(), cimat.data(), nmos, &lab, xy.data());

    // Perturb diagonal to destroy exact degeneracies.
    for (int i = 1; i <= lab; ++i)
        cimat[(i * (i + 1)) / 2] +=
            1e-10 * i * (i % 2 == 0 ? 1 : -1);

    {
        fprintf(stderr, "[MECI] diag1..%d:", lab);
        for (int ii = 1; ii <= lab; ++ii) fprintf(stderr, " %+.6f", diag[ii]);
        fprintf(stderr, " gse=%.6f\n", gse); fflush(stderr);
        if (lab >= 2) {
            fprintf(stderr, "[MECI] xy1211=%.6f xy1121=%.6f xy1111=%.6f xy1122=%.6f xy1212=%.6f xy2222=%.6f\n",
                    xy[xyidx(1,2,1,1,nmos)], xy[xyidx(1,1,2,1,nmos)],
                    xy[xyidx(1,1,1,1,nmos)], xy[xyidx(1,1,2,2,nmos)],
                    xy[xyidx(1,2,1,2,nmos)], xy[xyidx(2,2,2,2,nmos)]); fflush(stderr);
            fprintf(stderr, "[MECI] eiga1..2=%.6f %.6f nelec=%d ne=%d\n", eiga[1], eiga[2], nelec, ne); fflush(stderr);
        }
    }

    labsiz = std::min(lab, lroot + 10);
    if ((last == 3 && prnt2) || (root_requested != 1 && first1) ||
        root_ir != "XXXX")
        labsiz = lab;

    rsp(cimat.data(), lab, eig.data(), conf.data());
    {
        fprintf(stderr, "[MECI] conf16:"); for (int ii = 0; ii < 16; ++ii) fprintf(stderr, " %+.6f", conf[ii]); fprintf(stderr, "\n"); fflush(stderr);
        fprintf(stderr, "[MECI] eig4:"); for (int ii = 0; ii < 4; ++ii) fprintf(stderr, " %+.8f", eig[ii]); fprintf(stderr, "\n"); fflush(stderr);
        for (int cc = 0; cc < lab; ++cc) {
            double n2 = 0.0;
            for (int rr = 0; rr < lab; ++rr) n2 += conf[cc * lab + rr] * conf[cc * lab + rr];
            fprintf(stderr, "[MECI] col%d norm=%.8f\n", cc + 1, std::sqrt(n2)); fflush(stderr);
        }
    }

    // Ensure symmetry arrays are sized.
    if ((int)namo.size() < lab + 1) namo.resize(lab + 1);
    if ((int)jndex.size() < lab + 1) jndex.resize(lab + 1);

    // State analysis: spin quantum numbers.
    int iroot = 0;
    int jroot = 1;
    std::fill(cimat.begin(), cimat.end(), 0.1);
    std::fill(ispin.begin() + 1, ispin.begin() + lab + 1, 0);
    for (int i = 1; i <= labsiz; ++i) {
        if (last == 0 && root_ir == "XXXX") {
            namo[i] = " ";
            jndex[i] = i;
        }
        double x = 0.5 * xx;
        int ii = (i - 1) * lab;
        for (int j = 1; j <= lab; ++j) {
            int ji = ii + (j - 1);
            x -= conf[ji] * conf[ji] * spin[j] * 0.25;
            int k = ispqr[j][1];
            if (k == 1) continue;
            for (int l = 1; l <= k - 1; ++l)
                x += conf[ji] * conf[ii + ispqr[j][l + 1] - 1] * 2.0;
        }
        double y = (-1.0 + std::sqrt(1.0 + 4.0 * x)) * 0.5;
        ispin[i] = (int)std::floor(y * 2.0 + 1.0 + 0.5);
        if (std::fabs(ispin[i] - y * 2.0 - 1.0) > 0.2)
            ispin[i] = 10;
        else
            ispin[i] = std::min(ispin[i], 10);
        cimat[lab + 1] += 1.0;
    }

    // Assign quantum numbers, group by irrep and spin.
    int qn_temp[21][21];
    for (int a = 0; a < 21; ++a)
        for (int b = 0; b < 21; ++b) qn_temp[a][b] = 0;
    std::vector<std::string> srep(21, " ");
    int irep[21];
    for (int a = 0; a < 21; ++a) irep[a] = 0;
    int isize = std::max(maxci, std::max(nmos + 3, lab));
    nfa.assign(isize + 1, 0);
    for (int i = 1; i <= lab; ++i) nfa[i] = jndex[i];
    int jcnt = 0, lcnt = 0;
    for (int i = 1; i <= lab; ++i) {
        if (i > 1) {
            if (nfa[i] == nfa[i - 1] && namo[i] == namo[i - 1]) {
                jndex[i] = jndex[i - 1];
                if (smult < 0 || ispin[i] == smult) ++iroot;
                continue;
            }
        }
        int k = 1;
        for (; k <= jcnt; ++k)
            if (srep[k] == namo[i]) break;
        if (k > jcnt) {
            ++jcnt;
            srep[jcnt] = namo[i];
        }
        int m = 1;
        for (; m <= lcnt; ++m)
            if (irep[m] == ispin[i]) break;
        if (m > lcnt) {
            ++lcnt;
            irep[lcnt] = ispin[i];
        }
        if (k > 20 || m > 20) {
            mopend("Dimension error in MECI.F90");
            return 0.0;
        }
        qn_temp[k][m]++;
        jndex[i] = qn_temp[k][m];
        if (root_requested < 1 && (smult < 0 || ispin[i] == smult)) {
            ++iroot;
            if ((iroot >= lroot && jroot <= lroot && root_ir == "XXXX") ||
                (lroot == jndex[i] && root_ir == namo[i])) {
                root_requested = i;
                double meci_val = eig[i - 1];
                int mcnt = 0;
                int kk = i;
                for (; kk <= std::min(lab, i + 22); ++kk) {
                    if (std::fabs(eig[kk - 1] - eig[i - 1]) > 1e-2) break;
                    for (int p = 0; p < lab; ++p)
                        vectci[mcnt + p + 1] = conf[lab * (kk - 1) + p];
                    mcnt += lab;
                }
                if (kk > std::min(lab, i + 22)) kk = std::min(lab, i + 4) + 1;
                nstate = kk - i;
                (void)meci_val;
                fprintf(stderr, "[MECI] root_requested=%d nstate=%d vectci0..7=%+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
                        i, nstate, vectci[0], vectci[1], vectci[2], vectci[3],
                        vectci[4], vectci[5], vectci[6], vectci[7]); fflush(stderr);
            }
            jroot = iroot;
        }
    }
    if (root_requested < 1) {
        mopend("ROOT REQUESTED DOES NOT EXIST IN C.I.");
        return 0.0;
    }
    state_Irred_Rep = namo[root_requested];
    state_spin = tspin_tbl[ispin[root_requested] - 1];
    state_QN = jndex[root_requested];

    first1 = false;
    return eig[root_requested - 1];
}
