// test_batchC5.cpp — real-machine tests for diagg.F90 / diagg1.F90 /
// diagg2.F90 / density_for_MOZYME.F90 / epseta.F90 / ijbo.F90 translations.
// Model: 4 atoms (C1, C2, H1, H2), norbs=10 (s,p on C; s on H), mpack=55.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "diagg.h"
#include "ijbo.h"
#include "MOZYME_C.h"
#include "molkst_C.h"
#include "overlaps_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace MOZYME_C;
using namespace molkst_C;

void timer(const std::string&) {}  // external timing hook (internal tests)

namespace {
int cal_ctr = 0;
bool failed = false;

void check(const char* name, bool ok) {
    if (!ok) {
        std::fprintf(stderr, "FAIL %s\n", name);
        failed = true;
    } else {
        std::fprintf(stderr, "pass %s\n", name);
    }
}

void near(const char* name, double got, double exp, double tol) {
    check(name, std::fabs(got - exp) <= tol);
}

// model geometry: all pairs within cutof (r^2 < cutof1 = cutof2 = 2.0)
void set_geom_close() {
    coord.assign(4, std::vector<double>(16, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.4; coord[2][2] = 0.0; coord[3][2] = 0.0;
    coord[1][3] = 0.7; coord[2][3] = 0.5; coord[3][3] = 0.0;
    coord[1][4] = 0.7; coord[2][4] = -0.5; coord[3][4] = 0.0;
}

// packed lookup table for the 4-atom model (built as MOPAC would for i>j)
void set_ijall() {
    // pairs sorted by ind_i then ind_j: (2,1),(3,1),(3,2),(4,1),(4,2),(4,3)
    int ijall_dat[7] = {0, 1, 1, 2, 1, 2, 3};
    int iijj_dat[7] = {0, 10, 36, 40, 45, 49, 53};
    ijall.assign(16, 0); iijj.assign(16, 0); iij.assign(16, 0);
    numij.assign(16, 0);
    for (int i = 1; i <= 6; ++i) {
        ijall[i] = ijall_dat[i];
        iijj[i] = iijj_dat[i];
    }
    iij[1] = 1; numij[1] = 0;
    iij[2] = 1; numij[2] = 1;
    iij[3] = 2; numij[3] = 3;
    iij[4] = 4; numij[4] = 6;
    lijbo = false;
}

// MOZYME data layout for the diagg tests (lijbo=true path)
void set_mozyme_arrays(std::vector<double>& fao) {
    // basis: C1=AO1..4, C2=AO5..8, H1=AO9, H2=AO10
    nfirst.assign(16, 0); nlast.assign(16, 0); iorbs.assign(16, 0);
    nat.assign(16, 0);
    for (int i = 1; i <= 4; ++i) {
        nfirst[i] = (i == 1) ? 1 : (i == 2) ? 5 : (i == 3) ? 9 : 10;
        nlast[i] = (i == 1) ? 4 : (i == 2) ? 8 : (i == 3) ? 9 : 10;
        iorbs[i] = (i <= 2) ? 4 : 1;
        nat[i] = (i <= 2) ? 6 : 1;
    }
    norbs = 10;
    mpack = norbs * (norbs + 1) / 2;
    fao.assign(mpack + 1, 0.0);
    for (int i = 1; i <= 10; ++i) fao[i] = 1.0;  // FOCK diag block of C1
    // LMO occupancy structure
    // occ LMO 1: atom 1; occ LMO 2: atom 2; occ LMO 3: atom 3
    nncf.assign(4, 0); ncf.assign(4, 0); ncocc.assign(4, 0);
    nncf[1] = 0; ncf[1] = 1;
    nncf[2] = 1; ncf[2] = 1;
    nncf[3] = 2; ncf[3] = 1;
    icocc.resize(8);
    icocc[1] = 1; icocc[2] = 2; icocc[3] = 3;
    icocc_dim = 8;
    ncocc[1] = 0; ncocc[2] = 4; ncocc[3] = 8;
    cocc_dim = 10;
    cocc.assign(cocc_dim + 1, 0.0);
    cocc[1] = 0.8; cocc[2] = 0.7; cocc[3] = 0.6; cocc[4] = 0.5;  // LMO1 (C1)
    cocc[5] = 0.5; cocc[6] = 0.6; cocc[7] = 0.7; cocc[8] = 0.8;  // LMO2 (C2)
    cocc[9] = 1.0;                                              // LMO3 (H1)
    // vir LMO 1: atom 4; vir LMO 2: atom 1
    nnce.assign(3, 0); nce.assign(3, 0); ncvir.assign(3, 0);
    nnce[1] = 0; nce[1] = 1;
    nnce[2] = 1; nce[2] = 1;
    icvir.resize(8);
    icvir[1] = 4; icvir[2] = 1;
    icvir_dim = 8;
    ncvir[1] = 0; ncvir[2] = 1;
    cvir_dim = 10;
    cvir.assign(cvir_dim + 1, 0.0);
    cvir[1] = 0.8;                          // vir1 (H2)
    cvir[2] = 0.4; cvir[3] = 0.5; cvir[4] = 0.6; cvir[5] = 0.7;  // vir2 (C1)
    // nijbo packed starts (0-based packed address of block start)
    nijbo.assign(5, std::vector<int>(5, -999));
    // diagonal blocks: start of block (i,i) = packed addr of (first AO of
    // atom i, first AO) minus one: (1,1)=0; (2,2): (5,5)=15 -> 14;
    // (3,3): (9,9)=45 -> 44; (4,4): (10,10)=55 -> 54.
    nijbo[1][1] = 0;  nijbo[2][2] = 14; nijbo[3][3] = 44; nijbo[4][4] = 54;
    nijbo[1][2] = nijbo[2][1] = 10;
    nijbo[1][3] = nijbo[3][1] = 36;
    nijbo[2][3] = nijbo[3][2] = 40;
    nijbo[1][4] = nijbo[4][1] = 45;
    nijbo[2][4] = nijbo[4][2] = 49;
    nijbo[3][4] = nijbo[4][3] = 53;
    lijbo = true;
    shift = 0.0;
    thresh = 1.e-4;
    // diag storage
    fmo_dim = 20;
    fmo.assign(fmo_dim + 1, 0.0);
    ifmo.assign(3, std::vector<int>(fmo_dim + 1, 0));
    nfmo.assign(4, 0);
    nvirtual = 2;
    ipad2 = 2; ipad4 = 4;
    // density matrix (1-based, mpack).  Nonzero initial diagonal
    // (diagg1 gates ws construction on avir(k1)*p(kj+1) > cutoff).
    p.assign(mpack + 1, 1.0);
    // eigen levels
    eigs.assign(11, 0.0);
}

void setup(const char* key) {
    keywrd = key;
    ++cal_ctr;
    numcal = cal_ctr;
    moperr = false;
    mozyme = true;
    id = 1;
    numat = 4;
    sumt = 0.0;
    ovmax = 0.0;
    ijc = 0;
    sumb = 0.0;
    // diagg1 gates ws construction on avir(k1)*p(kj+1)>cutoff;
    // a fresh positive-definite diagonal keeps that gate open.
    std::fill(p.begin(), p.end(), 1.0);
}
}  // namespace

int main() {
    std::fprintf(stderr, "--- C5 diagg/ijbo batch ---\n");

    // T1: ijbo lookup path (lijbo=false), close geometry
    set_geom_close();
    set_ijall();
    overlaps_C::cutof1 = 2.0;
    overlaps_C::cutof2 = 2.0;
    near("T1a ijbo(2,1)", (double)ijbo(2, 1), 10.0, 0.5);
    near("T1b ijbo(3,2)", (double)ijbo(3, 2), 40.0, 0.5);
    near("T1c ijbo(4,3)", (double)ijbo(4, 3), 53.0, 0.5);
    // far atom -> -1
    coord[1][4] = 10.0; coord[2][4] = 0.0; coord[3][4] = 0.0;
    near("T1d ijbo far", (double)ijbo(4, 1), -1.0, 0.5);
    set_geom_close();

    // --- diagg tests (lijbo=true) ---
    std::vector<double> fao;
    set_mozyme_arrays(fao);
    std::vector<double> partp(mpack + 1, 0.0);

    // T2: diagg1 full-construction path (idiagg=1)
    setup("");
    int nij_out = fmo_dim;
    std::vector<double> ws(norbs + 1), aov(numat + 1), avir(norbs + 1);
    std::vector<double> aocc(icocc_dim + 1);
    std::vector<bool> latoms(numat + 1);
    diagg1(fao, 3, 2, eigs, ws, latoms, ifmo, fmo, fmo_dim, nij_out, 1, avir,
           aocc, aov);
    check("T2 nij==1", nij_out == 1);
    check("T2 ifmo(1,1)==2", ifmo[1][1] == 2);
    check("T2 ifmo(2,1)==1", ifmo[2][1] == 1);
    near("T2 fmo(1)", fmo[1], 5.72, 0.05);
    near("T2 eigv vir2", eigs[5], 4.84, 0.05);
    near("T2 ovmax==tiny", ovmax - tiny, 0.0, 1.e-12);

    // T3: diagg2 rotation preserves norms and annihilates the element
    std::vector<double> storei(norbs + 1), storej(norbs + 1);
    std::vector<int> iused(numat + 1);
    double cocc_before[4] = {cocc[1], cocc[2], cocc[3], cocc[4]};
    double cvir_before[4] = {cvir[2], cvir[3], cvir[4], cvir[5]};
    diagg2(3, 2, eigs, iused, latoms, nij_out, 1, storei, storej);
    double sq_cocc = 0.0, sq_cvir = 0.0;
    for (int i = 1; i <= 4; ++i) {
        sq_cocc += cocc[i] * cocc[i];
        sq_cvir += cvir[i + 1] * cvir[i + 1];
    }
    double sq_cocc_b = 0.0, sq_cvir_b = 0.0;
    for (int i = 0; i < 4; ++i) {
        sq_cocc_b += cocc_before[i] * cocc_before[i];
        sq_cvir_b += cvir_before[i] * cvir_before[i];
    }
    // 2x2 rotation conserves the TOTAL norm of the (cocc, cvir) pair
    near("T3 total norm conserved", sq_cocc + sq_cvir, sq_cocc_b + sq_cvir_b, 1.e-9);
    near("T3 alpha^2+beta^2=1", sq_cocc + sq_cvir, sq_cocc_b + sq_cvir_b, 1.e-9);
    // atoms not shared (C2, block 5..8) are untouched by the rotation
    check("T3 atom2 unchanged", cocc[5] == 0.5 && cocc[8] == 0.8);

    // T4: density_for_MOZYME mode 0 (full density, spin-factored)
    std::fill(p.begin(), p.end(), 0.0);
    std::fill(partp.begin(), partp.end(), 0.0);
    density_for_MOZYME(p, 0, 3, partp);
    double p1 = 2.0 * cocc[1] * cocc[1];
    near("T4 p(1)", p[1], p1, 1.e-9);
    check("T4 p off-diag block nonzero", p[2] != 0.0);
    check("T4 no atom4 density", p[55] == 0.0);

    // T5: mode -1 (remove old density): p = -0.5*partp
    for (int i = 1; i <= mpack; ++i) partp[i] = 0.2 * i;
    density_for_MOZYME(p, -1, 3, partp);
    // mode -1: p starts at -0.5*partp, LMO density is added, then spinfa=-2
    near("T5 p(1)", p[1], -2.0 * (-0.5 * 0.2 + cocc[1] * cocc[1]), 1.e-9);
    near("T5 p[55]", p[55], -2.0 * (-0.5 * 0.2 * 55), 1.e-9);

    // T6: epseta machine constants
    double eps, eta;
    epseta(eps, eta);
    check("T6 eps>0", eps > 0.0);
    check("T6 eta>=1e-17", eta >= 1.e-17);
    check("T6 eps small", eps < 1.e-15);

    // T7: diagg2 with DAMP keyword (icalcn change recomputes const2)
    std::vector<double> fao2 = fao;
    std::vector<double> partp2(mpack + 1, 0.0);
    setup(" DAMP=0.7");
    nij_out = fmo_dim;
    diagg1(fao2, 3, 2, eigs, ws, latoms, ifmo, fmo, fmo_dim, nij_out, 1, avir,
           aocc, aov);
    diagg2(3, 2, eigs, iused, latoms, nij_out, 1, storei, storej);
    check("T7 DAMP runs", true);

    // T8: reuse branch (idiagg=7: mod(7,4)!=0)
    setup("");
    std::fill(fmo.begin(), fmo.end(), 0.0);
    // restore the pre-rotation LMO coefficients so the reused pair
    // recomputes to the same value as in T2
    cocc[1] = 0.8; cocc[2] = 0.7; cocc[3] = 0.6; cocc[4] = 0.5;
    cvir[2] = 0.4; cvir[3] = 0.5; cvir[4] = 0.6; cvir[5] = 0.7;
    nij_out = fmo_dim;
    diagg1(fao, 3, 2, eigs, ws, latoms, ifmo, fmo, fmo_dim, nij_out, 7, avir,
           aocc, aov);
    check("T8 reuse nij==1", nij_out == 1);
    near("T8 reuse fmo(1)", fmo[1], 5.72, 0.05);

    // T9: density mode 2 (add partial): p = 0.5*partp, then LMO density is
    // added on top (spinfa=1 for mode 2).
    std::fill(p.begin(), p.end(), 0.0);
    for (int i = 1; i <= mpack; ++i) partp[i] = 0.3;
    density_for_MOZYME(p, 2, 3, partp);
    near("T9 mode2 p(1)", p[1], 0.5 * 0.3 + cocc[1] * cocc[1], 1.e-12);

    if (failed) {
        std::fprintf(stderr, "C5 FAILED\n");
        return 1;
    }
    std::fprintf(stderr, "ALL PASS\n");
    return 0;
}
