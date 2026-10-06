// test_batchC6.cpp — verification for batch C6 (geo_ref).
// T1 geo_diff numerics; T2 l_control 3-arg add/remove; T3 geo_ref 0SCF full
// run; T4 geo_ref with a missing file -> mopend.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "big_swap.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "geo_diff.h"
#include "geo_ref.h"
#include "maps_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

// Hooks for geo_ref (stubs until the full implementations are ported).
void dock(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, double&) {}
void add_path(std::string&) {}
void pdbout(int) {}

static int cal_ctr = 0;

static void setup() {
    ++cal_ctr;
    numcal = 100 + cal_ctr;  // avoid the ARC-detection branch in getgeo
    moperr = false;
    mozyme = false;
    keywrd = " ";
    line = " ";
    koment = " ";
    title = " ";
    ncomments = 0;
    numat = 0;
    natoms = 0;
    id = 1;
    nvar = 0;
    maxtxt = 0;
    numat_old = 0;
    ir = 5;
    iw = 6;
    output_fn = "test_dat.out";
    job_fn = "test_dat.mop";
    input_fn = "test_dat.dat";
    geo_dat_name = " ";
    geo_ref_name = " ";
    arc_hof_1 = 0.0;
    arc_hof_2 = 0.0;
    density = 0.0;
    na.assign(20, 0);
    nb.assign(20, 0);
    nc.assign(20, 0);
    nat.assign(20, 0);
    labels.assign(20, 0);
    atmass.assign(20, 0.0);
    txtatm.assign(20, " ");
    txtatm1.assign(20, " ");
    geo.assign(4, std::vector<double>(20, 0.0));
    coord.assign(4, std::vector<double>(20, 0.0));
    geoa.assign(4, std::vector<double>(20, 0.0));
    c.assign(4, std::vector<double>(20, 0.0));
    xparam.assign(80, 0.0);
    loc.assign(4, std::vector<int>(80, 0));
    simbol.assign(30, "---------");
    refkey.assign(7, " ");
    all_comments.assign(20, " ");
}

static void expect_close(const char* what, double got, double want, double tol) {
    if (std::fabs(got - want) > tol) {
        std::fprintf(stderr, "FAIL %s: got %.8f want %.8f\n", what, got, want);
        std::exit(1);
    }
}

// ---------------------------------------------------------------------------
static void T1_geo_diff() {
    setup();
    numat = 2;
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    geo[1][2] = 1.0; geo[2][2] = 1.0; geo[3][2] = 1.0;
    geoa[1][1] = 0.0; geoa[2][1] = 0.0; geoa[3][1] = 0.0;
    geoa[1][2] = 4.0; geoa[2][2] = 4.0; geoa[3][2] = 4.0;
    double sum = 0.0, rms = 0.0;
    geo_diff(sum, rms, false);
    expect_close("T1 sum", sum, std::sqrt(27.0), 1.e-9);
    expect_close("T1 rms", rms, 27.0, 1.e-9);
    std::fprintf(stderr, "T1 ok\n");
}

// ---------------------------------------------------------------------------
static void T2_l_control() {
    setup();
    keywrd = " XX     ";
    l_control("YY", 2, 1);  // add
    std::string k = keywrd;
    while (!k.empty() && k.back() == ' ') k.pop_back();
    assert(k == " XX YY");
    l_control("YY", 2, -1);  // remove
    k = keywrd;
    while (!k.empty() && k.back() == ' ') k.pop_back();
    assert(k == " XX");
    // HTML add/remove round trip (as used by geo_ref around pdbout).
    // keywrd is fixed-length with trailing blanks in Fortran; emulate that.
    keywrd = " XX HTML" + std::string(20, ' ');
    l_control("HTML", 4, -1);
    k = keywrd;
    while (!k.empty() && k.back() == ' ') k.pop_back();
    assert(k == " XX");
    l_control("HTML", 4, 1);
    k = keywrd;
    while (!k.empty() && k.back() == ' ') k.pop_back();
    assert(k == " XX HTML");
    // Two-arg convenience form still works.
    keywrd = " XX     ";
    l_control("ZZ", 1);
    k = keywrd;
    while (!k.empty() && k.back() == ' ') k.pop_back();
    assert(k == " XX ZZ");
    std::fprintf(stderr, "T2 ok\n");
}

// ---------------------------------------------------------------------------
static void write_ref_files() {
    // Full reference file read by geo_ref (unit 99): 3 header lines then atoms.
    {
        std::ofstream f("test_ref.mop");
        f << "GEO_REF test reference\n";
        f << "* comment\n";
        f << " second line of refkey\n";
        f << "C    0.0000    0.0000    0.0000   0   0   0\n";
        f << "C    1.4000    0.0000    0.0000   0   0   0\n";
        f << "H    0.7000    0.5000    0.0000   0   0   0\n";
        f << "H    0.7000   -0.5000    0.0000   0   0   0\n";
        f << "\n";
    }
    // stdin stream for getgeo (which reads the geometry after the header).
    {
        std::ofstream f("test_ref_stdin.txt");
        f << "C    0.0000    0.0000    0.0000   0   0   0\n";
        f << "C    1.4000    0.0000    0.0000   0   0   0\n";
        f << "H    0.7000    0.5000    0.0000   0   0   0\n";
        f << "H    0.7000   -0.5000    0.0000   0   0   0\n";
        f << "\n";
    }
}

static void T3_geo_ref_full() {
    setup();
    write_ref_files();
    // Data-set geometry: deliberately different from the reference.
    numat = 4;
    natoms = 4;
    id = 1;
    nvar = 15;
    keywrd = " GEO_REF=\"test_ref.mop\" 0SCF";
    refkey[1] = "GEO_REF=\"test_ref.mop\" 0SCF";
    refkey[2] = " ";
    refkey[3] = " ";
    refkey[4] = " ";
    refkey[5] = " ";
    refkey[6] = " ";
    for (int i = 1; i <= 4; ++i) {
        nat[i] = (i <= 2) ? 6 : 1;
        labels[i] = nat[i];
    }
    atmass[1] = 12.0; atmass[2] = 12.0; atmass[3] = 1.0; atmass[4] = 1.0;
    geo[1][1] = 0.0;  geo[2][1] = 0.0;  geo[3][1] = 0.0;
    geo[1][2] = 3.0;  geo[2][2] = 0.0;  geo[3][2] = 0.0;
    geo[1][3] = 0.5;  geo[2][3] = 0.5;  geo[3][3] = 0.0;
    geo[1][4] = 0.5;  geo[2][4] = -0.5; geo[3][4] = 0.0;
    coord = geo;
    // loc: parameter i -> atom loc[1][i], coordinate loc[2][i] (1=x,2=y,3=z).
    for (int i = 1; i <= nvar; ++i) {
        loc[1][i] = ((i - 1) % 4) + 1;
        loc[2][i] = ((i - 1) / 4) % 3 + 1;
    }
    std::freopen("test_ref_stdin.txt", "r", stdin);
    std::freopen("test_c6_stdout.txt", "w", stdout);
    geo_ref();
    std::fflush(stdout);
    std::fflush(stdin);

    // 0SCF completion signal.
    assert(moperr == true);
    // gmetry() resets id=0 on a new calculation (no translation vectors here),
    // so numat stays 4 (F90: gmetry sets id; geo_ref does numat -= id).
    assert(numat == 4);
    std::fprintf(stderr, "natoms after geo_ref = %d\n", natoms);
    assert(natoms == 4);
    assert(numat_old == 0);
    // Input geometry restored (coord was saved before docking).
    expect_close("T3 geo1x", geo[1][1], 0.0, 1.e-9);
    expect_close("T3 geo2x", geo[1][2], 3.0, 1.e-9);
    expect_close("T3 geo3y", geo[2][3], 0.5, 1.e-9);
    expect_close("T3 geo4y", geo[2][4], -0.5, 1.e-9);
    // xparam restored from geo via loc (loc[1][i]=atom, loc[2][i]=coord).
    // Test loc mapping: i=1->atom1.x=0, i=2->atom2.x=3, i=3->atom3.x=0.5,
    // i=4->atom4.x=0.5, i=5->atom1.y=0, i=6->atom2.y=0.
    expect_close("T3 xp1", xparam[1], 0.0, 1.e-9);
    expect_close("T3 xp2", xparam[2], 3.0, 1.e-9);
    expect_close("T3 xp5", xparam[5], 0.0, 1.e-9);
    expect_close("T3 xp6", xparam[6], 0.0, 1.e-9);
    // refkey[1] must have been swapped back.
    std::string r1 = refkey[1];
    assert(r1.find("GEO_REF") != std::string::npos);

    // stdout: "After docking" + reference geometry written by geout(99).
    {
        std::ifstream f("test_c6_stdout.txt");
        std::string content((std::istreambuf_iterator<char>(f)),
                            std::istreambuf_iterator<char>());
        assert(content.find("After docking") != std::string::npos);
        assert(content.find("0.0000") != std::string::npos);   // geout writes 4 decimals
        assert(content.find("1.4000") != std::string::npos);
    }
    // The .new file must exist (rotated reference written).
    {
        std::ifstream f("test_dat.new");
        assert(f.good());
    }
    std::fprintf(stderr, "T3 ok\n");
}

// ---------------------------------------------------------------------------
static void T4_geo_ref_missing() {
    setup();
    numat = 4;
    natoms = 4;
    id = 1;
    nvar = 15;
    keywrd = " GEO_REF=\"nonexist_ref.mop\" 0SCF";
    refkey[1] = "GEO_REF=\"nonexist_ref.mop\" 0SCF";
    for (int i = 1; i <= 4; ++i) {
        nat[i] = (i <= 2) ? 6 : 1;
        labels[i] = nat[i];
    }
    for (int i = 1; i <= 4; ++i) {
        geo[1][i] = 0.0; geo[2][i] = 0.0; geo[3][i] = 0.0;
    }
    geo_ref();
    assert(moperr == true);
    assert(numat == 4);  // aborted before numat -= id
    std::fprintf(stderr, "T4 ok\n");
}

int main() {
    T1_geo_diff();
    T2_l_control();
    T3_geo_ref_full();
    T4_geo_ref_missing();
    std::fprintf(stderr, "ALL PASS\n");
    return 0;
}
