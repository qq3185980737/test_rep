// test_m10_vecprt.cpp — M10 vecprt group:
// vecprt + vecprt_for_MOZYME.
// Acceptance: ALL PASS + ASan clean + zero warnings.
#define _CRT_SECURE_NO_WARNINGS
#include "vecprt.h"
#include "vecprt_for_MOZYME.h"

#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using common_arrays_C::nat;
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using common_arrays_C::l_atom;
using molkst_C::numat;
using molkst_C::mozyme;
using molkst_C::gui;
using elemts_C::elemnt;

// ijbo stub (host real implementation excluded from this link).
int ijbo(int, int) { return 0; }

static int fails = 0;
#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        if (cond) {                                                         \
            std::fprintf(stderr, "PASS: %s\n", msg);                        \
        } else {                                                            \
            std::fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__);    \
            ++fails;                                                        \
        }                                                                   \
    } while (0)

// t1: generic array (numat==0) — packed 3x3 lower triangle.
void t_vecprt_generic() {
    numat = 0;
    mozyme = false;
    double a[6] = {1.0, 0.5, 0.2, 0.3, 0.4, 0.1};
    std::fflush(stdout);
    FILE* cap = std::freopen("_vecprt_cap1.txt", "w", stdout);
    if (cap) {
        vecprt(a, 3);
        std::fflush(stdout);
        std::fclose(stdout);
    }
    std::string out;
    FILE* f = std::fopen("_vecprt_cap1.txt", "r");
    if (f) {
        char buf[2048];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
    }
    CHECK(out.find("1.000000") != std::string::npos,
          "vecprt: diag 1.000000 printed");
    CHECK(out.find("0.500000") != std::string::npos,
          "vecprt: off-diag 0.500000 printed");
    CHECK(out.find("0.100000") != std::string::npos,
          "vecprt: last 0.100000 printed");
}

// t2: atomic-orbital branch (nlast[numat] == numm).
void t_vecprt_orbitals() {
    numat = 1;
    mozyme = false;
    nat.assign(4, 0);
    nat[1] = 1;
    nfirst.assign(4, 0);
    nlast.assign(4, 0);
    nfirst[1] = 1;
    nlast[1] = 2;
    double a[3] = {1.0, 0.3, 0.2};
    std::fflush(stdout);
    FILE* cap = std::freopen("_vecprt_cap2.txt", "w", stdout);
    if (cap) {
        vecprt(a, 2);
        std::fflush(stdout);
        std::fclose(stdout);
    }
    std::string out;
    FILE* f = std::fopen("_vecprt_cap2.txt", "r");
    if (f) {
        char buf[2048];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
    }
    CHECK(out.find(" H") != std::string::npos, "vecprt: element H shown");
    CHECK(out.find(" S") != std::string::npos, "vecprt: orbital S shown");
    CHECK(out.find("0.300000") != std::string::npos,
          "vecprt: ao off-diag printed");
}

// t3: scaling — large diagonal terms trigger the multiplier message.
void t_vecprt_scale() {
    numat = 0;
    mozyme = false;
    double a[6] = {1.0e5, 0.5, 0.2, 0.3, 0.4, 0.1};
    std::fflush(stdout);
    FILE* cap = std::freopen("_vecprt_cap3.txt", "w", stdout);
    if (cap) {
        vecprt(a, 3);
        std::fflush(stdout);
        std::fclose(stdout);
    }
    std::string out;
    FILE* f = std::fopen("_vecprt_cap3.txt", "r");
    if (f) {
        char buf[2048];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
    }
    CHECK(out.find("Diagonal Terms should be Multiplied by") !=
              std::string::npos,
          "vecprt: scale message present");
    // After scaling + restoring, the input array is unchanged.
    CHECK(std::fabs(a[0] - 1.0e5) < 1e-6, "vecprt: input array restored");
}

// t4: MOZYME variant, Option (1) — l_atom selection over atoms.
void t_vecprt_mozyme() {
    numat = 2;
    mozyme = true;
    gui = false;
    nat.assign(4, 0);
    nat[1] = 1;
    nat[2] = 8;
    l_atom.assign(4, false);
    l_atom[1] = true;
    l_atom[2] = true;
    // 2x2 packed: a11=1, a21=0.3, a22=2 (0-based aa)
    double aa[3] = {1.0, 0.3, 2.0};
    std::fflush(stdout);
    FILE* cap = std::freopen("_vecprt_cap4.txt", "w", stdout);
    if (cap) {
        vecprt_for_MOZYME(aa, 2);
        std::fflush(stdout);
        std::fclose(stdout);
    }
    std::string out;
    FILE* f = std::fopen("_vecprt_cap4.txt", "r");
    if (f) {
        char buf[2048];
        while (std::fgets(buf, sizeof(buf), f)) out += buf;
        std::fclose(f);
    }
    CHECK(out.find(" H") != std::string::npos, "vecprt_for_MOZYME: H shown");
    CHECK(out.find(" O") != std::string::npos, "vecprt_for_MOZYME: O shown");
    CHECK(out.find("2.000000") != std::string::npos,
          "vecprt_for_MOZYME: diag printed");
    mozyme = false;
}

int main() {
    std::fprintf(stderr, "== M10 vecprt group ==\n");
    std::fprintf(stderr, "[t1 generic]\n");
    t_vecprt_generic();
    std::fprintf(stderr, "[t2 orbitals]\n");
    t_vecprt_orbitals();
    std::fprintf(stderr, "[t3 scale]\n");
    t_vecprt_scale();
    std::fprintf(stderr, "[t4 mozyme]\n");
    t_vecprt_mozyme();
    if (fails == 0) {
        std::fprintf(stderr, "ALL PASS: M10 vecprt group\n");
        return 0;
    }
    std::fprintf(stderr, "FAILURES: %d\n", fails);
    return 1;
}
