// test_batchM08g.cpp  M08 meci driver end-to-end + mecip density correction
#include <cmath>
#include <cstdio>
#include <vector>
#include "common_arrays_C.h"
#include "meci_C.h"
#include "meci.h"
#include "mecip.h"
#include "mndod_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

namespace cosmo_C { extern bool iseps; }

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08g batch - meci driver + mecip\n");

    // ---- Part 1: mecip density correction (hand-computed) ----
    {
        using namespace common_arrays_C;
        using namespace meci_C;
        using namespace molkst_C;
        norbs = 2; nmos = 2; lab = 2; nstate = 1; nelec = 0;
        occa.assign({0.0, 2.0, 0.0});
        microa.assign(3, std::vector<int>(3, 0));
        microb.assign(3, std::vector<int>(3, 0));
        microa[1][1] = 2; microa[2][1] = 0;
        microa[1][2] = 1; microa[2][2] = 1;
        microb[1][2] = 1; microb[2][2] = 1;
        nalmat.assign({0, 2, 2});
        vectci.assign(3, 0.0);
        vectci[1] = 0.9; vectci[2] = 0.1;
        c.assign(3, std::vector<double>(3, 0.0));
        c[1][1] = 0.6; c[2][1] = 0.8; c[1][2] = 0.8; c[2][2] = -0.6;
        p.assign(4, 0.0);
        p[1] = 0.36; p[2] = 0.48; p[3] = 0.64;
        deltap.assign(3, std::vector<double>(3, 0.0));
        mecip();
        std::printf("    mecip p = %.6f %.6f %.6f\n", p[1], p[2], p[3]);
        CHK_D(p[1], -0.4768, 1e-9, "mecip p(1,1)");
        CHK_D(p[2], -0.6624, 1e-9, "mecip p(2,1)");
        CHK_D(p[3], -0.8632, 1e-9, "mecip p(2,2)");
    }

    // ---- Part 2: meci driver, 2 MO active space, full CI ----
    {
        using namespace common_arrays_C;
        using namespace meci_C;
        using namespace molkst_C;
        // Reset CI state.
        spin.clear(); nalmat.clear(); eig.clear(); conf.clear();
        vectci.clear(); ispin.clear(); ispqr.clear(); occa.clear();
        deltap.clear(); nfa.clear(); rjkaa.clear(); rjkab.clear();
        eiga.clear(); cimat.clear();

        norbs = 4; numat = 2; lm61 = 8;
        nfirst.assign({0, 1, 3});
        nlast.assign({0, 2, 4});
        nclose = 1; nopen = 3; fract = 2.0; nelecs = 2;
        keywrd = " MECI";
        last = 0; numcal = 7;
        maxci = 70; cdiagi = 0.0;
        eigs.assign(5, 0.0);
        eigs[1] = 0.5; eigs[2] = 1.0; eigs[3] = 1.5; eigs[4] = 2.0;
        c.assign(5, std::vector<double>(5, 0.0));
        for (int i = 1; i <= 4; ++i) c[i][i] = 1.0;
        w.assign(128, 0.0);
        w[1] = 0.4; w[2] = 0.3; w[3] = 0.2;
        mndod_C::fx.assign(30, 0.0);
        mndod_C::fx[0] = 1; mndod_C::fx[1] = 1; mndod_C::fx[2] = 2; mndod_C::fx[3] = 6; mndod_C::fx[4] = 24; mndod_C::fx[5] = 120;
        cosmo_C::iseps = false;
        molkst_C::msdel = 0;
        symmetry_C::namo.assign(8, " ");
        symmetry_C::jndex.assign(8, 0);

        double e0 = meci();
        CHECK(std::isfinite(e0), "meci returns finite energy");
        CHECK(lab == 1, "lab==1 (single microstate)");
        CHECK(occa[1] == 1.0 && occa[2] == 1.0, "occa={1,1}");
        CHECK(nelec == 0, "nelec==0");
        CHECK(root_requested == 1, "root_requested==1");
        CHECK(nstate == 1, "nstate==1");
        CHECK(ispin[1] == 1, "ispin(1)==1 (singlet)");
        CHECK(symmetry_C::jndex[1] == 1, "jndex(1)==1");
        // eig[1] must equal the returned CI energy.
        CHECK(std::fabs(e0 - eig[1]) < 1e-12, "meci()==eig(1)");
        // rsp normalizes the single vector: conf(0)=1.
        CHECK(std::fabs(conf[0] - 1.0) < 1e-12, "conf(1,1)=1");

        // Second call: icalcn==numcal path (first1=false).
        double e1 = meci();
        CHECK(std::isfinite(e1), "meci second call finite");
        CHECK(e1 == e0, "meci second call reproducible");
    }

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
