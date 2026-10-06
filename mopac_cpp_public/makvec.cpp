// makvec.cpp — C++ translation of MOPAC 2016 "makvec.F90".
// Construct starting occupied/virtual localized MOs on atoms.

#include "makvec.h"

#include <cmath>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

// Forward decls.
int ijbo(int i, int j);
void timer(const std::string& msg);
void memory_error(const char* name);
void buildf(std::vector<double>& f, const std::vector<double>& partf, int zero);
void hybrid(double* catom);
void mbonds(int locc, int lvir, const double* f, double* catom,
            const int* nfirst, int ii, int jj, bool* ui, bool* uj,
            bool& lok, const int* iorbs, double* cocc, double* cvir,
            int cocc_dim, int cvir_dim, int numat, int norbs, int morb,
            int mpack1);
void mlmo(int& locc, int& lvir, int ii, int jj, int& nf_loc, int& ne,
          int& nocc, int& nvir, int* iz, int* ib, int* nce, int* ncf,
          int* ncocc, int* ncvir, const int* iorbs, int* icocc, int* icvir,
          double* cocc, double* cvir);

void makvec() {
    using namespace MOZYME_C;
    using namespace common_arrays_C;
    using namespace molkst_C;

    bool times = (keywrd.find(" TIMES") != std::string::npos);

    p.assign(p.size(), 0.0);
    int m = 0;
    for (int i = 1; i <= numat; ++i) {
        int j = ijbo(i, i);
        for (int k = 1; k <= iorbs[i]; ++k) {
            m = m + 1;
            j = j + k;
            p[j] = pdiag[m];
        }
    }
    if (times) timer(" After entry to Makvec");
    buildf(f, partf, 0);
    if (times) timer(" After BUILDF in MAKV");

    // Temporaries.
    // use global MOZYME_C::iz, ib (Fortran local arrays mirror them)
    std::vector<double> catom((morb * norbs) + 1, 0.0);
    std::vector<char> u(norbs + 1, 0);

    moperr = false;
    hybrid(catom.data());

    int nocc = 0, nvir = 0, nf_loc = 0, ne = 0, locc = 0, lvir = 0;
    for (int i = 1; i <= norbs; ++i) u[i] = 0;  // .false.

    for (int nLewis = 1; nLewis <= Lewis_tot; ++nLewis) {
        int ii = Lewis_elem[0][nLewis];
        int jj = Lewis_elem[1][nLewis];
        if (ii > 0 && jj > 0) {
            int ni = nfirst[ii];
            int nj = nfirst[jj];
            bool lok = false;
            mbonds(locc, lvir, f.data(), catom.data(), nfirst.data(), ii, jj,
                   (bool*)u.data() + ni - 1, (bool*)u.data() + nj - 1, lok,
                   iorbs.data(), cocc.data(), cvir.data(), cocc_dim, cvir_dim,
                   numat, norbs, morb, mpack);
            if (lok) {
                nncf[nocc + 1] = nf_loc;
                nnce[nvir + 1] = ne;
                mlmo(locc, lvir, ii, jj, nf_loc, ne, nocc, nvir,
                     MOZYME_C::iz.data(), MOZYME_C::ib.data(), nce.data(), ncf.data(),
                     ncocc.data(), ncvir.data(), iorbs.data(),
                     icocc.data(), icvir.data(), cocc.data(), cvir.data());
            } else {
                ii = 0;
            }
        } else if (ii > 0) {
            // Lone pair.
            for (int k = nlast[ii]; k >= nfirst[ii]; --k) {
                if (!u[k]) {
                    u[k] = 1;
                    m = locc;
                    for (int m1 = 1; m1 <= iorbs[ii]; ++m1) {
                        m = m + 1;
                        // catom(m1,k): column-major catom(morb,norbs)
                        cocc[m] = catom[(k - 1) * morb + (m1 - 1) + 1];
                    }
                    nncf[nocc + 1] = nf_loc;
                    mlmo(locc, lvir, ii, 0, nf_loc, ne, nocc, nvir,
                         MOZYME_C::iz.data(), MOZYME_C::ib.data(), nce.data(), ncf.data(),
                         ncocc.data(), ncvir.data(), iorbs.data(),
                         icocc.data(), icvir.data(), cocc.data(), cvir.data());
                    break;
                }
            }
        } else {
            // Virtual lone pair.
            for (int k = nlast[jj]; k >= nfirst[jj]; --k) {
                if (!u[k]) {
                    u[k] = 1;
                    m = lvir;
                    for (int m1 = 1; m1 <= iorbs[jj]; ++m1) {
                        m = m + 1;
                        cvir[m] = catom[(k - 1) * morb + (m1 - 1) + 1];
                    }
                    nnce[nvir + 1] = ne;
                    mlmo(locc, lvir, 0, jj, nf_loc, ne, nocc, nvir,
                         MOZYME_C::iz.data(), MOZYME_C::ib.data(), nce.data(), ncf.data(),
                         ncocc.data(), ncvir.data(), iorbs.data(),
                         icocc.data(), icvir.data(), cocc.data(), cvir.data());
                    break;
                }
            }
        }
    }
}
