// big_swap.cpp — C++ translation of MOPAC 2016 "big_swap.F90" (1320 lines).
//
// Protein transition-state location module, faithful port:
//   Locate_TS_for_Proteins, lbfgs_TS, compfg_TS, get_pars, big_swap
//   (+ internal copy_i_1/copy_i_2/copy_r_1/copy_r_2), build_active_site,
//   select_opt, Refine_TS_for_Proteins, l_control.
//
// Adaptation notes:
//  - 1-based Fortran indexing kept throughout (vectors n+1, matrices rows n+1).
//  - F90 write(iw,...) -> stdout; endfile/backspace -> no-ops.
//  - F90 open/write/delete of scratch files uses C stdio (remove() for delete).
//  - get_pars reads calibration data from the file named by "calib.dat" when
//    present, else from stdin-style unit-33 semantics are approximated.

#include "big_swap.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "blas1.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "compfg.h"
#include "cosmo_C.h"
#include "ef_C.h"
#include "ef.h"
#include "elemts_C.h"
#include "geo_diff.h"
#include "geout.h"
#include "lbfgs.h"
#include "memory_error.h"
#include "moldat.h"
#include "molkst_C.h"
#include "mopend.h"
#include "MOZYME_C.h"
#include "prttim.h"
#include "run_mopac_deps_stubs.h"
#include "reada.h"
#include "second.h"
#include "set_up_dentate.h"
#include "set_up_RAPID.h"
#include "timer.h"
#include "to_screen.h"
#include "upcase.h"
#include "write_cell.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace cosmo_C;
using namespace ef_C;
using namespace elemts_C;
using namespace molkst_C;
using namespace MOZYME_C;

// Fortran NINT: nearest integer (round half away from zero).
static inline int nint(double x) {
  return (x < 0.0) ? (int)(x - 0.5) : (int)(x + 0.5);
}

namespace bs_impl {

// 1-based pointer BLAS helpers (element i at p[1+(i-1)*inc]).
inline double ddot1b(int n, const double* x, int incx, const double* y, int incy) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[(size_t)i * incx] * y[(size_t)i * incy];
  return s;
}
inline void dcopy1b(int n, const double* x, int incx, double* y, int incy) {
  for (int i = 0; i < n; ++i) y[(size_t)i * incy] = x[(size_t)i * incx];
}
inline void daxpy1b(int n, double a, const double* x, int incx, double* y, int incy) {
  for (int i = 0; i < n; ++i) y[(size_t)i * incy] += a * x[(size_t)i * incx];
}
// F90 dot(a, b, n) = sum_i a(i)*b(i), 1-based pointers.
inline double dot1b(int n, const double* a, const double* b) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += a[i] * b[i];
  return s;
}

}  // namespace bs_impl

using bs_impl::ddot1b;
using bs_impl::dcopy1b;
using bs_impl::daxpy1b;
using bs_impl::dot1b;

// ============================================================================
// l_control(txt, nt, mode): add (mode=1) or remove (mode=-1) the word(s) in
// txt from keywrd.  If a word already exists when adding, the old occurrence
// is deleted first.  nt is the declared length of txt.
// ============================================================================
void l_control(const std::string& txt_in, int nt, int mode) {
  std::string local_txt = txt_in;
  while (true) {
    size_t s0 = local_txt.find_first_not_of(' ');
    if (s0 == std::string::npos) break;
    local_txt = local_txt.substr(s0);
    int limit = (int)std::min((size_t)nt, local_txt.size());
    if (limit < 1) break;
    char ch = ' ';
    int mt = 1;
    for (; mt <= limit; ++mt) {
      char c = local_txt[mt - 1];
      if (c == ch) break;
      if (c == '"') ch = '"';
    }
    if (ch == ' ') mt = mt - 1;
    if (mt < 1) break;
    std::string word = local_txt.substr(0, mt);
    if (mt + 1 >= (int)local_txt.size()) local_txt.clear();
    else local_txt = local_txt.substr(mt + 1);
    int i;
    if (word[0] >= '0' && word[0] <= '9') {
      i = mt + 1;
    } else {
      i = 1;
      for (; i <= mt; ++i) {
        char c = word[i - 1];
        if (!((c >= 'A' && c <= 'Z') || c == '_' || c == '-')) break;
      }
    }
    int trim_len = i - 1;
    std::string stem = word.substr(0, trim_len);
    while (true) {
      std::string needle = " " + stem;
      size_t idx = keywrd.find(needle);
      if (idx == std::string::npos) break;
      size_t j = idx + 1;
      while (j < keywrd.size() && keywrd[j] != ' ') ++j;
      std::string tail = (j + 1 < keywrd.size()) ? keywrd.substr(j + 1) : std::string();
      keywrd = keywrd.substr(0, idx) + tail;
    }
    if (mode == 1) {
      // Fortran: keywrd is a fixed-length string padded with blanks, so the
      // insert point found by index(keywrd, blank(:mt+2)) is the padded tail
      // and existing keywords are preserved.  C++ keywrd is a variable-length
      // string, so append the word at the end instead of slicing keywrd.
      if (keywrd.empty()) {
        keywrd = word;
      } else if (keywrd.back() == ' ') {
        keywrd += word;
      } else {
        keywrd += " " + word;
      }
    }
    if (local_txt.find_first_not_of(' ') == std::string::npos) break;
  }
}

// Two-argument convenience form (nt inferred from the string length).
void l_control(const std::string& txt, int mode) {
  l_control(txt, (int)txt.size(), mode);
}

// ============================================================================
// big_swap — store (mode=0) or extract (mode=1) all arrays of system 1 or 2.
// ============================================================================
namespace bs_impl {
static void copy_i_1(const std::vector<int>& from, std::vector<int>& to) {
  if (!from.empty()) to = from;
}
static void copy_i_2(const std::vector<std::vector<int>>& from,
                     std::vector<std::vector<int>>& to, int dim_1) {
  if (!from.empty()) {
    int n2 = (int)from.size() / dim_1;
    to.assign((size_t)dim_1 + 1, std::vector<int>((size_t)n2, 0));
    for (int a = 1; a <= dim_1; ++a)
      for (int b = 1; b < n2; ++b) to[a][b] = from[a][b];
  }
}
static void copy_r_1(const std::vector<double>& from, std::vector<double>& to) {
  if (!from.empty()) to = from;
}
static void copy_r_2(const std::vector<std::vector<double>>& from,
                     std::vector<std::vector<double>>& to, int dim_1) {
  if (!from.empty() && !from[0].empty()) {
    int n2 = (int)from[0].size();
    to.assign((size_t)dim_1 + 1, std::vector<double>((size_t)n2, 0.0));
    for (int a = 1; a <= dim_1; ++a)
      for (int b = 1; b < n2; ++b) to[a][b] = from[a][b];
  }
}
}  // namespace bs_impl

void big_swap(int mode, int system) {
  static int icocc_dim_1 = 0, icvir_dim_1 = 0, cvir_dim_1 = 0, cocc_dim_1 = 0;
  static int mpack_1 = 0, norred_1 = 0, nelred_1 = 0;
  static int icocc_dim_2 = 0, icvir_dim_2 = 0, cvir_dim_2 = 0, cocc_dim_2 = 0;
  static int mpack_2 = 0, norred_2 = 0, nelred_2 = 0;
  static double refnuc_1 = 0.0, refnuc_2 = 0.0, solv_energy_1 = 0.0, solv_energy_2 = 0.0;
  using bs_impl::copy_i_1;
  using bs_impl::copy_i_2;
  using bs_impl::copy_r_1;
  using bs_impl::copy_r_2;
  if (mode == 0) {
    // Store system
    if (system == 1) {
      if (!nbonds.empty()) {
        nbonds_1 = nbonds;
        ibonds_1 = ibonds;
      }
      icocc_dim_1 = icocc_dim;
      icvir_dim_1 = icvir_dim;
      cvir_dim_1 = cvir_dim;
      cocc_dim_1 = cocc_dim;
      norred_1 = norred;
      nelred_1 = nelred;
      refnuc_1 = refnuc;
      solv_energy_1 = solv_energy;
      mpack_1 = mpack;
      copy_r_2(geo, geo_1, 3);
      copy_r_1(dxyz, dxyz_1);
      copy_i_1(icocc, icocc_1);
      copy_i_1(icvir, icvir_1);
      copy_i_1(ncocc, ncocc_1);
      copy_i_1(ncvir, ncvir_1);
      copy_i_1(nncf, nncf_1);
      copy_i_1(nnce, nnce_1);
      copy_i_1(ncf, ncf_1);
      copy_i_1(nce, nce_1);
      copy_i_1(iij, iij_1);
      copy_i_1(iijj, iijj_1);
      copy_i_1(ijall, ijall_1);
      copy_i_1(numij, numij_1);
      copy_i_1(iorbs, iorbs_1);
      copy_i_2(nijbo, nijbo_1, numat);
      copy_r_1(cocc, cocc_1);
      copy_r_1(cvir, cvir_1);
      copy_r_1(xparam, xparam_1);
      copy_r_1(partf, partf_1);
      copy_r_1(partp, partp_1);
      copy_r_1(parth, parth_1);
      copy_r_1(f, f_1);
      copy_r_1(p, p_1);
    } else {
      if (!nbonds.empty()) {
        nbonds_2 = nbonds;
        ibonds_2 = ibonds;
      }
      icocc_dim_2 = icocc_dim;
      icvir_dim_2 = icvir_dim;
      cvir_dim_2 = cvir_dim;
      cocc_dim_2 = cocc_dim;
      norred_2 = norred;
      nelred_2 = nelred;
      refnuc_2 = refnuc;
      solv_energy_2 = solv_energy;
      mpack_2 = mpack;
      copy_r_2(geo, geo_2, 3);
      copy_r_1(dxyz, dxyz_2);
      copy_i_1(icocc, icocc_2);
      copy_i_1(icvir, icvir_2);
      copy_i_1(ncocc, ncocc_2);
      copy_i_1(ncvir, ncvir_2);
      copy_i_1(nncf, nncf_2);
      copy_i_1(nnce, nnce_2);
      copy_i_1(ncf, ncf_2);
      copy_i_1(nce, nce_2);
      copy_i_1(iij, iij_2);
      copy_i_1(iijj, iijj_2);
      copy_i_1(ijall, ijall_2);
      copy_i_1(numij, numij_2);
      copy_i_1(iorbs, iorbs_2);
      copy_i_2(nijbo, nijbo_2, numat);
      copy_r_1(cocc, cocc_2);
      copy_r_1(cvir, cvir_2);
      copy_r_1(xparam, xparam_2);
      copy_r_1(partf, partf_2);
      copy_r_1(partp, partp_2);
      copy_r_1(parth, parth_2);
      copy_r_1(f, f_2);
      copy_r_1(p, p_2);
    }
  } else {
    // Extract system
    if (system == 1) {
      for (int i = 1; i <= numat; ++i) nbonds[i] = nbonds_1[i];
      for (int j = 1; j <= 15; ++j)
        for (int i = 1; i <= numat; ++i) ibonds[j][i] = ibonds_1[j][i];
      if (icocc_dim_1 > 0) icocc_dim = icocc_dim_1;
      if (icvir_dim_1 > 0) icvir_dim = icvir_dim_1;
      if (cvir_dim_1 > 0) cvir_dim = cvir_dim_1;
      if (cocc_dim_1 > 0) cocc_dim = cocc_dim_1;
      if (norred_1 > 0) norred = norred_1;
      if (nelred_1 > 0) nelred = nelred_1;
      if (refnuc_1 > 0) refnuc = refnuc_1;
      if (std::fabs(solv_energy_1) > 0.1) solv_energy = solv_energy_1;
      if (mpack_1 > 0) mpack = mpack_1;
      copy_r_2(geo_1, geo, 3);
      copy_r_1(dxyz_1, dxyz);
      copy_i_1(icocc_1, icocc);
      copy_i_1(icvir_1, icvir);
      copy_i_1(ncocc_1, ncocc);
      copy_i_1(ncvir_1, ncvir);
      copy_i_1(nncf_1, nncf);
      copy_i_1(nnce_1, nnce);
      copy_i_1(ncf_1, ncf);
      copy_i_1(nce_1, nce);
      copy_i_1(iij_1, iij);
      copy_i_1(iijj_1, iijj);
      copy_i_1(ijall_1, ijall);
      copy_i_1(numij_1, numij);
      copy_i_1(iorbs_1, iorbs);
      copy_i_2(nijbo_1, nijbo, numat);
      copy_r_1(cocc_1, cocc);
      copy_r_1(cvir_1, cvir);
      copy_r_1(xparam_1, xparam);
      copy_r_1(partf_1, partf);
      copy_r_1(partp_1, partp);
      copy_r_1(parth_1, parth);
      copy_r_1(f_1, f);
      copy_r_1(p_1, p);
      pa = p;
      for (size_t k = 0; k < pa.size(); ++k) pa[k] *= 0.5;
      pb = pa;
    } else {
      for (int i = 1; i <= numat; ++i) nbonds[i] = nbonds_2[i];
      for (int j = 1; j <= 15; ++j)
        for (int i = 1; i <= numat; ++i) ibonds[j][i] = ibonds_2[j][i];
      if (icocc_dim_2 > 0) icocc_dim = icocc_dim_2;
      if (icvir_dim_2 > 0) icvir_dim = icvir_dim_2;
      if (cvir_dim_2 > 0) cvir_dim = cvir_dim_2;
      if (cocc_dim_2 > 0) cocc_dim = cocc_dim_2;
      if (norred_2 > 0) norred = norred_2;
      if (nelred_2 > 0) nelred = nelred_2;
      if (refnuc_2 > 0) refnuc = refnuc_2;
      if (std::fabs(solv_energy_2) > 0.1) solv_energy = solv_energy_2;
      if (mpack_2 > 0) mpack = mpack_2;
      copy_r_2(geo_2, geo, 3);
      copy_r_1(dxyz_2, dxyz);
      copy_i_1(icocc_2, icocc);
      copy_i_1(icvir_2, icvir);
      copy_i_1(ncocc_2, ncocc);
      copy_i_1(ncvir_2, ncvir);
      copy_i_1(nncf_2, nncf);
      copy_i_1(nnce_2, nnce);
      copy_i_1(ncf_2, ncf);
      copy_i_1(nce_2, nce);
      copy_i_1(iij_2, iij);
      copy_i_1(iijj_2, iijj);
      copy_i_1(ijall_2, ijall);
      copy_i_1(numij_2, numij);
      copy_i_1(iorbs_2, iorbs);
      copy_i_2(nijbo_2, nijbo, numat);
      copy_r_1(cocc_2, cocc);
      copy_r_1(cvir_2, cvir);
      copy_r_1(xparam_2, xparam);
      copy_r_1(partf_2, partf);
      copy_r_1(partp_2, partp);
      copy_r_1(parth_2, parth);
      copy_r_1(f_2, f);
      copy_r_1(p_2, p);
      pa = p;
      for (size_t k = 0; k < pa.size(); ++k) pa[k] *= 0.5;
      pb = pa;
    }
    if ((int)coord.size() < 4)
      coord.assign(4, std::vector<double>((size_t)numat + 1, 0.0));
    for (int j = 1; j <= 3; ++j)
      for (int i = 1; i <= numat; ++i) coord[j-1][i] = geo[j][i];
  }
}


// ============================================================================
// Locate_TS_for_Proteins — multi-step transition-state search using two
// stationary points (data-set geometry vs. GEO_REF geometry).
// ============================================================================
void Locate_TS_for_Proteins() {
  // 1-based Fortran-style indexing for keywrd parsing.
  auto chr = [&](int p1) -> char { return keywrd[p1 - 1]; };

  std::vector<int> active_site(201, 0);
  int ninsite[4] = {0, 0, 0, 0};
  if ((int)geoa.size() < 3) geoa.assign(3, std::vector<double>((size_t)numat + 1, 0.0));
  for (int j = 1; j <= 3; ++j)
    for (int i = 1; i <= numat; ++i) geoa[j][i] = geo_1[j][i];
  {
    std::vector<int> nsi(4, 0);
    build_active_site(active_site, nsi);
    ninsite[1] = nsi[1];
    ninsite[2] = nsi[2];
    ninsite[3] = nsi[3];
  }
  if (moperr) return;

  int nset = 1;
  int shell = 0;
  int loop = 0, nloop = 0;
  int store_mpack = 0, store_n2elec = 0;
  std::string line1;
  std::string region;   // hoisted: goto label_99 jumps past later inits
  bool extra_print = false;
  bool exists = false;
  bool l_refine = true;
  char num = '4';
  double stresses[41], gradients[41];
  for (int i = 0; i < 41; ++i) { stresses[i] = 0.0; gradients[i] = 0.0; }
  int big_nvar = 0;
  std::vector<double> big_xparam;
  // Hoisted declarations (C2362: goto label_99 must not skip initialization).
  int p0 = 0, i = 0, j = 0, k = 0;

  if (keywrd.find(" LOCATE-TS(SET") != std::string::npos) {
    shell = ninsite[1];
    grad.assign((size_t)3 * numat + 1, 0.0);
    goto label_99;
  }
  for (int iq = 1; iq <= 40; ++iq)
    gradients[iq] = std::max(4.0, std::min(20.0, std::sqrt(numat * 0.5)));

  // Parse parameters after " LOCATE-TS".
  p0 = (int)keywrd.find(" LOCATE-TS");               // 0-based keyword position
  i = p0 + 11;                                        // 1-based: just past keyword
  j = i + 1;
  while (j <= std::min(240, i + 100)) {
    if (chr(j) == ' ') break;
    ++j;
  }
  // param region is keywrd[i+1 .. j-1] in 1-based; build 0-based substr
  region = keywrd.substr(i, j - i);
  k = (int)region.find("C:");
  if (k >= 0) {
    extra_print = true;
    i = i + k + 1;
    nloop = 0;
    while (true) {
      k = (int)chr(i) - (int)'0';
      if (k >= 0 && k < 10) {
        ++nloop;
        stresses[nloop] = reada(keywrd, i);
        ++i;
        while (true) {
          char c = chr(i);
          if (c != '.' && (c < '0' || c > '9')) break;
          ++i;
        }
        ++i;
      } else {
        break;
      }
      if (i > j) break;
    }
  } else {
    extra_print = false;
    stresses[1] = 3.0; stresses[2] = 30.0; stresses[3] = 30.0;
    stresses[4] = 30.0; stresses[5] = 30.0;
    nloop = 5;
  }

  // Second scan: "SET" selection.
  i = p0 + 11;
  j = i + 1;
  while (j <= (int)keywrd.size()) {
    if (chr(j) == ' ') break;
    ++j;
  }
  region = keywrd.substr(i, j - i);
  k = (int)region.find("SET");
  if (k >= 0) {
    nset = (int)nint(reada(keywrd, i + k + 3));
    std::printf("\n          Set%2d selected by keywork 'SET' within keyword 'LOCATE-TS'\n", nset);
    if (nset > 2) {
      std::printf("\n          (This is greater than 2, so re-set to 2)\n");
      nset = 2;
    }
  } else {
    if ((int)region.find("C:") < 0) {
      std::printf("\n          By default, set 1 will be used.  To change default, see keyword 'LOCATE-TS'\n");
    } else {
      l_refine = false;
    }
  }

  big_nvar = 2 * nvar;
  big_xparam.assign((size_t)big_nvar + 1, 0.0);
  if (nset > 0) {
    shell = ninsite[nset];
    select_opt(shell, active_site);
    big_nvar = 2 * nvar;
    for (i = 1; i <= nvar; ++i) {
      k = loc[1][i];
      int l = loc[2][i];
      big_xparam[i] = geo_1[l][k];
      big_xparam[i + nvar] = geo_2[l][k];
    }
  } else {
    shell = 0;
  }

  if (nloop == 0) {
    // Use average of the two geometries.
    for (i = 1; i <= nvar; ++i) {
      xparam[i] = 0.5 * (big_xparam[i] + big_xparam[i + nvar]);
      k = loc[1][i];
      int l = loc[2][i];
      geo[l][k] = xparam[i];
    }
    if (extra_print) {
      std::printf("\n          Average of the input and reference geometries\n");
      geout(iw);
    }
  } else {
    // Locate the transition state by climbing the barrier.
    density = 0.0;
    big_swap(1, 2);
    ++numcal;
    ++step_num;
    set_up_rapid("ON");
    set_up_rapid("OFF");
    big_swap(0, 2);
    for (int jj = 1; jj <= 3; ++jj)
      for (int ii = 1; ii <= numat; ++ii) geoa[jj][ii] = geo[jj][ii];

    big_swap(1, 1);
    ++numcal;
    ++step_num;
    set_up_rapid("ON");
    set_up_rapid("OFF");
    for (loop = 1; loop <= nloop; ++loop) {
      density = stresses[loop];
      if (density > 99.9499) num = '6';
      else if (density > 9.9499) num = '5';
      else num = '4';
      std::printf("\n          Constraining constant: %.2f Kcal/mol/Angstrom^2\n", density);
      std::string ltmp;
      {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "C: %.2f", density);
        ltmp = buf;
      }
      char l2[64];
      std::snprintf(l2, sizeof(l2), "GNORM=%.1f", gradients[loop]);
      ltmp += "   ";
      ltmp += l2;
      std::string l3 = l2;
      l_control(l3, (int)l3.size(), 1);
      extra_print = (extra_print || loop == nloop);
      lbfgs_TS(big_xparam, big_nvar, escf, extra_print);
    }
  }

  if (shell == 0) {
    if (nset == 0) mopend("Gradient minimization not requested");
    else mopend("No bonds made or broken in active site");
    return;
  }
label_99:
  if (!l_refine) goto label_98;
  use_ref_geo = false;

  density = 0.0;
  for (int jj = 1; jj <= 3; ++jj)
    for (int ii = 1; ii <= numat; ++ii) geoa[jj][ii] = geo[jj][ii];
  {
    std::string lt = "LET DDMIN=0.D0 GEO-OK";
    l_control(lt, (int)lt.size(), 1);
  }
  ++numcal;
  store_mpack = mpack;
  store_n2elec = n2elec;
  moldat(1);
  mpack = store_mpack;
  n2elec = store_n2elec;
  for (int jj = 1; jj <= 3; ++jj)
    for (int ii = 1; ii <= numat; ++ii) lopt[jj][ii] = 1;
  for (i = 1; i <= shell; ++i) {
    int ai = active_site[i];
    for (int jj = 1; jj <= 3; ++jj) lopt[jj][ai] = 0;
  }
  Refine_TS_for_Proteins();
  // close(iarc) -> no-op

label_98:
  line = input_fn.substr(0, input_fn.size() - 4) + "den";
  if (!line.empty()) {
    if (std::remove(line.c_str()) == 0) { /* deleted */ }
  }
}


// ============================================================================
// lbfgs_TS — dual-geometry L-BFGS optimization (climb the barrier between the
// two stationary points under a GEO_REF stress constraint).
// ============================================================================
void lbfgs_TS(std::vector<double>& big_xparam, int big_nvar, double& escf_tot,
              bool extra_print) {
  int m = 12;
  nstep = 0;
  int niwa = 3 * big_nvar;
  int nwa = 2 * big_nvar * m + 4 * big_nvar + 11 * m * m + 8 * m;
  std::vector<double> bot((size_t)big_nvar + 1, 0.0), gold((size_t)big_nvar + 1, 0.0);
  std::vector<double> top((size_t)big_nvar + 1, 0.0), xold((size_t)big_nvar + 1, 0.0);
  std::vector<double> wa((size_t)nwa + 1, 0.0);
  std::vector<int> iwa((size_t)niwa + 1, 0);
  std::vector<int> nbd((size_t)big_nvar + 1, 0);
  std::vector<double> big_grad((size_t)big_nvar + 1, 0.0);
  std::vector<double> store_big_grad((size_t)big_nvar + 1, 0.0);

  std::string task = " Unused";
  std::string csave = " Unused";
  std::vector<int> lsave(5, 0);
  std::vector<int> isave(45, 0);
  std::vector<double> dsave(30, 0.0);

  bool times = (keywrd.find(" TIMES") != std::string::npos);
  double tolerg;
  if (keywrd.find("GNORM=") != std::string::npos) {
    tolerg = reada(keywrd, (int)keywrd.find("GNORM=") + 7);
    if (keywrd.find(" LET") == std::string::npos && tolerg < 1.e-2) {
      std::printf("\n  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01\n");
      tolerg = 1.e-2;
    }
  } else {
    tolerg = 1.0;
    if (id != 0) tolerg = id * 2.0 - 1.0;
    if (keywrd.find(" PREC") != std::string::npos) tolerg = tolerg * 0.2;
  }
  double oldstp[13];
  for (int i = 1; i <= 10; ++i) oldstp[i] = 1.0;
  double tlast = tleft;
  double tx2 = second(2);
  double tx1 = tx2;
  for (int i = 1; i <= big_nvar; ++i) nbd[i] = 0;
  task = "START";
  int jcyc = 0;
  double cycmx = 0.0;
  tlast = tleft;
  int itry1 = 0;
  double absmin = 1.e6;
  double tstep = 0.0;
  double tt0 = 0.0;
  double sum = 0.0, rms = 0.0, sum1 = 0.0, sum2 = 0.0;
  double e_stress = 0.0, escf1 = 0.0, escf2 = 0.0;
  char txt = ' ';
  double tprt = 0.0;
  bool first = true;
  int nflush = 1;

  while (true) {
    if (times) timer(" Before SETULB");
    dcopy1b(big_nvar, &big_xparam[1], 1, &xold[1], 1);
    setulb(big_nvar, m, &big_xparam[1], &bot[1], &top[1], &nbd[1], escf_tot,
           &big_grad[1], 0.0, 0.0, wa, iwa, task, -1, csave, lsave, isave, dsave);
    if (moperr) goto label_99;
    if (times) timer(" AFTER SETULB");
    dsave[2] = dsave[2] + 1.e4;
    if (task.substr(0, 2) == "FG") {
      if (jcyc > 1) {
        sum = 0.0;
        for (int i = 1; i <= big_nvar; ++i) sum += (big_xparam[i] - xold[i]) * (big_xparam[i] - xold[i]);
        sum = std::sqrt(sum);
        int i = std::min(big_nvar, 10);
        double stepmx = std::min(1.0, std::sqrt(ddot1b(i, &oldstp[1], 1, &oldstp[1], 1) / i));
        if (sum > stepmx * 2.0) {
          double const_f = 2.0 * stepmx / sum;
          for (i = 1; i <= big_nvar; ++i)
            big_xparam[i] = const_f * big_xparam[i] + (1.0 - const_f) * xold[i];
          sum = 2.0 * stepmx;
        }
        oldstp[(jcyc % 10) + 1] = sum;
      }
      // Limit step to 0.2 Angstroms.
      sum = 0.2;
      for (int i = 1; i <= big_nvar; ++i) {
        double dlt = big_xparam[i] - xold[i];
        if (std::fabs(dlt) > sum)
          big_xparam[i] = xold[i] + std::max(-sum, std::min(sum, dlt));
      }
      if (first) {
        geo_diff(sum, rms, false);
        std::printf("\n          Current value of GEO_REF constraint:%20.2f  Kcal/mol/Angstrom^2\n", density);
        std::printf("          Distance between the geometries:%24.2f  Angstroms\n", sum);
      }
      compfg_TS(big_xparam, (jcyc % 222) == 0, escf1, escf2, true, big_grad, true);
      if (first) {
        first = false;
        for (int i = 1; i <= big_nvar; ++i) store_big_grad[i] = big_grad[i];
        e_stress = 0.0;
        int k = 0;
        for (int i = 1; i <= numat; ++i) {
          for (int l = 1; l <= 3; ++l) {
            ++k;
            big_grad[k] = big_grad[k] + (geo[l][i] - geoa[l][i]) * density * 2.0;
            big_grad[k + nvar] = big_grad[k + nvar] + (geoa[l][i] - geo[l][i]) * density * 2.0;
            e_stress = e_stress + (geo[l][i] - geoa[l][i]) * (geo[l][i] - geoa[l][i]) * density;
          }
        }
        std::printf("\n          Heat of formation of the first geometry:%17.3f Kcal/mol\n", escf1 - e_stress);
        std::printf("          Heat of formation of the second geometry:%16.3f Kcal/mol\n", escf2 - e_stress);
        std::printf("          Contribution to heat of formation due to stress:%8.2f  Kcal/mol\n", 2.0 * e_stress);
        sum1 = std::sqrt(dot1b(nvar, &big_grad[1], &big_grad[1]));
        sum2 = std::sqrt(dot1b(nvar, &big_grad[nvar + 1], &big_grad[nvar + 1]));
        std::printf("          Gradient arising from first geometry:%19.2f  Kcal/mol/Angstrom\n", sum1);
        std::printf("          Gradient arising from second geometry:%18.2f  Kcal/mol/Angstrom\n", sum2);
        sum = dot1b(nvar, &big_grad[1], &big_grad[nvar + 1]) / (sum1 * sum2);
        if (sum < 0.0)
          std::printf("          Angle between gradient vectors:%25.2f  degrees\n", std::acos(sum) * 57.2957795);
        else
          std::printf("          WARNING! - Angle between gradient vectors:%13.2f  degrees\n", std::acos(sum) * 57.2957795);
        for (int i = 1; i <= big_nvar; ++i) big_grad[i] = store_big_grad[i];
      }
      escf_tot = escf1 + escf2;
      if (moperr) goto label_99;
      if (absmin - escf_tot < 1.e-7) {
        ++itry1;
        if (itry1 > 900 || (gnorm < 1.0 && itry1 > 9)) {
          std::printf("\n\n HEAT OF FORMATION IS ESSENTIALLY STATIONARY\n");
          iflepo = 3;
          break;
        }
      } else {
        itry1 = 0;
        absmin = escf_tot;
      }
      if (times) timer(" AFTER COMPFG_TS");
      ++jcyc;
      tx2 = second(2);
      tstep = tx2 - tx1;
      cycmx = std::max(tstep, cycmx);
      tx1 = tx2;
      tleft = tleft - tstep;
      if (tlast - tleft > tdump) {
        tlast = tleft;
        tt0 = second(1) - time0;
        lbfsav(tt0, 1, wa, nwa, iwa, niwa, task, csave, lsave, isave, dsave, jcyc, escf_tot);
        if (moperr) goto label_99;
      }
      tleft = std::max(0.0, tleft);
      prttim(tleft, tprt, txt);
      gnorm = std::sqrt(ddot1b(big_nvar, &big_grad[1], 1, &big_grad[1], 1));
      if (id == 3) {
        write_cell(iw);
        write_cell(iw0);
      }
      ++nstep;
      char linebuf[160];
      std::snprintf(linebuf, sizeof(linebuf),
                    " CYCLE:%6d TIME:%8.3f TIME LEFT:%6.2f%c  GRAD.:%10.3f HEAT:%14.7e",
                    jcyc, std::min(tstep, 9999.99), tprt, txt, std::min(gnorm, 999999.999), escf_tot);
      line = linebuf;
      std::printf("%s\n", line.c_str());
      // endfile/backspace -> no-op
      if (chanel_C::log) std::printf("%s\n", line.c_str());
      to_screen(line);
      if (nflush != 0) {
        if (jcyc % nflush == 0) {
          // endfile/backspace -> no-op
        }
      }
      to_screen("To_file: Geometry optimizing");
      dcopy1b(big_nvar, &big_grad[1], 1, &gold[1], 1);
      if (gnorm < tolerg) {
        iflepo = 3;
        break;
      }
    } else if (task.substr(0, 5) != "NEW_X") {
      std::printf(" L-BFGS Message: %s\n", task.c_str());
      iflepo = 9;
      break;
    }
  }
  if (gnorm < tolerg)
    std::printf("\n      GRADIENT =%9.5f  IS LESS THAN CUTOFF =%9.5f\n\n", gnorm, tolerg);

label_99:
  geo_diff(sum, rms, false);
  e_stress = 0.0;
  {
    int k = 0;
    for (int i = 1; i <= numat; ++i) {
      for (int l = 1; l <= 3; ++l) {
        ++k;
        big_grad[k] = big_grad[k] + (geo[l][i] - geoa[l][i]) * density * 2.0;
        big_grad[k + nvar] = big_grad[k + nvar] + (geoa[l][i] - geo[l][i]) * density * 2.0;
        e_stress = e_stress + (geo[l][i] - geoa[l][i]) * (geo[l][i] - geoa[l][i]) * density;
      }
    }
  }
  if (density > 99.9499) txt = '5';
  else if (density > 9.9499) txt = '4';
  else txt = '3';

  line = refkey[1];
  upcase(line, (int)line.size());
  k = (int)line.find(" SETUP");
  if (k >= 0) {
    int l = k + 6;
    while (l < (int)line.size() && line[l] != ' ') ++l;
    refkey[1] = refkey[1].substr(0, k) + refkey[1].substr(l);
    line = refkey[1];
    upcase(line, (int)line.size());
  }
  k = (int)line.find("GEO_DAT=");
  if (k >= 0) {
    int kk = k + 8;                       // 1-based "GEO_DAT=" end
    int l = (int)line.find("\" ", kk);    // closing quote-blank
    if (l == std::string::npos) l = (int)line.size();
    refkey[1] = refkey[1].substr(0, kk - 1) + refkey[1].substr(l + 1);
    line = line.substr(0, kk - 1) + line.substr(l + 1);
  }

  int k_save = 0, l_save = 0;
  if (extra_print) {
    k = (int)line.find("GEO_REF=");
    k_save = k;
    if (k >= 0) {
      int kk = k + 8;
      int l = (int)line.find("\"", kk);
      l_save = l;
      line = input_fn.substr(0, input_fn.size() - 5);
      char buf[128];
      std::snprintf(buf, sizeof(buf), " %.1f first.mop", density);
      line += buf;
      int ii = (int)line.size() - 11;
      line[ii] = 'p';
      refkey[1] = refkey[1].substr(0, kk - 1) + "GEO_REF=\"" + line + "\"" +
                  ((l >= 0 && (size_t)l < refkey[1].size()) ? refkey[1].substr(l) : "");
      ii = (int)refkey[1].find("first.mop");
      if (ii >= 0)
        refkey[1] = refkey[1].substr(0, ii) + "second" + refkey[1].substr(ii + 5);
      add_path(line);
      std::printf("\n          First geometry (derived from data-set) after optimization subject to \n");
      std::printf("          GEO_REF constraint of %.1f Kcal/mol/Angstrom^2 towards the reference geometry written to file:\n", density);
      std::printf("          '%s'\n", line.c_str());
      geout(iarc);
      for (int jj = 1; jj <= 3; ++jj)
        for (int ii2 = 1; ii2 <= numat; ++ii2) geo[jj][ii2] = geoa[jj][ii2];
      line = refkey[1];
      upcase(line, (int)line.size());
      k = (int)line.find("GEO_REF=");
      if (k >= 0) {
        int kk = k + 8;
        int l = (int)line.find("\"", kk);
        line = input_fn.substr(0, input_fn.size() - 5);
        std::snprintf(buf, sizeof(buf), " %.1f second.mop", density);
        line += buf;
        ii = (int)line.size() - 12;
        line[ii] = 'p';
        refkey[1] = refkey[1].substr(0, kk - 1) + "GEO_REF=\"" + line + "\"" +
                    ((l >= 0 && (size_t)l < refkey[1].size()) ? refkey[1].substr(l) : "");
        ii = (int)refkey[1].find("second.mop");
        if (ii >= 0)
          refkey[1] = refkey[1].substr(0, ii) + "first" + refkey[1].substr(ii + 6);
        add_path(line);
        std::printf("\n          Second geometry (derived from reference geometry)after optimization subject to \n");
        std::printf("          GEO_REF constraint of %.1f Kcal/mol/Angstrom^2 towards the data-set geometry written to file:\n", density);
        std::printf("          '%s'\n", line.c_str());
        geout(iarc);
      }
    }
  }
  std::printf("\n          Job name:  %39s\n", ("'" + input_fn.substr(0, input_fn.size() - 5) + "'").c_str());
  std::printf("          Current value of GEO_REF constraint:%20.2f  Kcal/mol/Angstrom^2\n", density);
  std::printf("          Distance between the geometries:%24.2f  Angstroms\n", sum);
  for (int i = 1; i <= nvar; ++i) {
    int kk = loc[1][i];
    int ll = loc[2][i];
    geo[ll][kk] = big_xparam[i];
  }
  geo_diff(sum, rms, true);
  std::printf("\n          Heat of formation of the first geometry:%17.3f Kcal/mol\n", escf1 - e_stress);
  std::printf("          Heat of formation of the second geometry:%16.3f Kcal/mol\n", escf2 - e_stress);
  std::printf("          Contribution to heat of formation due to stress:%8.2f  Kcal/mol\n", 2.0 * e_stress);
  escf1 = std::sqrt(dot1b(nvar, &big_grad[1], &big_grad[1]));
  escf2 = std::sqrt(dot1b(nvar, &big_grad[nvar + 1], &big_grad[nvar + 1]));
  std::printf("          Gradient arising from first geometry:%19.2f  Kcal/mol/Angstrom\n", escf1);
  std::printf("          Gradient arising from second geometry:%18.2f  Kcal/mol/Angstrom\n", escf2);
  sum = dot1b(nvar, &big_grad[1], &big_grad[nvar + 1]) / (escf1 * escf2);
  if (sum < 0.0)
    std::printf("          Angle between gradient vectors:%25.2f  degrees\n", std::acos(sum) * 57.2957795);
  else
    std::printf("          WARNING! - Angle between gradient vectors:%12.2f  degrees\n", std::acos(sum) * 57.2957795);
  for (int i = 1; i <= nvar; ++i) {
    xparam[i] = 0.5 * (big_xparam[i] + big_xparam[i + nvar]);
    int kk = loc[1][i];
    int ll = loc[2][i];
    geo[ll][kk] = xparam[i];
  }
  if (extra_print) {
    line = input_fn.substr(0, input_fn.size() - 5);
    char buf[128];
    std::snprintf(buf, sizeof(buf), " %.1f average.mop", density);
    line += buf;
    int ii = (int)line.size() - 13;
    line[ii] = 'p';
    add_path(line);
    std::printf("\n          Average of first and second geometries after optimization subject to \n");
    std::printf("          GEO_REF constraint of %.1f Kcal/mol/Angstrom^2 towards the data-set geometry written to file:\n", density);
    std::printf("          '%s'\n", line.c_str());
    geout(iarc);
  }
}

// ============================================================================
// compfg_TS — evaluate heat of formation (and gradients) of both systems.
// ============================================================================
void compfg_TS(const std::vector<double>& big_xparam, bool int_flag, double& escf1,
               double& escf2, bool fulscf, std::vector<double>& big_grad, bool lgrad) {
  std::vector<double> xparam((size_t)nvar + 1, 0.0);
  for (int i = 1; i <= nvar; ++i) xparam[i] = big_xparam[i];
  big_swap(1, 1);
  compfg(xparam, int_flag, escf1, fulscf, grad, lgrad);
  big_swap(0, 1);
  for (int i = 1; i <= nvar; ++i) big_grad[i] = grad[i];
  for (int j = 1; j <= 3; ++j)
    for (int i = 1; i <= numat; ++i) geoa[j][i] = geo[j][i];
  big_swap(1, 2);
  for (int i = 1; i <= nvar; ++i) xparam[i] = big_xparam[nvar + i];
  compfg(xparam, int_flag, escf2, fulscf, grad, lgrad);
  for (int i = 1; i <= nvar; ++i) big_grad[nvar + i] = grad[i];
  big_swap(0, 2);
  for (int j = 1; j <= 3; ++j)
    for (int i = 1; i <= numat; ++i) geo[j][i] = geo_1[j][i];
  for (int j = 1; j <= 3; ++j)
    for (int i = 1; i <= numat; ++i) geoa[j][i] = geo_2[j][i];
}

// ============================================================================
// get_pars — calibration reader.  F90 reads unit 33 (externally attached);
// C++ reads the file "calib.dat" when present, else leaves nloop = 0.
// ============================================================================
void get_pars(std::vector<double>& stresses, std::vector<double>& gradients,
              std::vector<double>& relscf, std::vector<double>& cutoff, int& nloop) {
  FILE* f = std::fopen("calib.dat", "r");
  if (!f) {
    nloop = 0;
    return;
  }
  char buf[512];
  if (!std::fgets(buf, sizeof(buf), f)) { std::fclose(f); nloop = 0; return; }
  nloop = 0;
  while (std::fgets(buf, sizeof(buf), f)) {
    ++nloop;
    double s = 0.0, g = 0.0, r = 0.0, c = 0.0;
    if (std::sscanf(buf, "%lf %lf %lf %lf", &s, &g, &r, &c) != 4) break;
    stresses[nloop] = s;
    gradients[nloop] = g;
    relscf[nloop] = r;
    cutoff[nloop] = c;
  }
  --nloop;
  std::fclose(f);
}

// ============================================================================
// build_active_site — identify active-site atoms from topology differences
// between the two systems (data-set vs. reference), then expand by nearest
// neighbors to form Set 1 and Set 2.
// ============================================================================
void build_active_site(std::vector<int>& active_site, std::vector<int>& ninsite) {
  for (size_t i = 0; i < active_site.size(); ++i) active_site[i] = 0;
  ninsite[1] = 0; ninsite[2] = 0; ninsite[3] = 0;
  if (keywrd.find(" LOCATE-TS(SET") != std::string::npos) {
    int nsite = 0;
    for (int i = 1; i <= nvar; ++i) {
      int j = 1;
      for (; j <= nsite; ++j)
        if (loc[1][i] == active_site[j]) break;
      if (j > nsite) {
        nsite = j;
        active_site[nsite] = loc[1][i];
      }
    }
    ninsite[1] = nsite;
    return;
  }
  std::vector<int> nbonds_b((size_t)numat + 1, 0);
  std::vector<std::vector<int>> ibonds_b(16, std::vector<int>((size_t)numat + 1, 0));
  big_swap(1, 1);  // Extract system 1 - the input data set
  for (int i = 1; i <= numat; ++i) nbonds_b[i] = nbonds[i];
  for (int j = 1; j <= 15; ++j)
    for (int i = 1; i <= numat; ++i) ibonds_b[j][i] = ibonds[j][i];
  big_swap(1, 2);  // Extract system 2 - the reference geometry
  set_up_dentate();
  big_swap(0, 2);  // Store system 2 (adds nbonds/ibonds to the store)

  int nsite = 0;
  for (size_t i = 0; i < active_site.size(); ++i) active_site[i] = 0;
  for (int i = 1; i <= numat; ++i) {
    if (nbonds[i] != nbonds_b[i]) {
      int k = 1;
      for (; k <= nsite; ++k)
        if (active_site[k] == i) break;
      if (k <= nsite) continue;
      ++nsite;
      active_site[nsite] = i;
    } else {
      for (int j = 1; j <= nbonds[i]; ++j) {
        if (ibonds[j][i] != ibonds_b[j][i]) {
          int k = 1;
          for (; k <= nsite; ++k)
            if (active_site[k] == i) break;
          if (k <= nsite) continue;
          ++nsite;
          active_site[nsite] = i;
        }
      }
    }
  }
  if (nsite == 0) {
    mopend(" No atoms involved in bond-making or bond-breaking!");
    return;
  }
  int k = 0;
  for (int i = 1; i <= nsite; ++i) k = std::max(k, nbonds_b[active_site[i]]);
  char ch[8];
  std::snprintf(ch, sizeof(ch), "%d", k * 6 - 6);
  if (maxtxt == 26) {
    std::printf("\n%6s Set 1: Atoms involved in covalent bond-breaking and bond-making%s\n",
                "", "Connectivity of atoms");
    std::printf("%77s%*s%s\n", "Data-set", k * 6, "", "Reference");
  } else {
    std::printf("\n%6s Set 1: Atoms involved in covalent%s\n", "", "Connectivity of atoms");
    std::printf("%7s%s\n", "bond-breaking and bond-making            Data-set", "Reference");
  }
  for (int i = 1; i <= nsite; ++i) {
    char lbuf[256];
    int pos = 0;
    for (int j = 1; j <= nbonds_b[active_site[i]]; ++j)
      pos += std::snprintf(lbuf + pos, sizeof(lbuf) - pos, "%6d", ibonds_b[j][active_site[i]]);
    std::string lstr(lbuf, (size_t)pos);
    lstr.resize((size_t)k * 6 + 9, ' ');
    char rbuf[256];
    int pos2 = 0;
    for (int j = 1; j <= nbonds[active_site[i]]; ++j)
      pos2 += std::snprintf(rbuf + pos2, sizeof(rbuf) - pos2, "%6d", ibonds[j][active_site[i]]);
    std::string rpart(rbuf, (size_t)pos2);
    lstr.replace((size_t)k * 6 + 9, rpart.size(), rpart);
    if (maxtxt == 26) {
      std::printf("%4d   Atom:  %s%6d  PDB Label: (%s)%s\n",
                  i, elemnt[nat[active_site[i]]].c_str(), active_site[i],
                  txtatm[active_site[i]].c_str(), lstr.c_str());
    } else {
      std::printf("%12d   Atom:  %s%6d      %s\n",
                  i, elemnt[nat[active_site[i]]].c_str(), active_site[i], lstr.c_str());
    }
  }
  ninsite[1] = nsite;

  // Select neighbors (Set 2).
  int ii = nsite;
  for (int jj = 1; jj <= ii; ++jj) {
    int i = active_site[jj];
    for (int j = 1; j <= nbonds[i]; ++j) {
      int kk = ibonds[j][i];
      int l = 1;
      for (; l <= nsite; ++l)
        if (active_site[l] == kk) break;
      if (l <= nsite) continue;
      ++nsite;
      active_site[nsite] = kk;
    }
  }
  for (int jj = 1; jj <= ii; ++jj) {
    int i = active_site[jj];
    for (int j = 1; j <= nbonds_b[i]; ++j) {
      int kk = ibonds_b[j][i];
      int l = 1;
      for (; l <= nsite; ++l)
        if (active_site[l] == kk) break;
      if (l <= nsite) continue;
      ++nsite;
      active_site[nsite] = kk;
    }
  }
  std::printf("\n Set 2 : Atoms involved in covalent bond-breaking and bond-making, plus nearest neighbors\n");
  if (maxtxt == 26) {
    for (int i = 1; i <= nsite; ++i)
      std::printf("%4d   Atom:  %s%6d  PDB Label: (%s)\n",
                  i, elemnt[nat[active_site[i]]].c_str(), active_site[i],
                  txtatm[active_site[i]].c_str());
  } else {
    for (int i = 1; i <= nsite; ++i)
      std::printf("%12d   Atom:  %s%6d\n",
                  i, elemnt[nat[active_site[i]]].c_str(), active_site[i]);
  }
  ninsite[2] = nsite;
}

// ============================================================================
// select_opt — sort the first shell atoms of the active site ascending.
// The F90 "if (.false.)" block (variable re-selection) is dead code and is
// omitted; the sorting that feeds the active-site variable set is kept.
// ============================================================================
void select_opt(int shell, std::vector<int>& active_site) {
  std::vector<int> use((size_t)shell + 1, 0), tmp((size_t)shell + 1, 0);
  for (int i = 1; i <= shell; ++i) tmp[i] = active_site[i];
  for (int i = 1; i <= shell; ++i) {
    int k = 100000;
    int l = 1;
    for (int j = 1; j <= shell; ++j) {
      if (tmp[j] < k) {
        k = tmp[j];
        l = j;
      }
    }
    tmp[l] = 200000;
    use[i] = k;
  }
  for (int i = 1; i <= shell; ++i) active_site[i] = use[i];
}

// ============================================================================
// Refine_TS_for_Proteins — iterative transition-state refinement: energy
// minimization of everything except the active site, then gradient
// minimization of the active site, up to five cycles.
// ============================================================================
void Refine_TS_for_Proteins() {
  int i = (int)keywrd.find(" LOCATE-TS");
  if (i < 0) i = (int)keywrd.find(" REFINE-TS");
  if (i < 0) i = 0;
  int j = i + 12;
  while (j <= std::min(240, i + 100)) {
    if (keywrd[j] == ' ') break;
    ++j;
  }
  std::string region = (size_t)i < keywrd.size() ? keywrd.substr(i, (size_t)(j - i)) : "";
  bool extra_print = (region.find("C:") != std::string::npos);
  char num = (char)('1' + (int)std::log10(std::max(1.0, j * 1.01)));
  bool converged = false;
  double gnorm_lim = 0.0;

  for (int loop = 1; loop <= 5; ++loop) {
    if (loop > 1) {
      for (int i2 = 1; i2 <= numat; ++i2)
        for (int j2 = 1; j2 <= 3; ++j2) lopt[j2][i2] = 1 - lopt[j2][i2];
    }
    loc[1].assign(loc[1].size(), 0);
    loc[2].assign(loc[2].size(), 0);
    nvar = 0;
    for (int i2 = 1; i2 <= numat; ++i2) {
      for (int j2 = 1; j2 <= 3; ++j2) {
        if (lopt[j2][i2] == 1) {
          ++nvar;
          loc[1][nvar] = i2;
          loc[2][nvar] = j2;
          xparam[nvar] = geo[j2][i2];
        }
      }
    }
    l_control("TS", 2, -1);
    std::printf("  Loop:%4d  Energy minimization, excluding active site\n", loop);
    gnorm_lim = nint(std::pow(numat, 0.25) * 2.0 + 1.0);
    {
      char lbuf[64];
      std::snprintf(lbuf, sizeof(lbuf), "LET DDMIN=0 GNORM=%.1f", gnorm_lim);
      std::string ltmp = lbuf;
      l_control(ltmp, (int)ltmp.size(), 1);
    }
    ++numcal;
    if (nvar > 0) {
      lbfgs(&xparam[1], escf);
    } else {
      compfg(xparam, true, escf, true, grad, false);
      gnorm = 0.0;
    }
    if (gnorm < gnorm_lim) {
      int ki = (int)keywrd.find(" GNORM");
      double cut = (ki >= 0) ? reada(keywrd, ki + 6) : gnorm_lim;
      std::printf("\n      GRADIENT =%9.5f  IS LESS THAN CUTOFF =%9.5f\n\n", gnorm, cut);
    }
    converged = (nstep < 3);
    if (moperr) return;
    for (int i2 = 1; i2 <= numat; ++i2)
      for (int j2 = 1; j2 <= 3; ++j2) lopt[j2][i2] = 1 - lopt[j2][i2];
    loc[1].assign(loc[1].size(), 0);
    loc[2].assign(loc[2].size(), 0);
    nvar = 0;
    for (int i2 = 1; i2 <= numat; ++i2) {
      for (int j2 = 1; j2 <= 3; ++j2) {
        if (lopt[j2][i2] == 1) {
          ++nvar;
          loc[1][nvar] = i2;
          loc[2][nvar] = j2;
          xparam[nvar] = geo[j2][i2];
        }
      }
    }
    std::printf("  Loop:%4d  Gradient minimization of atoms in the active site\n", loop);
    line = input_fn.substr(0, input_fn.size() - 5);
    char numbuf[16];
    std::snprintf(numbuf, sizeof(numbuf), " Loop%d.mop", loop);
    line += numbuf;
    if (extra_print) {
      add_path(line);
      std::printf("\n          Transition state on cycle%2d written to file:\n", loop);
      std::printf("          '%s'\n", line.c_str());
      geout(iarc);
    }
    l_control("TS", 2, 1);
    l_control("GNORM=3", 7, 1);
    if (loop == 2) l_control("OLD_HESS", 8, 1);
    l_control("OLD_SCF", 7, 1);
    for (int j2 = 1; j2 <= 3; ++j2)
      for (int i2 = 1; i2 <= numat; ++i2) geoa[j2][i2] = geo[j2][i2];
    ef(xparam, escf);
    if (moperr && loop > 1) {
      std::printf("\n          Gradient minimization failed.  An attempt will be made to correct the error\n");
      moperr = false;
      ++numcal;
      l_control("OLD_HESS", 8, -1);
      l_control("OLD_SCF", 7, -1);
      for (int j2 = 1; j2 <= 3; ++j2)
        for (int i2 = 1; i2 <= numat; ++i2) geo[j2][i2] = geoa[j2][i2];
      for (int i2 = 1; i2 <= nvar; ++i2)
        xparam[i2] = geo[loc[2][i2]][loc[1][i2]];
      ef(xparam, escf);
      if (moperr) {
        std::printf("\n          Gradient minimization failed.  Best geometry at this point will be output\n");
        for (int j2 = 1; j2 <= 3; ++j2)
          for (int i2 = 1; i2 <= numat; ++i2) geo[j2][i2] = geoa[j2][i2];
        moperr = false;
        converged = true;
        return;
      } else {
        std::printf("\n          Gradient minimization succeeded.  Error removed.  Job continuing\n");
        l_control("OLD_HESS", 8, 1);
        l_control("OLD_SCF", 7, 1);
      }
    }
    converged = (converged && nstep < 3);
    if (converged) break;
  }
}
