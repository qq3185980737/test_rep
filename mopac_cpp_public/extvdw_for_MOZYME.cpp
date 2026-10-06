// extvdw_for_MOZYME.cpp — C++ translation of "extvdw_for_MOZYME"
// (geochk.F90 lines 1490-1579).  Fills radius(1..numat) from the reference
// table refvdw, honoring "METAL" and "VDWM(:X=value:...)" keyword overrides.
#include "extvdw_for_MOZYME.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "mod_atomradii.h"
#include "molkst_C.h"
#include "mopend.h"
#include "reada.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace mod_atomradii;
using namespace molkst_C;

// 1-based Fortran index.
static int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// Fortran len_trim: length without trailing blanks.
static std::size_t len_trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? 0 : e + 1;
}

static int nint(double x) { return static_cast<int>(std::floor(x + 0.5)); }

void extvdw_for_MOZYME(std::vector<double>& radius,
                       const std::vector<double>& refvdw) {
  std::vector<double> vdw(107, 0.0);  // 1-based, 106 entries
  std::string txt_rad(80, ' ');
  int i, j, k;

  // Modify radii of various atoms.  Keyword: METAL(...)
  // e.g. METAL(:H=1.0:Cl=1.7)
  i = index1(keywrd, " METAL");
  if (i != 0) {
    j = index1(keywrd.substr(i - 1), ")") + i;   // global end of METAL(...)
    txt_rad = keywrd.substr(i - 1, j - i + 1);   // keywrd(i:j)
    for (k = 1; k <= 83; ++k) {
      if (k >= (int)is_metal.size()) break;
      if (!is_metal[k])
        is_metal[k] = (index1(txt_rad, cap_elemnt[k]) != 0);
    }
  }

  // Keyword: VDWM(:symbol=value:...) e.g. VDWM(:H=1.0:Cl=1.7)
  i = index1(keywrd, " VDWM(");
  std::string line;
  if (i == 0) {
    line = " ";
  } else {
    // Normalize separators to ';' inside the VDWM(...) field.
    int p6 = i + 5;  // 0-based position of character after "VDWM("
    std::string& kw = keywrd;
    // keywrd is a std::string (not fixed length); guard against short keys.
    if ((int)kw.size() < p6) { line = " "; }
    else {
      if (kw[p6] != ';' && kw[p6] != ':')
        kw.insert(p6, 1, ';');
      j = index1(kw.substr(p6), ")") + p6 + 1;  // global 1-based end of VDWM(...)
      for (k = p6; k <= j; ++k) {               // normalize separators
        if (k >= (int)kw.size()) break;
        if (kw[k] == ':') kw[k] = ';';
        if (kw[k] == ',') kw[k] = ';';
      }
      line = kw.substr(p6, j - p6);             // "(:H=1.0;Cl=1.7)" per Fortran slice
    }
  }

  vdw = refvdw;   // vdw(:106) = refvdw(:106)
  if (line != " ") {
    for (i = 1; i <= 106; ++i) {
      j = 2;
      if ((int)cap_elemnt[i].size() >= 2 && cap_elemnt[i][1] == ' ') j = 1;
      std::string key = ";" + cap_elemnt[i].substr(0, j) + "=";
      k = index1(line, key);
      if (k > 0) vdw[i] = reada(line, k);
    }
  }

  // Verify that all radii that will be used are, in fact, set correctly.
  for (i = 1; i <= numat; ++i) {
    if (nat[i] > 102) continue;
    if (vdw[nat[i]] > 900.0) {
      line = "MISSING VAN DER WAALS RADIUS FOR " + cap_elemnt[nat[i]];
      mopend(line.substr(0, len_trim(line)));
      j = 2;
      if ((int)cap_elemnt[i].size() >= 2 && cap_elemnt[i][1] == ' ') j = 1;
      std::fprintf(stdout, "To correct this, add keyword 'VDWM(:%s=n.nn)'\n",
                   cap_elemnt[nat[i]].substr(0, j).c_str());
      return;
    }
  }

  // Flag which elements are to be treated as metals by giving them
  // negative atomic radii.
  for (i = 1; i <= 102; ++i) {
    if (i < (int)is_metal.size() && is_metal[i]) vdw[i] = -1.0;
  }
  for (i = 1; i <= numat; ++i) radius[i] = vdw[nat[i]];

  // Set VDW radius for individual atoms to -1 (from numeric METAL list).
  if (len_trim(txt_rad) != 0) {
    for (i = 1; i <= (int)len_trim(txt_rad); ++i) {
      int d = (unsigned char)txt_rad[i - 1] - (unsigned char)'0';
      if (d <= 9 && d > 0) {
        j = nint(reada(txt_rad.substr(i - 1), 1));
        radius[j] = -1.0;
        for (j = i; j <= (int)len_trim(txt_rad); ++j) {
          int d2 = (unsigned char)txt_rad[i - 1] - (unsigned char)'0';
          if (d2 > 9 || d2 < 0) break;
          txt_rad[i - 1] = ' ';
        }
      }
    }
  }
}
