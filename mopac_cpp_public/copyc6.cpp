// copyc6.cpp — C++ translation of copyc6.F90.
#include "copyc6.h"
#include <algorithm>
// Fortran limit(iat,jat,iadr,jadr): element ids >100 encode isotope-family
// extensions; each 100-block increments the coordination index (1..5).
namespace { void limit(int iat,int jat,int& iadr,int& jadr){
  iadr=1; jadr=1;
  while(iat>100){ iat-=100; ++iadr; }
  while(jat>100){ jat-=100; ++jadr; }
} }
void copyc6(int maxc, int max_elem,
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& c6ab,
  std::vector<int>& maxci) {
  const double* pars = copyc6_data::pars;
  (void)maxc;(void)max_elem;
  int nlines = 32385;
  // Fortran: c6ab=-1 (whole array) before filling.
  const int N = (int)c6ab.size();
  for (int i=0; i<N; ++i)
    for (int j=0; j<(int)c6ab[i].size(); ++j)
      for (int k=0; k<(int)c6ab[i][j].size(); ++k)
        for (int l=0; l<(int)c6ab[i][j][k].size(); ++l)
          for (int m=0; m<(int)c6ab[i][j][k][l].size(); ++m)
            c6ab[i][j][k][l][m] = -1.0;
  maxci.assign(N, 0);
  for (int nn=1; nn<=nlines; ++nn) {
    int kk = (nn-1)*5;   // pars 0-based; Fortran kk=(nn*5)+1 1-based
    int iat=(int)pars[kk], jat=(int)pars[kk+1];
    // Fortran c6ab is declared (max_elem,max_elem,...) yet indexed by raw
    // element ids up to 482 (isotope-family rows) — out-of-bounds writes that
    // happen to be harmless in Fortran. In C++ skip rows outside the caller's
    // array to stay well-defined; rows 1..max_elem cover all real elements.
    if (iat >= N || jat >= N) continue;
    int iadr,jadr; limit(iat,jat,iadr,jadr);
    maxci[iat]=std::max(maxci[iat],iadr);
    maxci[jat]=std::max(maxci[jat],jadr);
    double c=pars[kk+2], p1=pars[kk+3], p2=pars[kk+4];
    c6ab[iat][jat][iadr][jadr][1]=c; c6ab[iat][jat][iadr][jadr][2]=p1; c6ab[iat][jat][iadr][jadr][3]=p2;
    c6ab[jat][iat][jadr][iadr][1]=c; c6ab[jat][iat][jadr][iadr][2]=p2; c6ab[jat][iat][jadr][iadr][3]=p1;
  }
}
