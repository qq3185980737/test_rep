#include <cstdio>
#include "eigenvectors_LAPACK.h"
int main(){
  double xmat[6] = {3.0, 0.0, 1.0, 0.0, 0.0, 2.0};
  double vecs[9] = {}, vals[3] = {};
  eigenvectors_LAPACK(vecs, xmat, vals, 3);
  for (int c=0;c<3;++c){ printf("col%d: ", c); for(int r=0;r<3;++r) printf("%.6f ", vecs[r+c*3]); printf("\n"); }
  return 0;
}
