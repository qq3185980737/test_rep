// test_eig_lapack.cpp
#include <cstdio>
#include "eigenvectors_LAPACK.h"
int main() { double v[4],x[3],e[2]; eigenvectors_LAPACK(v,x,e,2); std::printf("eig_lapack links OK PASS\n"); return 0; }
