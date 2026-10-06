// deri21.h — C++ translation of MOPAC 2016 "deri21.F90".
#pragma once
// Least-squares (PCA) orthonormal basis extraction from set {A}.
// a: minear x nvar_nvo, column-major 0-based (F90 a(minear,nvar_nvo));
// vnert: nvar_nvo; pnert: nvar_nvo (sqrt of |eigenvalues|, ascending);
// b: minear x ncut, column-major 0-based output.
void deri21(const double* a, int nvar_nvo, int minear, int ifirst,
            double* vnert, double* pnert, double* b, int& ncut);
