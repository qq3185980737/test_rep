// linpack.h — C++ translation.
#pragma once
void dgedi(double* a, int lda, int n, const int* ipvt, double* det,
           double* work, int job);
void dgefa(double* a, int lda, int n, int* ipvt, int& info);
