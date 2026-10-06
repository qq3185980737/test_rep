// mod_calls_cublas.h — C++ translation.
#pragma once
extern "C" {
void call_asum_cublas(int n, const double* vecx, int incx, double* res);
}
