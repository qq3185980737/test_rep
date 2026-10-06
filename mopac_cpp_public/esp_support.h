// esp_support.h — helpers for new_esp (MOPAC 2016).
#pragma once
#include <vector>

// S^{-1/2} helper — direct translation of esp_utilities.F90
// "get_minus_point_five_overlap" (rebuilds the normalised overlap from the
// packed h matrix, diagonalises it, forms S^{-1/2}).  Fills s (1-based
// matrix) with S^{-1/2}.
void get_minus_point_five_overlap(std::vector<std::vector<double>>& s);

// Density-matrix construction from MO coefficients (density_for_GPU, GPU
// variant of densit): p(ij) = 2*C_closed*C_closed^T + fract*C_open*C_open^T.
void density_for_GPU(const double* vecs, double fract, int nclose, int nopen,
                     double two, int ij, int norbs, int mode,
                     std::vector<double>& p, int mode2);
