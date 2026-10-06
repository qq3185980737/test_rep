// sqrdc_sqrsl.h — C++ translation of LINPACK SQRDC / SQRSL (float, 1-based,
// column-major storage x[(j-1)*ldx + (i-1)] = element (i,j)).
#pragma once
#include <vector>

void sqrdc(std::vector<float>& x, int ldx, int n, int p,
           std::vector<float>& qraux, std::vector<int>& jpvt,
           std::vector<float>& work, int job);
void sqrsl(std::vector<float>& x, int ldx, int n, int k,
           std::vector<float>& qraux, std::vector<float>& y,
           std::vector<float>& qy, std::vector<float>& qty,
           std::vector<float>& b, std::vector<float>& rsd,
           std::vector<float>& xb, int job, int& info);
