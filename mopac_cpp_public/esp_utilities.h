// esp_utilities.h — C++ translation of MOPAC 2016 "esp_utilities.F90".
#pragma once
#include <vector>

void rys(double x, int nroots, double u[4], double w[4]);
void evec(std::vector<float>& aVector, double x, double y, double z,
          const std::vector<std::vector<double>>& coord, int numat);
void vint(double& xint, double& yint, double& zint, int ni, int nj,
          double x0, double y0, double z0,
          double xi, double yi, double zi,
          double xj, double yj, double zj, double t);

// BLAS (float)
void saxpy(int n, float sa, const std::vector<float>& sx, int incx,
           std::vector<float>& sy, int incy);
float sdot(int n, const std::vector<float>& sx, int incx,
           const std::vector<float>& sy, int incy);
void sscal(int n, float da, std::vector<float>& dx, int incx);
float snrm2(int n, const std::vector<float>& sx, int incx);
void sswap(int n, std::vector<float>& sx, int incx,
           std::vector<float>& sy, int incy);
void scopy(int n, const std::vector<float>& sx, int incx,
           std::vector<float>& sy, int incy);
