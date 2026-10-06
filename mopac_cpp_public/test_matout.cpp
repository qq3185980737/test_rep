// test_matout.cpp
#include <cstdio>
#include "matout.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
  chanel_C::iw = 6;
  molkst_C::numat = 1;
  common_arrays_C::nfirst.resize(2); common_arrays_C::nfirst[1] = 1;
  common_arrays_C::nlast.resize(2);  common_arrays_C::nlast[1] = 5;
  common_arrays_C::nat.resize(2);    common_arrays_C::nat[1] = 1;
  double a[10] = {1.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  double b[5] = {0.0, -12.3, -4.56, 0.0, 0.0};
  int nr = 5;
  matout(a, b, 2, nr, 5);
  std::printf("matout nr=%d PASS\n", nr);
  return 0;
}
