#include <cstdio>
#include <vector>
#include "dihed.h"
int main(){
  std::vector<std::vector<double>> xyz(4, std::vector<double>(5, 0.0));
  double a;
  xyz[1][1]=1;xyz[2][1]=0;xyz[3][1]=0; xyz[1][2]=0;xyz[2][2]=0;xyz[3][2]=0;
  xyz[1][3]=0;xyz[2][3]=1;xyz[3][3]=0; xyz[1][4]=0;xyz[2][4]=0;xyz[3][4]=1;
  dihed(xyz,1,2,3,4,a); std::printf("A = %.12f\n", a);
  xyz[1][1]=0;xyz[2][1]=0;xyz[3][1]=0; xyz[1][2]=1;xyz[2][2]=0;xyz[3][2]=0;
  xyz[1][3]=1;xyz[2][3]=1;xyz[3][3]=0; xyz[1][4]=1;xyz[2][4]=1;xyz[3][4]=1;
  dihed(xyz,1,2,3,4,a); std::printf("C = %.12f\n", a);
  return 0;
}
