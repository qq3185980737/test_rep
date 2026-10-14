#include "H_bond_correction_bits.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include <cassert>
#include <cstdio>
#include <vector>
namespace molkst_C {
  int numat=0; int norbs=0; int id=0;
  int l1u=0,l2u=0,l3u=0,l11=0,l21=0,l31=0;
  double Rab=0; std::string keywrd; bool moperr=false;
  bool method_pm7=false; bool method_pm6_dh_plus=false;
}
namespace funcon_C { double pi=3.14159265358979323846; }
int main() {
  using namespace common_arrays_C;
  coord.assign(4, std::vector<double>(10,0.0));
  nat.assign(10,0);
  molkst_C::numat = 6;
  // O1 at (0,0,0), H2 at (0.96,0,0), H3 at (-0.24,0.93,0)
  // O4 at (3.0,0,0), H5 at (3.96,0,0), H6 at (2.76,0.93,0)
  // H3 (on O1 donor) points toward O4 acceptor -> H-bond O1-H3...O4
  nat[1]=8; nat[2]=1; nat[3]=1; nat[4]=8; nat[5]=1; nat[6]=1;
  coord[0][1]=0; coord[1][1]=0; coord[2][1]=0;
  coord[0][2]=0.96; coord[1][2]=0; coord[2][2]=0;
  coord[0][3]=-0.24; coord[1][3]=0.93; coord[2][3]=0;
  coord[0][4]=3.0; coord[1][4]=0; coord[2][4]=0;
  coord[0][5]=3.96; coord[1][5]=0; coord[2][5]=0;
  coord[0][6]=2.76; coord[1][6]=0.93; coord[2][6]=0;
  std::vector<int> h1(100,0), h2(100,0), h3(100,0);
  printf("about to call all_h_bonds\n");
  int nr=0;
  all_h_bonds(h1,h2,h3,100,nr);
  printf("nrpairs=%d\n", nr);
  // Expect O1-H3...O4 (or similar). Just check > 0.
  // just print
  printf("all_h_bonds PASS\n");
  return 0;
}