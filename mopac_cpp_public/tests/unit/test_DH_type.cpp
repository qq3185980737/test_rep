// test_DH_type.cpp
#include <cstdio>
#include "H_bond_correction_PM6_DH_type.h"
int main() { bool l; double v=PM6_DH_H_bond_corrections(false,false); double c[94]; setup_DH_Plus(0,nullptr,nullptr,l,c); std::printf("DH_type links=%g PASS\n",v); return 0; }
