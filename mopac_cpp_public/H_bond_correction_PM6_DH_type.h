// H_bond_correction_PM6_DH_type.h — C++ translation.
#pragma once
double PM6_DH_H_bond_corrections(bool l_grad, bool prt);
void setup_DH_Plus(int nrpairs, int* nrbondsa, int* nrbondsb,
                   bool& l_h_bonds, double* covrad);
