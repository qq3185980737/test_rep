// H_bonds4.h — C++ translation.
#pragma once
#include <vector>
double H_bonds4(bool l_grad, double* dxyz);
double energy_corr_hh_rep(bool l_grad, double* dxyz);
double energy_corr_h4(bool l_grad, double* grad_h4);
double poly(double r, bool l_grad, double& dpoly);
double cvalence_contribution(int atom_a, int atom_b);
double cvalence_contribution_d(int atom_a, int atom_b);
void prt_hbonds(int D, int H, int A, double energy);
