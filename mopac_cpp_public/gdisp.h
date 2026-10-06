// gdisp.h — C++ translation of MOPAC 2016 "gdisp.F90".
#pragma once
#include <vector>
#include "dftd3_bits.h"

void ncoord(int natoms, const std::vector<double>& rcov,
            const std::vector<int>& nat,
            const std::vector<std::vector<double>>& xyz,
            std::vector<double>& cn);

void get_dC6_dCNij(int maxc, int max_elem, const c6ab_t& c6ab,
                   int mxci, int mxcj, double cni, double cnj,
                   int izi, int izj, double& c6check, double& dc6i, double& dc6j);

void gdisp(const std::vector<std::vector<double>>& xyz,
           const std::vector<std::vector<double>>& r0ab,
           double rs6, double alp6, const c6ab_t& c6ab, double s6,
           const std::vector<int>& mxc,
           const std::vector<double>& rcov,
           std::vector<std::vector<double>>& dxyz_temp);
