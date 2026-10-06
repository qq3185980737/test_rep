// deri2.h — C++ translation of MOPAC 2016 "deri2.F90".
#pragma once
#include <vector>

// Relaxation part of non-variational energy derivatives (skeleton; external
// linear-algebra / CI helpers not yet ported).
void deri2(int minear,
           std::vector<std::vector<double>>& f,
           std::vector<std::vector<double>>& fd,
           std::vector<std::vector<double>>& fci,
           int ninear, int nvar_nvo,
           std::vector<double>& dxyzr, double throld,
           std::vector<double>& diag, std::vector<double>& scalar,
           std::vector<double>& work);
