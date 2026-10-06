// ef.h — C++ translation of MOPAC 2016 "ef.F90".
#pragma once
#include <vector>

// P-RFO/QA optimizer driver (ef.F90 "ef").
void ef(std::vector<double>& xparam, double& funct);
// Hessian update: iupd=1 Powell, iupd=2 BFGS (ef.F90 "updhes").
void updhes(std::vector<double>& svec, std::vector<double>& tvec,
            const std::vector<double>& grad, int nvar, int iupd,
            std::vector<double>& hess, const std::vector<double>& oldf,
            const std::vector<double>& d);
