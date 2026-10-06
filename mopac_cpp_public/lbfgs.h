// lbfgs.h — C++ translation of MOPAC 2016 "lbfgs.F90".
#pragma once
#include <string>
#include <vector>

// Driver used by run_mopac dispatch ("LBFGS"). xparam is 1-based (size nvar+1).
void lbfgs(double* xparam, double& escf);

// L-BFGS-B restart-file save (mode=1) / restore (mode=0).
void lbfsav(double tt0, int mode, std::vector<double>& wa, int nwa,
            std::vector<int>& iwa, int niwa, std::string& task,
            std::string& csave, std::vector<int>& lsave, std::vector<int>& isave,
            std::vector<double>& dsave, int& nstep, double& escf);

// L-BFGS-B entry point (called by lbfgs_TS). All arrays are 1-based.
void setulb(int n, int m, double* x, double* l, double* u, int* nbd,
            double& f, double* g, double factr, double pgtol,
            std::vector<double>& wa, std::vector<int>& iwa, std::string& task,
            int iprint, std::string& csave, std::vector<int>& lsave,
            std::vector<int>& isave, std::vector<double>& dsave);
