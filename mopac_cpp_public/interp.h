// interp.h — C++ translation.
#pragma once
void interp(int np, int nq, int& mode, double e, double* fp, double* cp,
            double* theta, double* vec_interp, double* fock_interp,
            double* p_interp, double* h_interp, double* vecl, double& eold_ref);
void spline(double* x, double* f, double* df, double xhigh, double xlow,
            double& xmin, int npnts);
