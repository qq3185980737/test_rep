// compct.cpp — C++ translation of MOPAC 2016 "compct.F90".

#include "compct.h"

void compct(std::vector<int>& nncnew, std::vector<int>& ncnew,
            std::vector<int>& ncmnew, int itop, std::vector<int>& nc,
            std::vector<int>& ic, std::vector<int>& iws, int n01,
            std::vector<double>& c, int n02, int nmos, int idone,
            int& lb, int& mb, int icref, int inref) {
    (void)n01; (void)n02;
    int inic = icref, inco = inref;
    int inew = itop - 1;
    int ii = idone;
    for (int i = idone - 1; i >= 1; --i) {
        if (nc[i] != 0) {
            ii--;
            inew++;
            int natom = nc[i], ncoef = iws[i];
            inic -= natom;
            inco -= ncoef;
            int ioic = nncnew[inew], ioco = ncmnew[inew];
            for (int n = natom; n >= 1; --n) ic[inic + n] = ic[ioic + n];
            for (int n = ncoef; n >= 1; --n) c[inco + n] = c[ioco + n];
            ncnew[inew] = natom;
            nncnew[inew] = inic;
            ncmnew[inew] = inco;
            nc[ii] = natom;
            iws[ii] = ncoef;
            if (inew == nmos) break;
        }
    }
    for (int i = ii - 1; i >= 1; --i) nc[i] = 0;
    lb = inic;
    mb = inco;
}
