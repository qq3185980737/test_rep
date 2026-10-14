// test_frame.cpp
#include <cstdio>
#include <vector>
#include "frame.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
extern "C" void axis_(double&, double&, double&, double* rot) {
    rot[0]=1; rot[4]=1; rot[8]=1;
    rot[1]=0; rot[2]=0; rot[3]=0;
    rot[5]=0; rot[6]=0; rot[7]=0;
}
int main() {
    molkst_C::numat = 1;
    common_arrays_C::coord.assign(4, std::vector<double>(2, 0.0));
    common_arrays_C::atmass.assign(2, 1.0);
    int lin = 3*4/2;
    std::vector<double> fmat(lin+1, 0.0);
    frame(fmat, molkst_C::numat, 1);
    bool ok = (fmat[1] == 41000.0);
    std::printf("frame fmat[1]=%g %s\n", fmat[1], ok?"PASS":"FAIL");
    return ok?0:1;
}
