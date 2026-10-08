#include <cstdio>
#include <vector>
#include "partxy.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
int main(){
    molkst_C::lm61=8; molkst_C::numat=2; molkst_C::numcal=1;
    common_arrays_C::nfirst.assign({0,1,2}); common_arrays_C::nlast.assign({0,1,2});
    // direct partxy call replicating ijkl cij
    double cij[8]={0,0.72,1.28,0,0,0,0,0};
    double wcij[8]={};
    double w[16]={0,0.4,0.3,0.2,0,0,0,0,0,0,0,0,0,0,0,0};
    partxy(cij,wcij,w);
    std::printf("direct wcij: %f %f\n", wcij[1], wcij[2]);
    return 0;
}
