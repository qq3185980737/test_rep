#include <cstdio>
#include "partxy.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
int main(){
    molkst_C::lm61=8; molkst_C::numat=1; molkst_C::numcal=3;
    common_arrays_C::nfirst.assign({0,1}); common_arrays_C::nlast.assign({0,1});
    double c34[8]={0,2.0,0,0,0,0,0,0};
    double pq34[8]={};
    double w[16]={0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    partxy(c34,pq34,w);
    std::printf("single: pq34[0]=%.6f pq34[1]=%.6f\n", pq34[0], pq34[1]);
    return 0;
}
