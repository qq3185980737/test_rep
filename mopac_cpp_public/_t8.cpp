#include <cstdio>
#include "partxy.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
int main(){
    molkst_C::lm61=8; molkst_C::numat=2; molkst_C::numcal=1;
    common_arrays_C::nfirst.assign({0,1,2}); common_arrays_C::nlast.assign({0,1,2});
    double c34[8]={0,0.5,1.5,0,0,0,0,0};
    double pq34[8]={};
    double w[16]={0,0.4,0.3,0.2,0,0,0,0,0,0,0,0,0,0,0,0};
    partxy(c34,pq34,w);
    std::printf("pq34[0..2]: %.6f %.6f %.6f\n", pq34[0], pq34[1], pq34[2]);
    return 0;
}
