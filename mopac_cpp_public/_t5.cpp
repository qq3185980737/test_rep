#include <cstdio>
#include <vector>
#include "ijkl.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
namespace cosmo_C { extern bool iseps; }
int main(){
    molkst_C::norbs=2; molkst_C::numat=2; molkst_C::lm61=8;
    common_arrays_C::nfirst.assign({0,1,2}); common_arrays_C::nlast.assign({0,1,2});
    common_arrays_C::w.assign(16,0.0);
    common_arrays_C::w[1]=0.4; common_arrays_C::w[2]=0.3; common_arrays_C::w[3]=0.2;
    cosmo_C::iseps=false;
    double cp[2]={0.6,0.8}; double cf[4]={0.6,0.8,0.8,-0.6};
    std::vector<double> dijkl(2*1*1,0.0);
    double cij[8]={}, ckl[8]={}, wcij[8]={}; double xy[1]={0.0};
    ijkl(cp,cf,0,1,dijkl.data(),cij,ckl,wcij,xy);
    std::printf("cij: %f %f\n", cij[1], cij[2]);
    std::printf("wcij: %f %f\n", wcij[1], wcij[2]);
    std::printf("dijkl: %f %f\n", dijkl[0], dijkl[1]);
    std::printf("xy: %f\n", xy[0]);
    return 0;
}
