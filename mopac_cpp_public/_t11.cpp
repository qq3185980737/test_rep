#include <cstdio>
#include <vector>
#include "dijkl1.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
int main(){
    meci_C::nmos=2; molkst_C::norbs=2; molkst_C::numat=2;
    common_arrays_C::nfirst.assign({0,1,2}); common_arrays_C::nlast.assign({0,1,2});
    std::vector<std::vector<double>> cv; cv.push_back({0.6,0.8}); cv.push_back({0.8,-0.6});
    std::vector<double> w(8,0.0); w[1]=0.3;
    std::vector<double> cij(8,0.0), wcij(8,0.0), ckl(8,0.0), xy(16,0.0);
    dijkl1(cv,2,1,w,cij,wcij,ckl,xy);
    std::printf("wcij[0..3]: %.6f %.6f %.6f %.6f\n", wcij[0], wcij[1], wcij[2], wcij[3]);
    std::printf("xy[0,1,5,13,14]: %.6f %.6f %.6f %.6f %.6f\n", xy[0], xy[1], xy[5], xy[13], xy[14]);
    return 0;
}
