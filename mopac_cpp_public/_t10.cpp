#include <cstdio>
#include <vector>
#include "formxy.h"
int main(){
    // replicate partxy atom2 formxy call (kr=2 after atom1 1-center)
    std::vector<double> w(16,0.0); w[1]=0.4; w[2]=0.3; w[3]=0.2;
    double c34_[8]={0,0.5,1.5,0,0,0,0,0};
    double pq34_[8]={0,0,0,0,0,0,0,0};
    int kr=2;
    std::vector<double> wv(w.begin()+kr-1, w.begin()+kr-1+1*1+1);  // w[1..2]
    std::vector<double> ca(c34_+2-1, c34_+2-1+2);   // c34[1..2]
    std::vector<double> cb(c34_+1-1, c34_+1-1+2);   // c34[0..1]
    std::vector<double> wca(pq34_+2-1, pq34_+2-1+2); // pq34[1..2]
    std::vector<double> wcb(pq34_+1-1, pq34_+1-1+2); // pq34[0..1]
    std::printf("wv[0..2]: %f %f %f\n", wv[0], wv[1], wv[2]);
    std::printf("ca[0..1]: %f %f\n", ca[0], ca[1]);
    std::printf("cb[0..1]: %f %f\n", cb[0], cb[1]);
    int kr2=kr;
    formxy(wv, kr2, wca, wcb, ca, cb, 1, 1);
    std::printf("kr %d -> %d\n", kr, kr2);
    std::printf("wca: %f %f\n", wca[0], wca[1]);
    std::printf("wcb: %f %f\n", wcb[0], wcb[1]);
    return 0;
}
