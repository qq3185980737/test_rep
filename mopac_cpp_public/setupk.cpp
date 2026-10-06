// setupk.cpp
#include "setupk.h"
#include <vector>
#include <algorithm>
#include "molkst_C.h"
using molkst_C::numat;
namespace mozyme_C {
    std::vector<int> icocc, ncf, nncf, kopt;
}
void setupk(int nocc1) {
    using namespace mozyme_C;
    std::fill(kopt.begin(),kopt.end(),0);
    for (int i=1;i<=nocc1;++i) {
        int j=nncf[i];
        for (int k=1;k<=ncf[i];++k)
            kopt[icocc[k+j]]=1;
    }
    int l=0;
    for (int i=1;i<=numat;++i) {
        if (kopt[i]==1) { ++l; kopt[l]=i; }
    }
    if (l!=numat) kopt[l+1]=0;
}
