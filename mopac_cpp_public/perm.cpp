// perm.cpp — CI microstate generator (combinations of nels electrons in nmos MOs)
#include "perm.h"
#include <vector>
#include <algorithm>
#include <functional>
void perm(int* iperm, int nels, int nmos, int& nperms, int limci) {
    // iperm is (nmos, maxci*4) col-major; write into flat 1-based rows.
    int Upper = nmos - nels;
    int Lower = 1;
    int Prm = 1;
    int El = nels;
    // Occ(-1..nels); map Occ[k] -> idx k+1 (size nels+2)
    std::vector<int> Occ(nels+2, 0);
    for (int i=-1; i<=nels-1; ++i) Occ[i+1] = nmos + 1 - i;
    auto iat = [&](int i,int p){ return (p-1)*nmos + (i-1); }; // col-major (nmos,nperms)
    // recursive lambda
    std::function<void(int,int,int)> rperm;
    rperm = [&](int Lower_, int Upper_, int El_) {
        if (El_ != 0) {
            for (int j=Lower_; j<=Upper_; ++j) {
                Occ[El_+1] = j;
                rperm(Occ[El_+1]+1, Occ[El_-2+1]-2, El_-1);
            }
            Occ[El_+1] = Upper_ + 1;
            if (Upper_+1 > nmos || Prm > nperms) return;
        }
        for (int i=1;i<=nmos;++i) iperm[iat(i,Prm)] = 0;
        for (int i=1;i<=nels;++i) iperm[iat(Occ[i+1],Prm)] = 1;
        El = nels + 1;
        if (limci != 0 && Prm > 1) {
            int k = 0;
            for (int j=1;j<=nmos;++j) k += std::abs(iperm[iat(j,Prm)]-iperm[iat(j,1)]);
            if (k > limci) --Prm;
        }
        ++Prm;
    };
    rperm(Lower, Upper, El);
    nperms = Prm - 1;
}
