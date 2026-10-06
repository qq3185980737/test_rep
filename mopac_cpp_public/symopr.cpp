// symopr.cpp
#include "symopr.h"
void symopr(int numat, double* coord, int jump, double* r) {
    for (int i=1;i<=numat;++i) {
        double help[4] = {0.0, coord[(i-1)*3+0], coord[(i-1)*3+1], coord[(i-1)*3+2]};
        for (int j=1;j<=3;++j) {
            double s=0.0;
            if (jump>=0) {
                for (int k=1;k<=3;++k) s += r[(j-1)*3+(k-1)]*help[k];
            } else {
                for (int k=1;k<=3;++k) s += r[(k-1)*3+(j-1)]*help[k];
            }
            coord[(i-1)*3+(j-1)] = s;
        }
    }
}
