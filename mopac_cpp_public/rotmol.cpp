// rotmol.cpp
#include "rotmol.h"
extern void symopr(int numat, double* coord, int, double* r);
void rotmol(int numat, double* coord, double sina, double cosa, int i,
            int j, double* r) {
    symopr(numat, coord, -1, r);
    for (int k=1;k<=3;++k) {
        double buff = -sina*r[(i-1)*3+(k-1)] + cosa*r[(j-1)*3+(k-1)];
        r[(i-1)*3+(k-1)] = cosa*r[(i-1)*3+(k-1)] + sina*r[(j-1)*3+(k-1)];
        r[(j-1)*3+(k-1)] = buff;
    }
    symopr(numat, coord, 1, r);
}
