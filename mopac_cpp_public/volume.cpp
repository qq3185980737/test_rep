// volume.cpp
#include "volume.h"
#include <cmath>
double volume(double* v, int ndim) {
    double a=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (ndim==1) return a;
    double b=std::sqrt(v[3]*v[3]+v[4]*v[4]+v[5]*v[5]);
    double dx=v[0]-v[3],dy=v[1]-v[4],dz=v[2]-v[5];
    double g=std::sqrt(dx*dx+dy*dy+dz*dz);
    double c=(a*a+b*b-g*g)/(2.0*a*b);
    double sing=std::sqrt(1.0-c*c);
    if (ndim==2) return a*b*sing;
    return std::abs((v[1]*v[5]-v[2]*v[4])*v[6]+(v[2]*v[3]-v[0]*v[5])*v[7]+(v[0]*v[4]-v[1]*v[3])*v[8]);
}
