// set.cpp
#include "set.h"
#include <cmath>
namespace overlaps_C {
    thread_local int isp=0, ips=0;
    thread_local double sa=0, sb=0;
    thread_local double a[64]={}, b[64]={};
    // fact[18] is defined in overlaps_C.cpp (module overlaps_C).
}
static void aintgs(double x, int k) {
    using namespace overlaps_C;
    double c=std::exp(-x);
    a[1]=c/x;
    for (int i=1;i<=k;++i) a[i+1]=(a[i]*i+c)/x;
}
static void bintgs(double x, int k) {
    using namespace overlaps_C;
    int io=0;
    double absx=std::fabs(x);
    int last;
    if (absx>3.0) goto l40;
    if (absx>2.0) { if (k<=10) goto l40; last=15; goto l60; }
    if (absx>1.0) { if (k<=7) goto l40; last=12; goto l60; }
    if (absx>0.5) { if (k<=5) goto l40; last=7; goto l60; }
    if (absx<=1e-6) goto l90;
    last=6;
  l60:
    for (int i=io;i<=k;++i) {
        double y=0.0;
        for (int m=io;m<=last;++m) {
            double xf=(m!=0)?fact[m]:1.0;
            y=y+std::pow(-x,m)*(2*((m+i+1)%2))/(xf*(m+i+1));
        }
        b[i+1]=y;
    }
    return;
  l40:
    { double expx=std::exp(x), expmx=1.0/expx;
      b[1]=(expx-expmx)/x;
      for (int i=1;i<=k;++i) {
          double sgn=(i%2==0)?1.0:-1.0;
          b[i+1]=(i*b[i]+sgn*expx-expmx)/x;
      }
    }
    return;
  l90:
    for (int i=io;i<=k;++i) b[i+1]=(2*((i+1)%2))/(i+1.0);
}
void set_fn(double s1, double s2, int na, int nb, double rab, int ii) {
    using namespace overlaps_C;
    if (na<=nb) { isp=1; ips=2; sa=s1; sb=s2; }
    else { isp=2; ips=1; sa=s2; sb=s1; }
    int j=ii+2;
    if (ii>3) j-=1;
    double alpha=0.5*rab*(sa+sb);
    double beta=0.5*rab*(sb-sa);
    aintgs(alpha,j-1);
    bintgs(beta,j-1);
}
