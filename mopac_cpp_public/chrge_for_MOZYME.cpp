// chrge_for_MOZYME.cpp
#include "chrge_for_MOZYME.h"
extern int ijbo(int i, int j);
void chrge_for_MOZYME(const double* p, double* q, int numat, const int* iorbs) {
    for (int i=1;i<=numat;++i) {
        int ii=ijbo(i,i);
        double sum=0.0;
        for (int j=1;j<=iorbs[i];++j) { ii+=j; sum+=p[ii]; }
        q[i]=sum;
    }
}
