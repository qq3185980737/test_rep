// prttim.cpp
#include "prttim.h"
void prttim(double tleft, double& tprt, char& txt) {
    tprt=tleft; txt='S';
    if (tprt>=604800.0) { tprt/=604800.0; txt='W'; }
    else if (tprt>=86400.0) { tprt/=86400.0; txt='D'; }
    else if (tprt>=3600.0) { tprt/=3600.0; txt='H'; }
    else if (tprt>=60.0) { tprt/=60.0; txt='M'; }
}
