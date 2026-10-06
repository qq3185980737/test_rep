// molmec_C.h — C++ translation of MOPAC 2016 "molmec_C.F90".
#pragma once

namespace molmec_C {
// nhco(4,4000): hydrogen-bond / interaction list; second index 1-based -> [4001]
extern int nhco[4][4001];
extern int nnhco;
extern double htype;
}  // namespace molmec_C
