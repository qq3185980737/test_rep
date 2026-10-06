// mecih.cpp — build the packed canonical CI matrix CIMAT.
#include "mecih.h"

#include <cmath>
#include <vector>

#include "aabbcd.h"
#include "aababc.h"
#include "aabacd.h"
#include "babbcd.h"
#include "babbbc.h"
#include "meci_C.h"

using namespace meci_C;

void mecih(double* diag, double* cimat, int nmos, const int* /*lab_in*/,
           double* xy) {
    int ik = 0;
    for (iiloop = 1; iiloop <= lab; ++iiloop) {
        is = 2;
        for (jloop = 1; jloop <= iiloop; ++jloop) {
            ++ik;
            cimat[ik - 1] = 0.0;   // physical 0-based packed slot (rsp reads a[0]=slot1)
            int ix = 0, iy = 0;
            for (int j = 1; j <= nmos; ++j) {
                ix += std::abs(microa[j][iiloop] - microa[j][jloop]);
                iy += std::abs(microb[j][iiloop] - microb[j][jloop]);
            }
            if (ix + iy > 4 || nalmat[iiloop] != nalmat[jloop]) {
                if (iiloop == 2 && jloop == 1)
                    fprintf(stderr, "[MECIH] (2,1) SKIP ix=%d iy=%d nalma=%d nalmb=%d\n",
                            ix, iy, nalmat[iiloop], nalmat[jloop]); fflush(stderr);
                continue;
            }
            if (iiloop == 2 && jloop == 1)
                fprintf(stderr, "[MECIH] (2,1) PASS ix=%d iy=%d nalma=%d nalmb=%d\n",
                        ix, iy, nalmat[iiloop], nalmat[jloop]); fflush(stderr);

            std::vector<int> ma1(nmos + 1), mb1(nmos + 1), ma2(nmos + 1),
                mb2(nmos + 1);
            for (int j = 1; j <= nmos; ++j) {
                ma1[j] = microa[j][iiloop];
                mb1[j] = microb[j][iiloop];
                ma2[j] = microa[j][jloop];
                mb2[j] = microb[j][jloop];
            }

            if (ix + iy == 4) {
                if (ix == 0)
                    cimat[ik - 1] = babbcd(&ma1[0], &mb1[0], &ma2[0], &mb2[0],
                                           nmos, xy);
                else if (ix == 2)
                    cimat[ik - 1] = aabbcd(&ma1[0], &mb1[0], &ma2[0], &mb2[0],
                                           nmos, xy);
                else
                    cimat[ik - 1] = aabacd(&ma1[0], &mb1[0], &ma2[0], &mb2[0],
                                           nmos, xy);
            } else if (ix == 2) {
                cimat[ik - 1] = aababc(&ma1[0], &mb1[0], &ma2[0], nmos, xy);
                if (iiloop == 2 && jloop == 1)
                    fprintf(stderr, "[MECIH] (2,1) ix=%d iy=%d nalma=%d nalmb=%d val=%.8f\n",
                            ix, iy, nalmat[iiloop], nalmat[jloop], cimat[ik - 1]); fflush(stderr);
            } else if (iy == 2) {
                cimat[ik - 1] = babbbc(&ma1[0], &mb1[0], &mb2[0], nmos, xy);
                if (iiloop == 2 && jloop == 1)
                    fprintf(stderr, "[MECIH] (2,1) babbbc val=%.8f\n", cimat[ik - 1]); fflush(stderr);
            } else {
                cimat[ik - 1] = diag[iiloop];
            }
        }
        ispqr[iiloop][1] = is - 1;
    }
}
