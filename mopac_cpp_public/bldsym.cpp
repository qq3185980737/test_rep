// bldsym.cpp — C++ translation of MOPAC 2016 "bldsym.F90".

#include "bldsym.h"

#include <cmath>

#include "symmetry_C.h"
#include "mult33.h"

namespace {
// j(1..3, 1..20): diagonal axis signs / order and improper flag.
const int j[21][4] = {
    /*0*/ {0},
    /*1*/ {0, 1, -1, -1},
    /*2*/ {0, -1, 1, -1},
    /*3*/ {0, -1, -1, 1},
    /*4*/ {0, 1, 1, -1},
    /*5*/ {0, 1, -1, 1},
    /*6*/ {0, -1, 1, 1},
    /*7*/ {0, -1, -1, -1},
    /*8*/ {0, 3, 0, 1},
    /*9*/ {0, 4, 0, 1},
    /*10*/ {0, 5, 0, 1},
    /*11*/ {0, 6, 0, 1},
    /*12*/ {0, 7, 0, 1},
    /*13*/ {0, 8, 0, 1},
    /*14*/ {0, 4, 0, -1},
    /*15*/ {0, 6, 0, -1},
    /*16*/ {0, 8, 0, -1},
    /*17*/ {0, 10, 0, -1},
    /*18*/ {0, 12, 0, -1},
    /*19*/ {0, 5, 0, -1},
    /*20*/ {0, 0, 0, -1},
};
const double twopi = 6.2831853071796;
}

void bldsym(int ioper, int iplace) {
    using namespace symmetry_C;
    for (int i = 1; i <= 3; ++i) {
        for (int r = 1; r <= 3; ++r) elem[i][r][iplace] = 0.0;
        elem[i][i][iplace] = (double)j[ioper][i];
    }
    if (ioper != 20) {
        if (j[ioper][1] >= 2) {
            double angle = twopi / (double)j[ioper][1];
            elem[1][1][iplace] = std::cos(angle);
            elem[2][2][iplace] = elem[1][1][iplace];
            elem[2][1][iplace] = std::sin(angle);
            elem[1][2][iplace] = -elem[2][1][iplace];
        }
        if (ioper == 8 || ioper == 15) {
            // cub is stored with direct indexing (cub[1..3][1..3]); pack into a
            // compact column-major 3x3 for mult33 (Fortran layout).
            double cub3[9] = {};
            for (int a = 1; a <= 3; ++a)
                for (int b = 1; b <= 3; ++b) cub3[(b - 1) * 3 + (a - 1)] = cub[a][b];
            mult33(cub3, iplace);
        }
        return;
    }
    elem[1][2][iplace] = 1.0;
    elem[2][1][iplace] = 1.0;
}
