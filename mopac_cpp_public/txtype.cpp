// txtype.cpp
#include "txtype.h"
#include "common_arrays_C.h"
using namespace common_arrays_C;
void txtype(int& jj, int* jtype, char letter) {
    int j = 1;
    for (int i = 2; i <= jj; ++i) {
        bool found = false;
        for (int k = 1; k <= j; ++k) if (jtype[k] == jtype[i]) { found = true; break; }
        if (!found) { j++; jtype[j] = jtype[i]; }
    }
    jj = j;
    for (int loop = 1; loop <= 4; ++loop) {
        int m = 0, k = 0;
        for (int i = 1; i <= jj; ++i) if (nat[jtype[i]] != 1) { m++; k = i; }
        if (m == 1) {
            txtatm[jtype[k]][14] = letter;
            int jj2 = jtype[k];
            m = 0;
            for (int i = 1; i <= nbonds[jj2]; ++i) {
            int nbr = ibonds[i][jj2];
                if (txtatm[nbr][15] != ' ') {
                    txtatm[jj2][15] = txtatm[nbr][15];
                    m++;
                }
            }
            if (m != 1) txtatm[jj2][15] = ' ';
            if (txtatm[jj2].substr(17,3) == "TRP" && letter == 'H') txtatm[jj2][15] = '2';
        } else {
            m = 0;
            if (txtatm[j][17] == 'T' && txtatm[j].substr(17,3) == "TRP" && letter == 'Z') m = 1;
            for (int i = 1; i <= jj; ++i) {
                int jj2 = jtype[i];
                if (nat[jj2] != 1) {
                    m++;
                    txtatm[jj2][14] = letter;
                    char threshold = (char)((int)'0' + std::min(9, m));
                    if (txtatm[jj2].substr(17,3) != "UNK" || txtatm[jj2][15] == ' ' || txtatm[jj2][15] > threshold) {
                        txtatm[jj2][15] = threshold;
                    }
                }
            }
        }
    }
}
