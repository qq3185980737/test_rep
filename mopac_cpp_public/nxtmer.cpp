// nxtmer.cpp
#include "nxtmer.h"
#include "common_arrays_C.h"
using common_arrays_C::nat;
using common_arrays_C::ibonds;
using common_arrays_C::nbonds;
void nxtmer(int iatom, int* nbackb) {
    int ihcr=0, ico=0, ico1=0;
    int jatom=nbackb[3];
    for (int i=1;i<=nbonds[iatom];++i)
        if (nat[ibonds[i][iatom]]==6) ihcr=ibonds[i][iatom];
    int l=0,j=0,jj=0,jjofco=0,jofco=0,iii=0;
    for (int i=1;i<=nbonds[iatom];++i) {
        if (nat[ibonds[i][iatom]]==6) {
            l=ibonds[i][iatom]; jj=0; jjofco=0;
            for (int k=1;k<=nbonds[l];++k) {
                if (nat[ibonds[k][l]]==6) {
                    j=ibonds[k][l];
                    int go1000=0;
                    for (int ii=1;ii<=nbonds[j];++ii)
                        if (nat[ibonds[ii][j]]==8 && nbonds[ibonds[ii][j]]==1) { go1000=ii; break; }
                    if (!go1000) continue;
                    int ii=go1000;
                    jofco=ibonds[ii][j];
                    ico1=ico; ico=j; iii=ii; ihcr=l;
                    int hasN=0;
                    for (ii=1;ii<=nbonds[j];++ii)
                        if (nat[ibonds[ii][j]]==7) { hasN=ii; break; }
                    if (hasN) {
                        jatom=ibonds[hasN][j];
                        for (ii=1;ii<=nbonds[jatom];++ii)
                            if (nat[ibonds[ii][jatom]]==6) goto done;
                    } else {
                        int kk=0;
                        for (ii=1;ii<=nbonds[j];++ii)
                            if (nat[ibonds[ii][j]]==8) ++kk;
                        if (kk==2) { jatom=0; goto done; }
                        kk=0;
                        for (ii=1;ii<=nbonds[j];++ii)
                            if (nat[ibonds[ii][j]]==6) ++kk;
                        if (kk==1) { jj=j; jjofco=ibonds[iii][j]; }
                    }
                }
            }
            if (jj) j=jj;
            if (jjofco) jofco=jjofco;
        }
    }
done:
    if (ico1) {
        for (int i=1;i<=nbonds[ico1];++i)
            if (ibonds[i][ico1]==jofco) ico=ico1;
    }
    nbackb[0]=ihcr; nbackb[1]=ico; nbackb[2]=jofco; nbackb[3]=jatom;
}
