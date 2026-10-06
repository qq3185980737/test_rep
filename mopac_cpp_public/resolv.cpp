// resolv.cpp — resolve degenerate LMOs.
#include "resolv.h"
#include <cmath>
#include "common_arrays_C.h"
#include "molkst_C.h"
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using molkst_C::norbs;
using molkst_C::numat;
extern void rsp(double*, int, double*, double*);

void resolv(double* c, const double* cold, int mdim, const double* eig, int nocc) {
    double thresh = 1e-3, half = 0.5 - thresh;
    auto C=[&](int k,int l){ return c[(l-1)*mdim+(k-1)]; };
    auto CO=[&](int k,int l){ return cold[(l-1)*mdim+(k-1)]; };
    for (int loop=1; loop<=nocc; ++loop) {
        bool cand=true;
        for (int j=1;j<=numat;++j) {
            if (nlast[j]==nfirst[j]) continue;
            double suml=0;
            for (int k=nfirst[j];k<=nlast[j];++k) suml += C(k,loop)*C(k,loop);
            // F90: only LMOs that are (per atom) <= thresh or >= half are
            // candidates; a "moderate" contribution disqualifies.
            if (suml>thresh && suml<half) { cand=false; break; }
        }
        if (!cand) continue;
        int idegen[5]; int nsec=1; idegen[1]=loop;
        for (int i=loop+1;i<=nocc;++i) {
            bool sim=true;
            for (int j=1;j<=numat;++j) {
                double suml=0,sumi=0;
                for (int k=nfirst[j];k<=nlast[j];++k){
                    suml+=C(k,loop)*C(k,loop);
                    sumi+=C(k,i)*C(k,i);
                }
                if (std::fabs(suml-sumi)>thresh){sim=false;break;}
            }
            if (sim){ ++nsec; idegen[nsec]=i; }
        }
        if (nsec!=1) {
            double sec[55]={0}; double eigs[5]={0}; double vec[17]={0};
            int ij=0;
            for (int ii=1;ii<=nsec;++ii){
                int i=idegen[ii];
                for (int jj=1;jj<=ii;++jj){
                    int j=idegen[jj]; double sum=0;
                    for (int l=1;l<=nocc;++l){
                        double coi=0,coj=0;
                        for (int k=1;k<=norbs;++k){
                            coi+=CO(k,l)*C(k,i);
                            coj+=CO(k,l)*C(k,j);
                        }
                        sum += coi*eig[l-1]*coj;
                    }
                    // F90 packs sec(++ij); C++ rsp reads 0-based sec[0..], so store at sec[ij-1].
                    ++ij; sec[ij-1]=sum;
                }
            }
            rsp(sec,nsec,eigs,vec);
            if (nsec==2){
                int m1=idegen[1],m2=idegen[2];
                for (int i=1;i<=norbs;++i){
                    double s1=vec[0]*C(i,m1)+vec[1]*C(i,m2);
                    double s2=vec[2]*C(i,m1)+vec[3]*C(i,m2);
                    c[(m1-1)*mdim+(i-1)]=s1;
                    c[(m2-1)*mdim+(i-1)]=s2;
                }
            } else if (nsec==3){
                int m1=idegen[1],m2=idegen[2],m3=idegen[3];
                for (int i=1;i<=norbs;++i){
                    double s1=vec[0]*C(i,m1)+vec[1]*C(i,m2)+vec[2]*C(i,m3);
                    double s2=vec[3]*C(i,m1)+vec[4]*C(i,m2)+vec[5]*C(i,m3);
                    double s3=vec[6]*C(i,m1)+vec[7]*C(i,m2)+vec[8]*C(i,m3);
                    c[(m1-1)*mdim+(i-1)]=s1;
                    c[(m2-1)*mdim+(i-1)]=s2;
                    c[(m3-1)*mdim+(i-1)]=s3;
                }
            } else if (nsec==4){
                int m1=idegen[1],m2=idegen[2],m3=idegen[3],m4=idegen[4];
                for (int i=1;i<=norbs;++i){
                    double s1=vec[0]*C(i,m1)+vec[1]*C(i,m2)+vec[2]*C(i,m3)+vec[3]*C(i,m4);
                    double s2=vec[4]*C(i,m1)+vec[5]*C(i,m2)+vec[6]*C(i,m3)+vec[7]*C(i,m4);
                    double s3=vec[8]*C(i,m1)+vec[9]*C(i,m2)+vec[10]*C(i,m3)+vec[11]*C(i,m4);
                    double s4=vec[12]*C(i,m1)+vec[13]*C(i,m2)+vec[14]*C(i,m3)+vec[15]*C(i,m4);
                    c[(m1-1)*mdim+(i-1)]=s1;
                    c[(m2-1)*mdim+(i-1)]=s2;
                    c[(m3-1)*mdim+(i-1)]=s3;
                    c[(m4-1)*mdim+(i-1)]=s4;
                }
            }
        }
    }
}
