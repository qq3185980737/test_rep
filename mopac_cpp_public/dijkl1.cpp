// dijkl1.cpp
#include "dijkl1.h"
#include <vector>
#include "common_arrays_C.h"
#include "formxy.h"
#include "meci_C.h"
#include "molkst_C.h"
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using meci_C::nmos;
using molkst_C::numat;

void dijkl1(const std::vector<std::vector<double>>& cv, int n, int nati,
            const std::vector<double>& w,
            std::vector<double>& cij, std::vector<double>& wcij,
            std::vector<double>& ckl,
            std::vector<double>& xy) {
    (void)cv; (void)w; (void)n;
    static const int nb[9] = {1,0,0,10,0,0,0,0,45}; // 0..8
    // c(ip,i) accessed via cv[ip-1][i-1]
    auto cget=[&](int ip,int i)->double{ return cv[ip-1][i-1]; };
    int na = nmos;
    for (int i=1;i<=na;++i) {
        for (int j=1;j<=i;++j) {
            int ipq=0;
            for (int ii=1;ii<=numat;++ii) {
                if (ii==nati) continue;
                for (int ip=nfirst[ii];ip<=nlast[ii];++ip) {
                    int m=ip-nfirst[ii]+1;
                    if (m>0) {
                        for (int mm=1;mm<=m;++mm) {
                            int idx=nfirst[ii]+mm-1;
                            cij[ipq+mm]=cget(ip,i)*cget(idx,j)+cget(ip,j)*cget(idx,i);
                        }
                        ipq+=m;
                    }
                }
            }
            int i77=ipq+1;
            for (int ip=nfirst[nati];ip<=nlast[nati];++ip) {
                int m=ip-nfirst[nati]+1;
                if (m>0) {
                    for (int mm=1;mm<=m;++mm) {
                        int idx=nfirst[nati]+mm-1;
                        cij[ipq+mm]=cget(ip,i)*cget(idx,j)+cget(ip,j)*cget(idx,i);
                    }
                    ipq+=m;
                }
            }
            // wcij is 1-based semantic; clear physical [0..ipq-1] (F90 wcij(:ipq)=0).
            for (int q=0;q<ipq;++q) wcij[q]=0;
            int kr=1, js=1;
            int nbj=nlast[nati]-nfirst[nati];
            if (nbj>=0) {
                for (int ii=1;ii<=numat;++ii) {
                    if (ii==nati) continue;
                    int nbi=nlast[ii]-nfirst[ii];
                    if (nbi<0) continue;
                    // w is 1-based-padding (element n at physical n): view
                    // must start at w.begin()+kr (NOT +kr-1, which is element kr-1).
                    std::vector<double> wv(w.begin()+kr, w.begin()+kr+nb[nbj]*nb[nbi]+1);
                    std::vector<double> ca(cij.begin()+i77-1, cij.begin()+i77-1+nb[nbj]+1);
                    std::vector<double> cb(cij.begin()+js-1, cij.begin()+js-1+nb[nbi]+1);
                    std::vector<double> wca(nb[nbj]+1, 0.0);
                    std::vector<double> wcb(nb[nbi]+1, 0.0);
                    for (int t=1;t<=nb[nbj];++t) wca[t]=wcij[i77-1+t-1];
                    for (int t=1;t<=nb[nbi];++t) wcb[t]=wcij[js-1+t-1];
                    int kr2 = kr;
                    formxy(wv, kr2, wca, wcb, ca, cb, nb[nbj], nb[nbi]);
                    // formxy outputs are 1-based vectors (element 0 padding).
                    for (int t=1;t<=nb[nbj];++t) wcij[i77-1+(t-1)]=wca[t];
                    for (int t=1;t<=nb[nbi];++t) wcij[js-1+(t-1)]=wcb[t];
                    kr = kr2;
                    js += nb[nbi];
                }
            }
            for (int k=1;k<=i;++k) {
                int ll = (k==i)?j:k;
                for (int l=1;l<=ll;++l) {
                    ipq=0;
                    for (int ii=1;ii<=numat;++ii) {
                        if (ii==nati) continue;
                        for (int ip=nfirst[ii];ip<=nlast[ii];++ip) {
                            int m=ip-nfirst[ii]+1;
                            if (m>0) {
                                for (int mm=1;mm<=m;++mm) {
                                    int idx=nfirst[ii]+mm-1;
                                    ckl[ipq+mm]=cget(ip,k)*cget(idx,l)+cget(ip,l)*cget(idx,k);
                                }
                                ipq+=m;
                            }
                        }
                    }
                    for (int ip=nfirst[nati];ip<=nlast[nati];++ip) {
                        int m=ip-nfirst[nati]+1;
                        if (m>0) {
                            for (int mm=1;mm<=m;++mm) {
                                int idx=nfirst[nati]+mm-1;
                                ckl[ipq+mm]=cget(ip,k)*cget(idx,l)+cget(ip,l)*cget(idx,k);
                            }
                            ipq+=m;
                        }
                    }
                    double sum=0;
                    // ckl built 1-based physical (ckl[1]=first); wcij 1-based semantic (physical [0] = slot 1).
                    for (int q=1;q<=ipq;++q) sum+=ckl[q]*wcij[q-1];
                    // xy(nmos^4) col-major
                    auto xyidx=[&](int a,int b,int c,int d){
                        return ((d-1)*nmos*nmos*nmos+(c-1)*nmos*nmos+(b-1)*nmos+(a-1));
                    };
                    xy[xyidx(i,j,k,l)]=sum; xy[xyidx(i,j,l,k)]=sum;
                    xy[xyidx(j,i,k,l)]=sum; xy[xyidx(j,i,l,k)]=sum;
                    xy[xyidx(k,l,i,j)]=sum; xy[xyidx(k,l,j,i)]=sum;
                    xy[xyidx(l,k,i,j)]=sum; xy[xyidx(l,k,j,i)]=sum;
                }
            }
        }
    }
}
