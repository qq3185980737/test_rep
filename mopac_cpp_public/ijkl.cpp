// ijkl.cpp
#include "ijkl.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "partxy.h"
namespace cosmo_C { extern bool iseps; }
void ciint(double*, double*) {}

void ijkl(double* cp, double* cf, int nelec, int nmos,
          double* dijkl, double* cij, double* ckl, double* wcij, double* xy) {
    using namespace common_arrays_C;
    int norbs = molkst_C::norbs;
    int numat = molkst_C::numat;
    int ij = 0;
    for (int i=1;i<=nmos;++i) {
        for (int j=1;j<=i;++j) {
            ++ij;
            int ipq = 0;
            for (int ii=1;ii<=numat;++ii) {
                for (int ip=nfirst[ii];ip<=nlast[ii];++ip) {
                    int n = ip - nfirst[ii] + 1;
                    if (n > 0) {
                        for (int m=1;m<=n;++m) {
                            int idx = nfirst[ii]+m-1;
                            cij[ipq+m] = cp[(i-1)*norbs+ip-1]*cp[(j-1)*norbs+idx-1]
                                       + cp[(j-1)*norbs+ip-1]*cp[(i-1)*norbs+idx-1];
                        }
                        ipq = n + ipq;
                    }
                }
            }
            partxy(cij, wcij, w.data());
            for (int k=1;k<=norbs;++k) {
                for (int l=1;l<=nmos;++l) {
                    ipq = 0;
                    for (int ii=1;ii<=numat;++ii) {
                        for (int ip=nfirst[ii];ip<=nlast[ii];++ip) {
                            int n = ip - nfirst[ii] + 1;
                            if (n > 0) {
                                for (int m=1;m<=n;++m) {
                                    int idx = nfirst[ii]+m-1;
                                    ckl[ipq+m] = cf[(k-1)*norbs+ip-1]*cp[(l-1)*norbs+idx-1]
                                               + cp[(l-1)*norbs+ip-1]*cf[(k-1)*norbs+idx-1];
                                }
                                ipq = n + ipq;
                            }
                        }
                    }
                    double sum=0;
                    // wcij is 1-based semantic (partxy wrote physical [ls-1]);
                    // ckl built 1-based physical: ckl[ii]=ckl(ii).
                    for (int ii=1;ii<=ipq;++ii) sum += ckl[ii]*wcij[ii-1];
                    if (ij==2 && k==2 && l==2) {
                        fprintf(stderr, "[IJKL] ipq=%d sum=%.6f ckl:", ipq, sum);
                        for (int ii=1;ii<=ipq;++ii) fprintf(stderr, " %+.4f", ckl[ii]);
                        fprintf(stderr, "\n[IJKL] wcij:"); fflush(stderr);
                        for (int ii=1;ii<=ipq;++ii) fprintf(stderr, " %+.4f", wcij[ii-1]);
                        fprintf(stderr, "\n"); fflush(stderr);
                    }
                    // Fortran dijkl(norbs,nmos,mpair): linear = k + (l-1)*norbs + (ij-1)*norbs*nmos
                    dijkl[(ij-1)*norbs*nmos + (l-1)*norbs + k-1] = sum;
                }
            }
        }
    }
    int mpair = nmos*(nmos+1)/2;
    for (int k=1;k<=nmos;++k) {
        int kk = nelec + k;
        for (int l=1;l<=nmos;++l) {
            int ij2 = 0;
            for (int i=1;i<=nmos;++i) {
                for (int j=1;j<=i;++j) {
                    ++ij2;
                    double s = dijkl[(ij2-1)*norbs*nmos + (l-1)*norbs + kk-1];
                    // xy(nmos^4) col-major: xy(i,j,k,l)
                    auto xyidx = [&](int a,int b,int cc,int d){
                        return ((d-1)*nmos*nmos*nmos + (cc-1)*nmos*nmos + (b-1)*nmos + (a-1));
                    };
                    xy[xyidx(i,j,k,l)]=s; xy[xyidx(i,j,l,k)]=s;
                    xy[xyidx(j,i,k,l)]=s; xy[xyidx(j,i,l,k)]=s;
                    xy[xyidx(k,l,i,j)]=s; xy[xyidx(k,l,j,i)]=s;
                    xy[xyidx(l,k,i,j)]=s;
                }
            }
        }
    }
    if (cosmo_C::iseps) {
        for (int i=1;i<=nmos;++i) {
            int ipq=0;
            for (int ii=1;ii<=numat;++ii)
                for (int ip=nfirst[ii];ip<=nlast[ii];++ip) {
                    int n=ip-nfirst[ii]+1;
                    if (n>0) {
                        for (int m=1;m<=n;++m) {
                            int idx=nfirst[ii]+m-1;
                            cij[ipq+m]=cp[(i-1)*norbs+ip-1]*cp[(i-1)*norbs+idx-1];
                        }
                        ipq+=n;
                    }
                }
            ciint(cij,wcij);
            for (int l=1;l<=i;++l) {
                ipq=0;
                for (int ii=1;ii<=numat;++ii)
                    for (int ip=nfirst[ii];ip<=nlast[ii];++ip) {
                        int n=ip-nfirst[ii]+1;
                        if (n>0) {
                            for (int m=1;m<=n;++m) {
                                int idx=nfirst[ii]+m-1;
                                ckl[ipq+m]=cp[(l-1)*norbs+ip-1]*cp[(l-1)*norbs+idx-1];
                            }
                            ipq+=n;
                        }
                    }
            }
        }
    }
}
