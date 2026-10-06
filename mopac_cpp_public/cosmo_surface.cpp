// cosmo_surface.cpp — standalone surface-geometry helpers translated from
// MOPAC 2016 "cosmo.F90": dvfill (direction grids), ansude (analytic area of
// two intersecting spheres), surclo (closure of concave cavity regions).
// Kept separate from cosmo.cpp so consumers that only need the geometry
// helpers do not pull in the full COSMO solvation driver.
#include "cosmo.h"
#include "molkst_C.h"
#include "cosmo_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "mopend.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace molkst_C;
using namespace cosmo_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;

void ansude(double ra, double rb, double d, double rs,
            double& aar, double& abr, double& ara, double& arad,
            double& arb, double& arbd, double& rinc) {
    double qa = ra + rs, qb = rb + rs;
    double ca = (qa * qa + d * d - qb * qb) / (2.0 * qa * d);
    double cbv = (qb * qb + d * d - qa * qa) / (2.0 * qb * d);
    double sa = std::sqrt(1.0 - ca * ca);
    double sb = std::sqrt(1.0 - cbv * cbv);
    double ta = pi * sa, tb = pi * sb;
    double fza = (1.0 - std::cos(ta)) / 2.0;
    double fzb = (1.0 - std::cos(tb)) / 2.0;
    if (sa < 0 || sb < 0) fza = 1.0;
    if (sa < 0 || sb < 0) fzb = 1.0;
    double xa = fzb * rs * (ca + cbv);
    double xb = fza * rs * (ca + cbv);
    double ya = ra * sa - fzb * rb * sb;
    double yb = rb * sb - fza * ra * sa;
    double za = std::sqrt(xa * xa + ya * ya);
    double zb = std::sqrt(xb * xb + yb * yb);
    rinc = 0.5 * (za + zb) / std::sqrt(rs * rs * (ca + cbv) * (ca + cbv) +
                                        (ra * sa - rb * sb) * (ra * sa - rb * sb));
    ara = pi * ra * (2.0 * (1.0 + ca) * ra + sa * za);
    arb = pi * rb * (2.0 * (1.0 + cbv) * rb + sb * zb);
    aar = pi * ra * (sa * za);
    abr = pi * rb * (sb * zb);
    double cad = (qb * qb + d * d - qa * qa) / (2.0 * qa * d * d);
    double cbd = (qa * qa + d * d - qb * qb) / (2.0 * qb * d * d);
    double sad = -ca * cad / sa, sbd = -cbv * cbd / sb;
    double tad = pi * sad, tbd = pi * sbd;
    double fzad = std::sin(ta) * 0.5, fzbd = std::sin(tb) * 0.5;
    if (sa < 0 || sb < 0) fzad = 0.0;
    if (sa < 0 || sb < 0) fzbd = 0.0;
    double xad = rs * ((ca + cbv) * fzbd * tbd + fzb * (cad + cbd));
    double xbd = rs * ((ca + cbv) * fzad * tad + fza * (cad + cbd));
    double yad = ra * sad - fzbd * tbd * rb * sb - fzb * rb * sbd;
    double ybd = rb * sbd - fzad * tad * ra * sa - fza * ra * sad;
    double zad = (xa * xad + ya * yad) / za;
    double zbd = (xb * xbd + yb * ybd) / zb;
    arad = pi * ra * (sad * za + sa * zad + 2.0 * ra * cad);
    arbd = pi * rb * (sbd * zb + sb * zbd + 2.0 * rb * cbd);
}
void dvfill(int nppa, double* dirvec) {
    auto D=[&](int ix,int i) -> double& { return dirvec[(i-1)*4 + (ix-1)]; };
    auto Ds=[&](int ix,int i){ return &dirvec[(i-1)*4 + (ix-1)]; };
    int kset[31][3], fset[21][4];
    {
        int ks[60]={1,2,1,3,1,4,1,5,1,6,12,11,12,10,12,9,12,8,12,7,
                    2,3,3,4,4,5,5,6,6,2,7,8,8,9,9,10,10,11,11,7,2,
                    7,7,3,3,8,8,4,4,9,9,5,5,10,10,6,6,11,11,2};
        for(int i=0;i<60;++i){ int kk=i/2+1,kj=i%2; kset[kk][kj+1]=ks[i]; }
        int fs[60]={1,2,3,1,3,4,1,4,5,1,5,6,1,6,2,12,11,10,12,10,9,12,9,8,
                    12,8,7,12,7,11,2,3,7,3,4,8,4,5,9,5,6,10,6,2,11,7,8,3,
                    8,9,4,9,10,5,10,11,6,11,7,2};
        for(int i=0;i<60;++i){ int kk=i/3+1,kj=i%3; fset[kk][kj+1]=fs[i]; }
    }
    D(1,1)=-1; D(2,1)=0; D(3,1)=0;
    int nd=1;
    double r=std::sqrt(0.8), h=std::sqrt(0.2);
    for (int i=-1;i<=1;i+=2)
      for (int j=1;j<=5;++j) {
        ++nd;
        double beta = 1.0 + j*0.4*pi + (i+1)*0.1*pi;
        D(2,nd)=r*std::cos(beta); D(3,nd)=r*std::sin(beta); D(1,nd)=i*h;
      }
    D(1,12)=1; D(2,12)=0; D(3,12)=0;
    nd = 12;
    double cphi=std::cos(1.0), sphi=std::sin(1.0);
    for (int i=1;i<=12;++i) {
        double xx=D(1,i), yy=D(2,i);
        D(1,i)=cphi*xx+sphi*yy; D(2,i)=-sphi*xx+cphi*yy;
    }
    int m2=(nppa-2)/10;
    int m=(int)std::floor(std::sqrt((double)m2)+0.5);
    int k=0;
    if (m2 != m*m) { k=1; m2/=3; m=(int)std::floor(std::sqrt((double)m2)+0.5); }
    if (10*(k?3:1)*m*m+2 != nppa) return;
    for (int i=1;i<=30;++i) {
        int na=kset[i][1], nb=kset[i][2];
        for (int j=1;j<=m-1;++j) {
            ++nd;
            for (int ix=1;ix<=3;++ix)
                D(ix,nd)=D(ix,na)*(m-j)+D(ix,nb)*j;
        }
    }
    for (int i=1;i<=20;++i) {
        int na=fset[i][1], nb=fset[i][2], nc=fset[i][3];
        for (int j1=1;j1<=m-1;++j1)
          for (int j2=1;j2<=m-j1-1;++j2) {
            ++nd;
            for (int ix=1;ix<=3;++ix)
                D(ix,nd)=D(ix,na)*(m-j1-j2)+D(ix,nb)*j1+D(ix,nc)*j2;
          }
        if (k != 0) {
            double t=1.0/3;
            for (int j1=0;j1<=m-1;++j1)
              for (int j2=0;j2<=m-j1-1;++j2) {
                ++nd;
                for (int ix=1;ix<=3;++ix)
                    D(ix,nd)=D(ix,na)*(m-j1-j2-2*t)+D(ix,nb)*(j1+t)+D(ix,nc)*(j2+t);
              }
            t=2.0/3;
            for (int j1=0;j1<=m-2;++j1)
              for (int j2=0;j2<=m-j1-2;++j2) {
                ++nd;
                for (int ix=1;ix<=3;++ix)
                    D(ix,nd)=D(ix,na)*(m-j1-j2-2*t)+D(ix,nb)*(j1+t)+D(ix,nc)*(j2+t);
              }
        }
    }
    double sumar=0;
    for (int i=1;i<=nppa;++i) {
        double dist=0;
        for (int ix=1;ix<=3;++ix) dist += D(ix,i)*D(ix,i);
        dist = 1.0/std::sqrt(dist);
        double dist2 = (m*dist)*(m*dist);
        for (int ix=1;ix<=3;++ix) D(ix,i)*=dist;
        double ar = (i<=12) ? 5.0 : 6.0*dist2;
        D(4,i)=ar; sumar+=ar;
    }
    sumar = 4*pi/sumar;
    for (int i=1;i<=nppa;++i) D(4,i)*=sumar;
}

void addnuc() {
}

void surclo(double* coord, const int* nipa, const int* lipa, const bool* din,
            int dim_din, double** rsc, int* isort, int* ipsrs, int* nipsrs,
            const int* nat, const double* srad, int maxrs) {
    (void)dim_din;
    auto dot=[&](const double*a,const double*b){double s=0;for(int ix=0;ix<3;++ix)s+=a[ix]*b[ix];return s;};
    nipc = 0;
    int nrs = 0;
    for (int i = 1; i <= numat; ++i) nipsrs[i] = 0;
    int ilipa = 0;
    for (int ia = 1; ia <= numat - 1; ++ia) {
        if (din[ia]) { ilipa += nipa[ia]; continue; }
        double ra = srad[ia] + rsolv;
        std::array<double,3> xta1{coord[0*0+0],0,0};
        // xta(:,1) = coord(:,ia)
        std::array<double,3> xta;
        for (int ix = 0; ix < 3; ++ix) xta[ix] = coord[(ix)*0+ia]; // placeholder
        // use coord[ix][ia] via flat layout assumed coord[2][numat+1]
        for (int ix = 0; ix < 3; ++ix) xta[ix] = coord[ix*(numat+1)+ia];
    for (int iib = ilipa + 1; iib <= ilipa + nipa[ia]; ++iib) {
            int ib = lipa[iib];
            if (ib <= ia) continue;
            if (din[ib]) continue;
            double rb = srad[ib] + rsolv;
            double dab = 0.0;
            int nsab = 0; (void)nsab;
            std::array<double,3> xx;
            for (int ix = 0; ix < 3; ++ix) {
                double xb = coord[ix*(numat+1)+ib];
                double xa = coord[ix*(numat+1)+ia];
                xx[ix] = xb - xa;
                dab += xx[ix]*xx[ix];
            }
            dab = std::sqrt(dab);
            double aa,ab,aar,abr,aad,abd,rinc;
            ansude(ra-rsolv, rb-rsolv, dab, rsolv, aa,ab,aar,abr,aad,abd,rinc);
            double cosa = (ra*ra+dab*dab-rb*rb)/(2*dab*ra);
            double cosb  = (rb*rb+dab*dab-ra*ra)/(2*dab*rb);
            double sina = std::sqrt(std::max(0.0,1.0-cosa*cosa));
            double da = ra*cosa;
            double hh = ra*sina;
            double ddd = rsolv*(cosa+cosb)/dab;
            double fz1 = (1.0-std::cos(hh*pi/ra))/2.0;
            double fz2 = (1.0-std::cos(hh*pi/rb))/2.0;
            if (cosa*cosb < 0) fz1 = 1.0;
            if (cosa*cosb < 0) fz2 = 1.0;
            double yx1 = rsolv/ra, yx2 = rsolv/rb;
            std::array<double,3> xd;
            for (int ix=0;ix<3;++ix) xd[ix]=coord[ix*(numat+1)+ia]+da*xx[ix]/dab;
            std::array<double,3> rvx;
            rvx[0]=xx[1]*3.0-xx[2]*2.0;
            rvx[1]=xx[2]*1.0-xx[0]*3.0;
            rvx[2]=xx[0]*2.0-xx[1]*1.0;
            double dist=std::sqrt(dot(rvx.data(),rvx.data()));
            for (int ix=0;ix<3;++ix) rvx[ix]=hh*rvx[ix]/dist;
            std::array<double,3> rvy;
            rvy[0]=(xx[1]*rvx[2]-xx[2]*rvx[1])/dab;
            rvy[1]=(xx[2]*rvx[0]-xx[0]*rvx[2])/dab;
            rvy[2]=(xx[0]*rvx[1]-xx[1]*rvx[0])/dab;

            std::array<int,50> iset{};
            std::array<double,50> phiset{}, tarset{};
            int ntrp=0;
            for (int iic = ilipa+1; iic <= ilipa+nipa[ia]; ++iic) {
                int ic = lipa[iic];
                if (ic == ib) continue;
                double rc = srad[ic]+rsolv;
                double dabc=0, sp=0;
                std::array<double,3> xic;
                for (int ix=0;ix<3;++ix) {
                    double xxx = coord[ix*(numat+1)+ic]-xd[ix];
                    xic[ix]=coord[ix*(numat+1)+ic];
                    sp += xxx*xx[ix];
                    dabc += xxx*xxx;
                }
                dabc = std::sqrt(dabc);
                cosa = sp/dab/dabc;
                sina = std::sqrt(std::max(1e-28,1.0-cosa*cosa));
                double cj = (dabc*dabc+hh*hh-rc*rc)/(2*dabc*hh*sina);
                if (cj < 1.0) {} // ntrp2 unused downstream
                if (cj <= 1.0 && cj >= -1.0) {
                    double sj = std::sqrt(std::max(0.0,1.0-cj*cj));
                    std::array<double,3> tvx;
                    for (int ix=0;ix<3;++ix)
                        tvx[ix]=(xic[ix]-xd[ix])-cosa*dabc*xx[ix]/dab;
                    dist=std::sqrt(dot(tvx.data(),tvx.data()));
                    for (int ix=0;ix<3;++ix) tvx[ix]=hh*tvx[ix]/dist;
                    std::array<double,3> tvy;
                    tvy[0]=(xx[1]*tvx[2]-xx[2]*tvx[1])/dab;
                    tvy[1]=(xx[2]*tvx[0]-xx[0]*tvx[2])/dab;
                    tvy[2]=(xx[0]*tvx[1]-xx[1]*tvx[0])/dab;
                    for (int lsign : {-1,1}) {
                        int il = ntrp+1;
                        std::array<double,3> trp;
                        for (int ix=0;ix<3;++ix)
                            trp[ix]=xd[ix]+cj*tvx[ix]+sj*tvy[ix]*lsign;
                        bool bad=false;
                        for (int ik=ilipa+1;ik<=ilipa+nipa[ia];++ik) {
                            int k=lipa[ik];
                            if (k==ib||k==ic) continue;
                            double dabck=0;
                            for (int ix=0;ix<3;++ix) {
                                double dd=trp[ix]-coord[ix*(numat+1)+k];
                                dabck+=dd*dd;
                            }
                            dabck=std::sqrt(dabck);
                            if (dabck < srad[k]+rsolv) { bad=true; break; }
                        }
                        if (bad) continue;
                        ++ntrp;
                        double spx=0,spy=0;
                        for (int ix=0;ix<3;++ix) {
                            spx+=rvx[ix]*(trp[ix]-xd[ix]);
                            spy+=rvy[ix]*(trp[ix]-xd[ix]);
                        }
                        double phi=std::acos(spx/(hh*hh+1e-10));
                        if (spy<0) phi=-phi;
                        phiset[il]=phi;
                        sp=0;
                        for (int ix=0;ix<3;++ix)
                            sp+=(-spy*rvx[ix]+spx*rvy[ix])*(trp[ix]-xic[ix]);
                        iset[ntrp]= (sp<0)?-1:1;
                        double ee1[3],ee2[3],ee3[3];
                        double sp2=0,dac=0,dbc=0;
                        for (int ix=0;ix<3;++ix) {
                            ee1[ix]=trp[ix]+rsolv/srad[ia]*(coord[ix*(numat+1)+ia]-trp[ix]);
                            ee2[ix]=trp[ix]+rsolv/srad[ib]*(coord[ix*(numat+1)+ib]-trp[ix]);
                            ee3[ix]=trp[ix]+rsolv/srad[ic]*(xic[ix]-trp[ix]);
                            sp2+=(ee1[ix]-ee3[ix])*(ee2[ix]-ee3[ix]);
                            dac+=(ee1[ix]-ee3[ix])*(ee1[ix]-ee3[ix]);
                            dbc+=(ee2[ix]-ee3[ix])*(ee2[ix]-ee3[ix]);
                        }
                        tarset[il]=0.8*std::sqrt(std::max(0.0,dac*dbc-sp2*sp2))/12.0;
                    }
                }
            }
            if (ntrp%2 != 0) mopend("ODD NTRP");
            if (ntrp > 18) mopend("NTRP TOO LARGE");
            if (ntrp == 0) {
                phiset[1]=0; phiset[2]=(2.0*pi); tarset[2]=0;
                iset[1]=1; iset[2]=-1; ntrp=2;
            }
            // insertion sort by phiset
            bool again=true;
            while (again) {
                again=false;
                for (int l=2;l<=ntrp;++l) {
                    if (phiset[l]<phiset[l-1]) {
                        std::swap(phiset[l],phiset[l-1]);
                        std::swap(iset[l],iset[l-1]);
                        std::swap(tarset[l],tarset[l-1]);
                        again=true;
                    }
                }
                if (!again && iset[1]==-1) { phiset[1]+=(2.0*pi); again=true; }
            }
            double sumphi=0;

            int ips0=nrs; (void)ips0;
            for (int l=2;l<=ntrp;l+=2) {
                int k=l-1;
                double phiu=phiset[k], phio=phiset[l];
                int nsa=(int)((phio-phiu)/2/pi*20);
                nsa=std::max(nsa+1,2);
                sumphi += phio-phiu;
                double dp=(phio-phiu)/(nsa-1);
                for (int ich=1;ich<=2;++ich) {

                    int iat=(ich==1)?ib:ia;
                    double htr=(iat==ia)?aar/(2.0*pi):abr/(2.0*pi);
                    double fz=(ich==1)?fz1:fz2;
                    double yx=(ich==1)?yx1:yx2;
                    // third atom index for offset
                    int offi=(ich==1)?ia:ib;
                    for (int ja=ich;ja<=nsa;ja+=2) {

                        int jb=std::max(ja-1,1), jc=std::min(ja+1,nsa);
                        double xja[3],xjb[3],xjc[3];
                        for (int ix=0;ix<3;++ix) {
                            double phi=phiu+(ja-1)*dp;
                            double ca=xd[ix]+(std::cos(phi)*rvx[ix]+std::sin(phi)*rvy[ix])*fz;
                            ca=ca+(coord[ix*(numat+1)+(ich==1?ia:ib)]-ca)*yx;
                            xja[ix]=ca+(coord[ix*(numat+1)+offi]-ca)*ddd*(1.0-fz);
                        }
                        for (int ix=0;ix<3;++ix) {
                            double phi=phiu+(jb-1)*dp;
                            double ca=xd[ix]+std::cos(phi)*rvx[ix]+std::sin(phi)*rvy[ix];
                            xjb[ix]=ca+(coord[ix*(numat+1)+offi]-ca)*yx;
                        }
                        for (int ix=0;ix<3;++ix) {
                            double phi=phiu+(jc-1)*dp;
                            double ca=xd[ix]+std::cos(phi)*rvx[ix]+std::sin(phi)*rvy[ix];
                            xjc[ix]=ca+(coord[ix*(numat+1)+offi]-ca)*yx;
                        }
                        ++nrs;

                        double sp=0,d2=0,spn=0,spn2=0,dist2=0;
                        int iat2=(ich==2)?ia:ib;
                        ipsrs[nrs]=iat2;

                        nipsrs[iat2]=nipsrs[iat2]+1;
                        for (int ix=0;ix<3;++ix) {
                            rsc[nrs][ix]=(xja[ix]*0.5+xjb[ix]+xjc[ix])/2.5;
                            int i2=(ix+1)%3, i3=(i2+1)%3;
                            double cnrs=(xjc[i2]-xjb[i2])*(xja[i3]-xjb[i3])-(xja[i2]-xjb[i2])*(xjc[i3]-xjb[i3]);
                            dist2+=cnrs*cnrs;
                            spn+=cnrs*(rsc[nrs][ix]-coord[ix*(numat+1)+iat2]);
                            spn2+=cnrs*rsc[nrs][ix];
                        }
                        dist2=1.0/std::sqrt(dist2);
                        if (spn<0) dist2=-dist2;
                        rsc[nrs][3]=(jc-jb)*dp*htr;
                        if (ja==1)   rsc[nrs][3]+=tarset[k]*rinc;
                        if (ja==nsa) rsc[nrs][3]+=tarset[l]*rinc;
                        cosvol += rsc[nrs][3]*spn2*dist2;
                        area += rsc[nrs][3];
                    }
                }
            }
            if (sumphi > 1e-10) {
                ++nipc;
                isude[0][nipc]=ia; isude[1][nipc]=ib;
                sumphi/=(2.0*pi);
                sude[0][nipc]=aad*sumphi;
                sude[1][nipc]=abd*sumphi;
            }
        }
        ilipa += nipa[ia];
    }
    if (nrs > maxrs) mopend("NRS .GT. MAXRS IN SURCLO");

    // sort ring segments w.r.t. atoms
    int isum=0;
    for (int iat=1;iat<=numat;++iat) {
        int isum2=isum+nipsrs[iat];
        nipsrs[iat]=isum; isum=isum2;
    }
    for (int i=1;i<=nrs;++i) {
        int iat=ipsrs[i];
        int isum2=nipsrs[iat]+1;
        nipsrs[iat]=isum2; isort[i]=isum2;
    }
    for (int i=1;i<=nrs;++i) {
        while (isort[i]!=i) {
            int is=isort[i];
            for (int ix=0;ix<4;++ix) std::swap(rsc[i][ix],rsc[is][ix]);
            std::swap(ipsrs[i],ipsrs[is]);
            isort[i]=isort[is]; isort[is]=is;
        }
    }
    // match each ring segment to nearest primary segment of same atom
    for (int i=1;i<=nps;++i) nipsrs[i]=0;
    int iat0=0, ips1=0, ips0=0;
    double d2max=0;
    for (int i=1;i<=nrs;++i) {
        int iat=ipsrs[i];
        if (iat>iat0) {
            iat0=iat;
            d2max=16*srad[iat]*srad[iat]/n0[1];
            if (nat[iat]==1) d2max*=n0[1]/n0[2];
            int ips;
            for (ips=ips1+1;ips<=nps;++ips) if (iatsp[ips]==iat) break;
            ips0=ips;
            for (;ips<nps;++ips) if (iatsp[ips+1]>iat) break;
            ips1=ips;
        }
        double d2min=1e6; int ipsmin=ips0;
        for (int ips=ips0;ips<=ips1;++ips) {
            double d2=0;
            for (int ix=0;ix<3;++ix) {
                double dd=cosurf[ix+1][ips]-rsc[i][ix];
                d2+=dd*dd;
            }
            if (d2<d2min) { d2min=d2; ipsmin=ips; }
        }
        if (d2min>d2max) {
            ++ips1;
            for (int ip1=nps;ip1>=ips1;--ip1) {
                for (int ix=0;ix<4;++ix) cosurf[ix+1][ip1+1]=cosurf[ix+1][ip1];
                iatsp[ip1+1]=iatsp[ip1];
                nsetf[ip1+1]=nsetf[ip1];
                nar_csm[ip1+1]=nar_csm[ip1];
                nipsrs[ip1+1]=nipsrs[ip1];
            }
            ipsrs[i]=ips1; ++nps;
            for (int ix=0;ix<4;++ix) cosurf[ix+1][ips1]=rsc[i][ix];
            iatsp[ips1]=iatsp[ips1-1];
            nar_csm[ips1]=0; nipsrs[ips1]=1;
            nsetf[ips1]=nsetf[ips1+1];
        } else {
            ipsrs[i]=ipsmin; ++nipsrs[ipsmin];
            double arseg=cosurf[4][ipsmin], arsegn=arseg+rsc[i][3];
            for (int ix=0;ix<3;++ix)
                cosurf[ix+1][ipsmin]=(arseg*cosurf[ix+1][ipsmin]+rsc[i][3]*rsc[i][ix])/arsegn;
            cosurf[4][ipsmin]=arsegn;
        }
    }
}