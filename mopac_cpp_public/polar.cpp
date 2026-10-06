// polar.cpp - TDHF polarizability driver (MOPAC2016 polar.F90 -> C++)
// Fortran 1-based indexing retained (vectors padded by one element).
#include "polar.h"
#include "polar_helpers.h"
#include "parameters_C.h"
#include "chanel_C.h"
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>

// ---- module globals (extern) ----
namespace common_arrays_C {
extern std::vector<int> nat, na, labels;
extern std::vector<std::vector<double>> geo, coord, c;
extern std::vector<std::vector<double>> tvec;
}
extern int norbs;   // global (molkst_C.cpp:9)
namespace parameters_C {
extern double dd[N_PARAM + 1], gss[N_PARAM + 1], gsp[N_PARAM + 1], gpp[N_PARAM + 1];
extern double gp2[N_PARAM + 1], hsp[N_PARAM + 1], polvol[N_PARAM + 1];
}
namespace funcon_C { extern double ev, a0; }
namespace elemts_C { extern std::vector<std::string> elemnt; }
namespace molkst_C {
extern int numat, norbs, natoms, ndep, nvar, last, id;
extern bool limscf, moperr;
extern std::string keywrd;
}
namespace polar_C { extern double omega; }

namespace common_arrays_C {
extern std::vector<double> eigs;
extern std::vector<int> nfirst, nlast;
}
namespace molkst_C { extern int nclose; }
namespace polar_C { extern double alpavg; }
using common_arrays_C::eigs; using common_arrays_C::nfirst; using common_arrays_C::nlast;
using common_arrays_C::nat; using common_arrays_C::na; using common_arrays_C::labels;
using common_arrays_C::geo; using common_arrays_C::coord; using common_arrays_C::c;
using common_arrays_C::tvec;
using molkst_C::nclose; using molkst_C::numat; using molkst_C::keywrd; using polar_C::alpavg;
using funcon_C::ev; using funcon_C::a0;
using polar_C::omega;
extern double second(int);
extern void hmuf(std::vector<std::vector<double>>&, int, std::vector<std::vector<double>>&,
                 std::vector<int>&, std::vector<int>&, std::vector<int>&, int, int);
extern void makeuf(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&,
                   std::vector<std::vector<double>>&, std::vector<double>&, bool&, int, int,
                   double&, double);
extern void densf(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&,
                  std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, int, int,
                  std::vector<double>&);
extern void ffreq2(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&,
                   std::vector<double>&);
extern void ffreq1(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&,
                   std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, int);
extern void dawrit(std::vector<std::vector<double>>&, int, int);
extern void mopend(const std::string&);
extern void to_screen(const std::string&);


void bdenin(std::vector<std::vector<double>>& bdcon, std::vector<std::vector<double>>& ua,
            std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& c,
            int norbs, int nclose){
    std::vector<double> w1(norbs+1), w2(norbs+1), w3(norbs+1), w4(norbs+1);
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j){
        for(int l=nclose+1;l<=norbs;++l){
            double sum1=0.0, sum2=0.0;
            for(int lp=1;lp<=nclose;++lp){ sum1+=ub[l][lp]*c[j][lp]; sum2+=ua[l][lp]*c[j][lp]; }
            w1[l]=sum1; w2[l]=sum2;
        }
        for(int k=1;k<=nclose;++k){
            double sum3=0.0;
            for(int l=nclose+1;l<=norbs;++l) sum3+=ua[k][l]*w1[l]+ub[k][l]*w2[l];
            w3[k]=sum3;
        }
        for(int l=1;l<=nclose;++l){
            double sum1=0.0, sum2=0.0;
            for(int lp=nclose+1;lp<=norbs;++lp){ sum1+=ub[l][lp]*c[j][lp]; sum2+=ua[l][lp]*c[j][lp]; }
            w1[l]=sum1; w2[l]=sum2;
        }
        for(int k=nclose+1;k<=norbs;++k){
            double sum4=0.0;
            for(int l=1;l<=nclose;++l) sum4+=ua[k][l]*w1[l]+ub[k][l]*w2[l];
            w4[k]=sum4;
        }
        double s3=0.0, s4=0.0;
        for(int k=1;k<=nclose;++k) s3+=w3[k]*c[i][k];
        for(int k=nclose+1;k<=norbs;++k) s4+=w4[k]*c[i][k];
        bdcon[i][j]=s3-s4;
    }
}

void bdenup(std::vector<std::vector<double>>& bdcon, std::vector<std::vector<double>>& uab,
            std::vector<std::vector<double>>& c, std::vector<std::vector<double>>& d,
            std::vector<std::vector<double>>& da, int norbs, int nclose){
    std::vector<double> w1(norbs+1);
    zerom(d, norbs);
    for(int j=1;j<=norbs;++j) for(int k=1;k<=norbs;++k){
        double sum=0.0;
        for(int l=1;l<=nclose;++l) sum+=uab[k][l]*c[j][l];
        da[k][j]=sum;
    }
    for(int i=1;i<=norbs;++i){
        for(int k=1;k<=norbs;++k){
            double sum=0.0;
            for(int l=1;l<=nclose;++l) sum+=c[i][l]*uab[l][k];
            w1[k]=sum;
        }
        for(int j=1;j<=norbs;++j){
            double s1=0.0, s2=0.0;
            for(int k=1;k<=norbs;++k){ s1+=c[i][k]*da[k][j]; s2+=w1[k]*c[j][k]; }
            d[i][j]=2.0*(s1-s2+bdcon[i][j]);
        }
    }
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) da[i][j]=d[i][j]*0.5;
}


void bmakuf(std::vector<std::vector<double>>& ua, std::vector<std::vector<double>>& ub,
            std::vector<std::vector<double>>& uab, std::vector<std::vector<double>>& t,
            std::vector<std::vector<double>>& uold1, std::vector<std::vector<double>>& gab,
            std::vector<double>& eigs_, bool& last, int norbs, int nclose,
            double& diff, int iwflb, double& maxu, double btol){
    double hartr=ev;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=i;++j){
        double sum=0.0;
        int kll=1, kul=0;
        if(i<=nclose){ kll=nclose+1; kul=norbs; }
        else if(i>nclose && j>nclose){ kll=1; kul=nclose; }
        for(int k=kll;k<=kul;++k) sum+=ua[i][k]*ub[k][j]+ub[i][k]*ua[k][j];
        uab[i][j]=sum*0.5;
        uab[j][i]=sum*0.5;
    }
    for(int k=nclose+1;k<=norbs;++k) for(int l=1;l<=nclose;++l){
        if(iwflb==2){
            uab[k][l]=hartr*((gab[k][l]+t[k][l])/((eigs_[l]-eigs_[k])-omega));
            uab[l][k]=hartr*((gab[l][k]+t[l][k])/((eigs_[k]-eigs_[l])-omega));
        } else if(iwflb==3){
            uab[k][l]=hartr*((gab[k][l]+t[k][l])/(eigs_[l]-eigs_[k]));
            uab[l][k]=hartr*((gab[l][k]+t[l][k])/(eigs_[k]-eigs_[l]));
        } else {
            uab[k][l]=hartr*((gab[k][l]+t[k][l])/((eigs_[l]-eigs_[k])-2.0*omega));
            uab[l][k]=hartr*((gab[l][k]+t[l][k])/((eigs_[k]-eigs_[l])-2.0*omega));
        }
    }
    diff=0.0;
    maxu=-1000.0;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j){
        double udif=uab[i][j]-uold1[i][j];
        diff=std::max(std::abs(udif),diff);
        maxu=std::max(uab[i][j],maxu);
    }
    if(diff<btol) last=true;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) uold1[i][j]=uab[i][j];
}

void epsab(std::vector<std::vector<double>>& eigsab, std::vector<double>& eigs_,
           std::vector<std::vector<double>>& gab, std::vector<std::vector<double>>& ga,
           std::vector<std::vector<double>>& gb, std::vector<std::vector<double>>& ua,
           std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& uab,
           std::vector<std::vector<double>>& udms, int norbs, int nclose, int iwflb){
    double hartr=ev;
    zerom(eigsab, norbs);
    zerom(udms, norbs);
    double omval=2.0*omega;
    if(iwflb==3) omval=0.0;
    else if(iwflb==2) omval=omega;
    for(int i=1;i<=nclose;++i) for(int j=1;j<=nclose;++j){
        double s1=0.0;
        for(int k=nclose+1;k<=norbs;++k) s1+=ga[i][k]*ub[k][j]+gb[i][k]*ua[k][j];
        eigsab[i][j]=gab[i][j]+s1+uab[i][j]*(eigs_[i]-eigs_[j]+omval)/hartr;
    }
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j){
        double s2=0.0;
        for(int k=1;k<=norbs;++k) s2+=ua[i][k]*ub[k][j]+ub[i][k]*ua[k][j];
        udms[i][j]=s2-uab[i][j];
    }
}

void fhpatn(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b,
            int norbs, int itw, double sign){
    if(itw==1 || itw==4){
        for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) a[i][j]=b[i][j];
    } else {
        for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) a[i][j]=sign*b[j][i];
    }
}

void nonor(std::vector<std::vector<double>>& u0x, std::vector<std::vector<double>>& u1y,
           std::vector<std::vector<double>>& u1z, std::vector<std::vector<double>>& u1x,
           std::vector<std::vector<double>>& u0y, std::vector<std::vector<double>>& u0z,
           std::vector<std::vector<double>>& g0x, std::vector<std::vector<double>>& g1y,
           std::vector<std::vector<double>>& g1z, std::vector<std::vector<double>>& g1x,
           std::vector<std::vector<double>>& g0y, std::vector<std::vector<double>>& g0z){
    int maxsq=norbs*norbs;
    double bavx=0.0, bavy=0.0, bavz=0.0;
    daread(u0x, maxsq, 2);
    daread(u0y, maxsq, 3);
    daread(u0z, maxsq, 4);
    daread(g0x, maxsq, 5);
    daread(g0y, maxsq, 6);
    daread(g0z, maxsq, 7);
    daread(u1x, maxsq, 8);
    daread(u1y, maxsq, 9);
    daread(u1z, maxsq, 10);
    daread(g1x, maxsq, 11);
    daread(g1y, maxsq, 12);
    daread(g1z, maxsq, 13);
    double bxxx,byxx,bzxx,bxxy,byxy,bzxy,bxxz,byxz,bzxz,bxyx,byyx,bzyx,bxyy,byyy,bzyy;
    double bxyz,byyz,bzyz,bxzx,byzx,bzzx,bxzy,byzy,bzzy,bxzz,byzz,bzzz;
    betal1(u0x,g0x,u1x,g1x,u1x,g1x,nclose,norbs,bxxx); bavx+=3.0*bxxx;
    betal1(u0y,g0y,u1x,g1x,u1x,g1x,nclose,norbs,byxx); bavy+=byxx;
    betal1(u0z,g0z,u1x,g1x,u1x,g1x,nclose,norbs,bzxx); bavz+=bzxx;
    betal1(u0x,g0x,u1x,g1x,u1y,g1y,nclose,norbs,bxxy); bavy+=bxxy;
    betal1(u0y,g0y,u1x,g1x,u1y,g1y,nclose,norbs,byxy); bavx+=byxy;
    betal1(u0z,g0z,u1x,g1x,u1y,g1y,nclose,norbs,bzxy);
    betal1(u0x,g0x,u1x,g1x,u1z,g1z,nclose,norbs,bxxz); bavz+=bxxz;
    betal1(u0y,g0y,u1x,g1x,u1z,g1z,nclose,norbs,byxz);
    betal1(u0z,g0z,u1x,g1x,u1z,g1z,nclose,norbs,bzxz); bavx+=bzxz;
    betal1(u0x,g0x,u1y,g1y,u1x,g1x,nclose,norbs,bxyx); bavy+=bxyx;
    betal1(u0y,g0y,u1y,g1y,u1x,g1x,nclose,norbs,byyx); bavx+=byyx;
    betal1(u0z,g0z,u1y,g1y,u1x,g1x,nclose,norbs,bzyx);
    betal1(u0x,g0x,u1y,g1y,u1y,g1y,nclose,norbs,bxyy); bavx+=bxyy;
    betal1(u0y,g0y,u1y,g1y,u1y,g1y,nclose,norbs,byyy); bavy+=3.0*byyy;
    betal1(u0z,g0z,u1y,g1y,u1y,g1y,nclose,norbs,bzyy); bavz+=bzyy;
    betal1(u0x,g0x,u1y,g1y,u1z,g1z,nclose,norbs,bxyz);
    betal1(u0y,g0y,u1y,g1y,u1z,g1z,nclose,norbs,byyz); bavz+=byyz;
    betal1(u0z,g0z,u1y,g1y,u1z,g1z,nclose,norbs,bzyz); bavy+=bzyz;
    betal1(u0x,g0x,u1z,g1z,u1x,g1x,nclose,norbs,bxzx); bavz+=bxzx;
    betal1(u0y,g0y,u1z,g1z,u1x,g1x,nclose,norbs,byzx);
    betal1(u0z,g0z,u1z,g1z,u1x,g1x,nclose,norbs,bzzx); bavx+=bzzx;
    betal1(u0x,g0x,u1z,g1z,u1y,g1y,nclose,norbs,bxzy);
    betal1(u0y,g0y,u1z,g1z,u1y,g1y,nclose,norbs,byzy); bavz+=byzy;
    betal1(u0z,g0z,u1z,g1z,u1y,g1y,nclose,norbs,bzzy); bavy+=bzzy;
    betal1(u0x,g0x,u1z,g1z,u1z,g1z,nclose,norbs,bxzz); bavx+=bxzz;
    betal1(u0y,g0y,u1z,g1z,u1z,g1z,nclose,norbs,byzz); bavy+=byzz;
    betal1(u0z,g0z,u1z,g1z,u1z,g1z,nclose,norbs,bzzz); bavz+=3.0*bzzz;
    bavx/=5.0; bavy/=5.0; bavz/=5.0;
    double bvec=std::sqrt(bavx*bavx+bavy*bavy+bavz*bavz);
    std::printf("\n\n BETA (OPTICAL RECTIFICATION) \n");
    std::printf("\n\n  BXXX  %15.8e  BYXX %15.8e  BZXX %15.8e\n  BXXY  %15.8e  BYXY %15.8e  BZXY %15.8e\n  BXXZ  %15.8e  BYXZ %15.8e  BZXZ %15.8e\n  BXYX  %15.8e  BYYX %15.8e  BZYX %15.8e\n  BXYY  %15.8e  BYYY %15.8e  BZYY %15.8e\n  BXYZ  %15.8e  BYYZ %15.8e  BZYZ %15.8e\n  BXZX  %15.8e  BYZX %15.8e  BZZX %15.8e\n  BXZY  %15.8e  BYZY %15.8e  BZZY %15.8e\n  BXZZ  %15.8e  BYZZ %15.8e  BZZZ %15.8e\n",
        bxxx,byxx,bzxx,bxxy,byxy,bzxy,bxxz,byxz,bzxz,bxyx,byyx,bzyx,bxyy,byyy,bzyy,
        bxyz,byyz,bzyz,bxzx,byzx,bzzx,bxzy,byzy,bzzy,bxzz,byzz,bzzz);
    std::printf("\n\n AVERAGE BETAX VALUE AT %10.5fEV = %15.5f a.u.\n", omega, bavx);
    std::printf(" AVERAGE BETAY VALUE AT %10.5fEV = %15.5f a.u.\n", omega, bavy);
    std::printf(" AVERAGE BETAZ VALUE AT %10.5fEV = %15.5f a.u.\n\n", omega, bavz);
    std::printf("\n\n AVERAGE BETA(OR) VALUE AT %10.5fEV = %15.5f a.u.\n\n", omega, bvec);
}

void nonope(std::vector<std::vector<double>>& u0x, std::vector<std::vector<double>>& u1y,
            std::vector<std::vector<double>>& u1z, std::vector<std::vector<double>>& u1x,
            std::vector<std::vector<double>>& u0y, std::vector<std::vector<double>>& u0z,
            std::vector<std::vector<double>>& g0x, std::vector<std::vector<double>>& g1y,
            std::vector<std::vector<double>>& g1z, std::vector<std::vector<double>>& g1x,
            std::vector<std::vector<double>>& g0y, std::vector<std::vector<double>>& g0z){
    const double au_to_esu=8.639418e-33;
    int maxsq=norbs*norbs;
    double bavx=0.0, bavy=0.0, bavz=0.0;
    daread(u0x, maxsq, 2);
    daread(u0y, maxsq, 3);
    daread(u0z, maxsq, 4);
    daread(g0x, maxsq, 5);
    daread(g0y, maxsq, 6);
    daread(g0z, maxsq, 7);
    daread(u1x, maxsq, 8);
    daread(u1y, maxsq, 9);
    daread(u1z, maxsq, 10);
    daread(g1x, maxsq, 11);
    daread(g1y, maxsq, 12);
    daread(g1z, maxsq, 13);
    double bxxx,byxx,bzxx,bxxy,byxy,bzxy,bxxz,byxz,bzxz,bxyx,byyx,bzyx,bxyy,byyy,bzyy;
    double bxyz,byyz,bzyz,bxzx,byzx,bzzx,bxzy,byzy,bzzy,bxzz,byzz,bzzz;
    betall(u1x,g1x,u0x,g0x,u1x,g1x,nclose,norbs,bxxx); bavx+=3.0*bxxx;
    betall(u1y,g1y,u0x,g0x,u1x,g1x,nclose,norbs,byxx); bavy+=byxx;
    betall(u1z,g1z,u0x,g0x,u1x,g1x,nclose,norbs,bzxx); bavz+=bzxx;
    betall(u1x,g1x,u0x,g0x,u1y,g1y,nclose,norbs,bxxy); bavy+=bxxy;
    betall(u1y,g1y,u0x,g0x,u1y,g1y,nclose,norbs,byxy); bavx+=byxy;
    betall(u1z,g1z,u0x,g0x,u1y,g1y,nclose,norbs,bzxy);
    betall(u1x,g1x,u0x,g0x,u1z,g1z,nclose,norbs,bxxz); bavz+=bxxz;
    betall(u1y,g1y,u0x,g0x,u1z,g1z,nclose,norbs,byxz);
    betall(u1z,g1z,u0x,g0x,u1z,g1z,nclose,norbs,bzxz); bavx+=bzxz;
    betall(u1x,g1x,u0y,g0y,u1x,g1x,nclose,norbs,bxyx); bavy+=bxyx;
    betall(u1y,g1y,u0y,g0y,u1x,g1x,nclose,norbs,byyx); bavx+=byyx;
    betall(u1z,g1z,u0y,g0y,u1x,g1x,nclose,norbs,bzyx);
    betall(u1x,g1x,u0y,g0y,u1y,g1y,nclose,norbs,bxyy); bavx+=bxyy;
    betall(u1y,g1y,u0y,g0y,u1y,g1y,nclose,norbs,byyy); bavy+=3.0*byyy;
    betall(u1z,g1z,u0y,g0y,u1y,g1y,nclose,norbs,bzyy); bavz+=bzyy;
    betall(u1x,g1x,u0y,g0y,u1z,g1z,nclose,norbs,bxyz);
    betall(u1y,g1y,u0y,g0y,u1z,g1z,nclose,norbs,byyz); bavz+=byyz;
    betall(u1z,g1z,u0y,g0y,u1z,g1z,nclose,norbs,bzyz); bavy+=bzyz;
    betall(u1x,g1x,u0z,g0z,u1x,g1x,nclose,norbs,bxzx); bavz+=bxzx;
    betall(u1y,g1y,u0z,g0z,u1x,g1x,nclose,norbs,byzx);
    betall(u1z,g1z,u0z,g0z,u1x,g1x,nclose,norbs,bzzx); bavx+=bzzx;
    betall(u1x,g1x,u0z,g0z,u1y,g1y,nclose,norbs,bxzy);
    betall(u1y,g1y,u0z,g0z,u1y,g1y,nclose,norbs,byzy); bavz+=byzy;
    betall(u1z,g1z,u0z,g0z,u1y,g1y,nclose,norbs,bzzy); bavy+=bzzy;
    betall(u1x,g1x,u0z,g0z,u1z,g1z,nclose,norbs,bxzz); bavx+=bxzz;
    betall(u1y,g1y,u0z,g0z,u1z,g1z,nclose,norbs,byzz); bavy+=byzz;
    betall(u1z,g1z,u0z,g0z,u1z,g1z,nclose,norbs,bzzz); bavz+=3.0*bzzz;
    bavx/=5.0; bavy/=5.0; bavz/=5.0;
    double bvec=std::sqrt(bavx*bavx+bavy*bavy+bavz*bavz);
    std::printf("  BETA (ELECTOPTIC POCKELS EFFECT) \n");
    std::printf("\n\n  BXXX  %15.8e  BYXX %15.8e  BZXX %15.8e\n  BXXY  %15.8e  BYXY %15.8e  BZXY %15.8e\n  BXXZ  %15.8e  BYXZ %15.8e  BZXZ %15.8e\n  BXYX  %15.8e  BYYX %15.8e  BZYX %15.8e\n  BXYY  %15.8e  BYYY %15.8e  BZYY %15.8e\n  BXYZ  %15.8e  BYYZ %15.8e  BZYZ %15.8e\n  BXZX  %15.8e  BYZX %15.8e  BZZX %15.8e\n  BXZY  %15.8e  BYZY %15.8e  BZZY %15.8e\n  BXZZ  %15.8e  BYZZ %15.8e  BZZZ %15.8e\n",
        bxxx,byxx,bzxx,bxxy,byxy,bzxy,bxxz,byxz,bzxz,bxyx,byyx,bzyx,bxyy,byyy,bzyy,
        bxyz,byyz,bzyz,bxzx,byzx,bzzx,bxzy,byzy,bzzy,bxzz,byzz,bzzz);
    std::printf("\n");
    if(bvec<1e10){
        std::printf(" AVERAGE BETAX      VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n", omega, bavx, bavx*au_to_esu);
        std::printf(" AVERAGE BETAY      VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n", omega, bavy, bavy*au_to_esu);
        std::printf(" AVERAGE BETAZ      VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n\n\n", omega, bavz, bavz*au_to_esu);
    } else {
        std::printf(" AVERAGE BETAX      VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n", omega, bavx, bavx*au_to_esu);
        std::printf(" AVERAGE BETAY      VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n", omega, bavy, bavy*au_to_esu);
        std::printf(" AVERAGE BETAZ      VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n\n\n", omega, bavz, bavz*au_to_esu);
    }
    std::printf("\n");
}

void nonbet(std::vector<std::vector<double>>& u1x, std::vector<std::vector<double>>& u1y,
            std::vector<std::vector<double>>& u1z, std::vector<std::vector<double>>& u2x,
            std::vector<std::vector<double>>& u2y, std::vector<std::vector<double>>& u2z,
            std::vector<std::vector<double>>& g1x, std::vector<std::vector<double>>& g1y,
            std::vector<std::vector<double>>& g1z, std::vector<std::vector<double>>& g2x,
            std::vector<std::vector<double>>& g2y, std::vector<std::vector<double>>& g2z){
    int maxsq=norbs*norbs;
    double bavx=0.0, bavy=0.0, bavz=0.0;
    daread(u1x, maxsq, 8);
    daread(u1y, maxsq, 9);
    daread(u1z, maxsq, 10);
    daread(g1x, maxsq, 11);
    daread(g1y, maxsq, 12);
    daread(g1z, maxsq, 13);
    daread(u2x, maxsq, 14);
    daread(u2y, maxsq, 15);
    daread(u2z, maxsq, 16);
    daread(g2x, maxsq, 17);
    daread(g2y, maxsq, 18);
    daread(g2z, maxsq, 19);
    double bxxx,byxx,bzxx,bxxy,byxy,bzxy,bxxz,byxz,bzxz,bxyx,byyx,bzyx,bxyy,byyy,bzyy;
    double bxyz,byyz,bzyz,bxzx,byzx,bzzx,bxzy,byzy,bzzy,bxzz,byzz,bzzz;
    betcom(u1x,g1x,u2x,g2x,nclose,norbs,bxxx); bavx+=3.0*bxxx;
    betcom(u1x,g1x,u2y,g2y,nclose,norbs,byxx); bavy+=byxx;
    betcom(u1x,g1x,u2z,g2z,nclose,norbs,bzxx); bavz+=bzxx;
    betall(u2x,g2x,u1x,g1x,u1y,g1y,nclose,norbs,bxxy); bavy+=bxxy;
    betall(u2y,g2y,u1x,g1x,u1y,g1y,nclose,norbs,byxy); bavx+=byxy;
    betall(u2z,g2z,u1x,g1x,u1y,g1y,nclose,norbs,bzxy);
    betall(u2x,g2x,u1x,g1x,u1z,g1z,nclose,norbs,bxxz); bavz+=bxxz;
    betall(u2y,g2y,u1x,g1x,u1z,g1z,nclose,norbs,byxz);
    betall(u2z,g2z,u1x,g1x,u1z,g1z,nclose,norbs,bzxz); bavx+=bzxz;
    betall(u2x,g2x,u1y,g1y,u1x,g1x,nclose,norbs,bxyx); bavy+=bxyx;
    betall(u2y,g2y,u1y,g1y,u1x,g1x,nclose,norbs,byyx); bavx+=byyx;
    betall(u2z,g2z,u1y,g1y,u1x,g1x,nclose,norbs,bzyx);
    betcom(u1y,g1y,u2x,g2x,nclose,norbs,bxyy); bavx+=bxyy;
    betcom(u1y,g1y,u2y,g2y,nclose,norbs,byyy); bavy+=3.0*byyy;
    betcom(u1y,g1y,u2z,g2z,nclose,norbs,bzyy); bavz+=bzyy;
    betall(u2x,g2x,u1y,g1y,u1z,g1z,nclose,norbs,bxyz);
    betall(u2y,g2y,u1y,g1y,u1z,g1z,nclose,norbs,byyz); bavz+=byyz;
    betall(u2z,g2z,u1y,g1y,u1z,g1z,nclose,norbs,bzyz); bavy+=bzyz;
    betall(u2x,g2x,u1z,g1z,u1x,g1x,nclose,norbs,bxzx); bavz+=bxzx;
    betall(u2y,g2y,u1z,g1z,u1x,g1x,nclose,norbs,byzx);
    betall(u2z,g2z,u1z,g1z,u1x,g1x,nclose,norbs,bzzx); bavx+=bzzx;
    betall(u2x,g2x,u1z,g1z,u1y,g1y,nclose,norbs,bxzy);
    betall(u2y,g2y,u1z,g1z,u1y,g1y,nclose,norbs,byzy); bavz+=byzy;
    betall(u2z,g2z,u1z,g1z,u1y,g1y,nclose,norbs,bzzy); bavy+=bzzy;
    betcom(u1z,g1z,u2x,g2x,nclose,norbs,bxzz); bavx+=bxzz;
    betcom(u1z,g1z,u2y,g2y,nclose,norbs,byzz); bavy+=byzz;
    betcom(u1z,g1z,u2z,g2z,nclose,norbs,bzzz); bavz+=3.0*bzzz;
    bavx/=5.0; bavy/=5.0; bavz/=5.0;
    double bvec=std::sqrt(bavx*bavx+bavy*bavy+bavz*bavz);
    std::printf("\n\n BETA (SECOND HARMONIC GENERATION)\n\n");
    std::printf("\n\n  BXXX  %15.8e  BYXX %15.8e  BZXX %15.8e\n  BXXY  %15.8e  BYXY %15.8e  BZXY %15.8e\n  BXXZ  %15.8e  BYXZ %15.8e  BZXZ %15.8e\n  BXYX  %15.8e  BYYX %15.8e  BZYX %15.8e\n  BXYY  %15.8e  BYYY %15.8e  BZYY %15.8e\n  BXYZ  %15.8e  BYYZ %15.8e  BZYZ %15.8e\n  BXZX  %15.8e  BYZX %15.8e  BZZX %15.8e\n  BXZY  %15.8e  BYZY %15.8e  BZZY %15.8e\n  BXZZ  %15.8e  BYZZ %15.8e  BZZZ %15.8e\n",
        bxxx,byxx,bzxx,bxxy,byxy,bzxy,bxxz,byxz,bzxz,bxyx,byyx,bzyx,bxyy,byyy,bzyy,
        bxyz,byyz,bzyz,bxzx,byzx,bzzx,bxzy,byzy,bzzy,bxzz,byzz,bzzz);
    std::printf("\n");
    if(bvec<1e10){
        std::printf(" AVERAGE BETA X (SHG) VALUE AT%10.5f EV = %11.4f a.u.\n", omega, bavx);
        std::printf(" AVERAGE BETA Y (SHG) VALUE AT%10.5f EV = %11.4f a.u.\n", omega, bavy);
        std::printf(" AVERAGE BETA Z (SHG) VALUE AT%10.5f EV = %11.4f a.u.\n\n", omega, bavz);
        std::printf(" AVERAGE BETA   (SHG) VALUE AT%10.5f EV = %11.4f a.u.\n", omega, bvec);
    } else if(bvec<1e15){
        std::printf(" AVERAGE BETA X (SHG) VALUE AT%10.5f EV = %16.4f a.u.\n", omega, bavx);
        std::printf(" AVERAGE BETA Y (SHG) VALUE AT%10.5f EV = %16.4f a.u.\n", omega, bavy);
        std::printf(" AVERAGE BETA Z (SHG) VALUE AT%10.5f EV = %16.4f a.u.\n\n", omega, bavz);
        std::printf(" AVERAGE BETA   (SHG) VALUE AT%10.5f EV = %16.4f a.u.\n", omega, bvec);
    } else {
        std::printf(" AVERAGE BETA X (SHG) VALUE AT%10.5f EV = %21.4f a.u.\n", omega, bavx);
        std::printf(" AVERAGE BETA Y (SHG) VALUE AT%10.5f EV = %21.4f a.u.\n", omega, bavy);
        std::printf(" AVERAGE BETA Z (SHG) VALUE AT%10.5f EV = %21.4f a.u.\n\n", omega, bavz);
        std::printf(" AVERAGE BETA   (SHG) VALUE AT%10.5f EV = %21.4f a.u.\n", omega, bvec);
    }
    std::printf("\n");
}

void ngoke(int igam, std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& gd3,
           std::vector<std::vector<double>>& ud3, std::vector<std::vector<double>>& g1,
           std::vector<std::vector<double>>& u1, std::vector<std::vector<double>>& gs,
           std::vector<std::vector<double>>& usmd, std::vector<std::vector<double>>& eps,
           std::vector<std::vector<double>>& us){
    const int ida[16]={0,1,2,3,1,1,2,2,3,3,1,1,2,2,3,3};
    const int idb[16]={0,1,2,3,1,1,2,2,3,3,2,3,1,3,1,2};
    const int idc[16]={0,1,2,3,2,3,1,3,1,2,2,3,1,3,1,2};
    const int idd[16]={0,1,2,3,2,3,1,3,1,2,1,1,2,2,3,3};
    const int ip[4][4]={{0,0,0,0},{0,1,2,3},{0,2,4,5},{0,3,5,6}};
    const int ipair[4][4]={{0,0,0,0},{0,1,4,7},{0,2,5,8},{0,3,6,9}};
    const char* alab[4]={"","X","Y","Z"};
    double one=1.0;
    int msq=norbs*norbs;
    if(igam==3) std::printf("\n\n GAMMA (IDRI) AT %10.5f EV.\n\n", omega);
    else std::printf("\n\n GAMMA (OKE) AT %10.5f EV.\n\n", omega);
    int jgarc=10, juarc=7, jurec=1, jgrec=4;
    double gav=0.0;
    double gamma[16];
    for(int ie=1;ie<=15;++ie){
        int ia=ida[ie], ib=idb[ie], ic=idc[ie], id=idd[ie];
        int icd=ipair[ic][id], ibd=ipair[ib][id], ibc=ip[ib][ic];
        daread(x, msq, jgarc+ia);
        fhpatn(gd3, x, norbs, 2, one);
        daread(x, msq, juarc+ia);
        fhpatn(ud3, x, norbs, 2, (-one));
        double yy=0.0;
        for(int imove=1;imove<=3;++imove){
            int j2, j34, jg2rec, ju2rec, ju2mrc, jeprec;
            if(imove==2){ j2=ic; j34=ibd; jg2rec=82; ju2rec=73; ju2mrc=100; jeprec=91; }
            else if(imove==3){ j2=id; j34=ibc; jg2rec=31; ju2rec=25; ju2mrc=43; jeprec=37; }
            else { j2=ib; j34=icd; jg2rec=82; ju2rec=73; ju2mrc=100; jeprec=91; }
            if(imove==3) daread(u1, msq, juarc+j2);
            else daread(u1, msq, jurec+j2);
            if(imove==3) daread(g1, msq, jgarc+j2);
            else daread(g1, msq, jgrec+j2);
            daread(gs, msq, jg2rec+j34);
            daread(us, msq, ju2rec+j34);
            daread(usmd, msq, ju2mrc+j34);
            daread(eps, msq, jeprec+j34);
            yy+=trsub(ud3,g1,us,nclose,norbs,norbs);
            yy-=trsub(usmd,g1,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,g1,us,norbs,nclose,norbs);
            yy+=trsub(usmd,g1,ud3,norbs,nclose,norbs);
            yy+=trsub(ud3,gs,u1,nclose,norbs,norbs);
            yy+=trsub(u1,gs,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,eps,u1,norbs,nclose,norbs);
            yy-=trsub(u1,eps,ud3,norbs,nclose,norbs);
            yy+=trsub(u1,gd3,us,nclose,norbs,norbs);
            yy-=trsub(usmd,gd3,u1,nclose,norbs,norbs);
            yy-=trsub(u1,gd3,us,norbs,nclose,norbs);
            yy+=trsub(usmd,gd3,u1,norbs,nclose,norbs);
        }
        gamma[ie]=yy;
        if(ie<=3) gav+=3.0*yy;
        else if(ie>9) gav+=yy;
        else gav+=2.0*yy;
        std::printf(" GAMMA(%c,%c,%c,%c) = %14.7e\n", alab[ia][0], alab[ib][0], alab[ic][0], alab[id][0], gamma[ie]);
    }
    double gave=gav/15.0;
    std::printf("\n\n  AVERAGE GAMMA VALUE AT %10.5f = %14.7e a.u. = %14.7e ESU (X10-39)\n\n", omega, gave, gave/1.985425939776565);
}

void ngidri(std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& gd3,
            std::vector<std::vector<double>>& ud3, std::vector<std::vector<double>>& g1,
            std::vector<std::vector<double>>& u1, std::vector<std::vector<double>>& gs,
            std::vector<std::vector<double>>& usmd, std::vector<std::vector<double>>& eps,
            std::vector<std::vector<double>>& us){
    const int ida[16]={0,1,2,3,1,1,2,2,3,3,1,1,2,2,3,3};
    const int idb[16]={0,1,2,3,1,1,2,2,3,3,2,3,1,3,1,2};
    const int idc[16]={0,1,2,3,2,3,1,3,1,2,2,3,1,3,1,2};
    const int idd[16]={0,1,2,3,2,3,1,3,1,2,1,1,2,2,3,3};
    const int ip[4][4]={{0,0,0,0},{0,1,2,3},{0,2,4,5},{0,3,5,6}};
    const int ipair[4][4]={{0,0,0,0},{0,1,4,7},{0,2,5,8},{0,3,6,9}};
    const char* alab[4]={"","X","Y","Z"};
    double one=1.0;
    int msq=norbs*norbs;
    std::printf("\n\n GAMMA (IDRI) AT %10.5f EV.\n\n", omega);
    int jgarc=10, juarc=7, jurec=7, jgrec=10;
    double gav=0.0;
    double gamma[16];
    for(int ie=1;ie<=15;++ie){
        int ia=ida[ie], ib=idb[ie], ic=idc[ie], id=idd[ie];
        int icd=ipair[ic][id], ibd=ipair[ib][id], ibc=ip[ib][ic];
        daread(x, msq, jgarc+ia);
        fhpatn(gd3, x, norbs, 2, one);
        daread(x, msq, juarc+ia);
        fhpatn(ud3, x, norbs, 2, (-one));
        double yy=0.0;
        for(int imove=1;imove<=3;++imove){
            int j2, j34, jg2rec, ju2rec, ju2mrc, jeprec;
            if(imove==2){ j2=ic; j34=ibd; jg2rec=118; ju2rec=109; ju2mrc=136; jeprec=127; }
            else if(imove==3){ j2=id; j34=ibc; jg2rec=55; ju2rec=49; ju2mrc=67; jeprec=61; }
            else { j2=ib; j34=icd; jg2rec=118; ju2rec=109; ju2mrc=136; jeprec=127; }
            if(imove==3){
                daread(x, msq, jurec+j2);
                fhpatn(u1, x, norbs, 2, (-one));
            } else {
                daread(u1, msq, jurec+j2);
            }
            if(imove==3){
                daread(x, msq, jgrec+j2);
                fhpatn(g1, x, norbs, 2, one);
            } else {
                daread(g1, msq, jgrec+j2);
            }
            daread(gs, msq, jg2rec+j34);
            daread(us, msq, ju2rec+j34);
            daread(usmd, msq, ju2mrc+j34);
            daread(eps, msq, jeprec+j34);
            yy+=trsub(ud3,g1,us,nclose,norbs,norbs);
            yy-=trsub(usmd,g1,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,g1,us,norbs,nclose,norbs);
            yy+=trsub(usmd,g1,ud3,norbs,nclose,norbs);
            yy+=trsub(ud3,gs,u1,nclose,norbs,norbs);
            yy+=trsub(u1,gs,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,eps,u1,norbs,nclose,norbs);
            yy-=trsub(u1,eps,ud3,norbs,nclose,norbs);
            yy+=trsub(u1,gd3,us,nclose,norbs,norbs);
            yy-=trsub(usmd,gd3,u1,nclose,norbs,norbs);
            yy-=trsub(u1,gd3,us,norbs,nclose,norbs);
            yy+=trsub(usmd,gd3,u1,norbs,nclose,norbs);
        }
        gamma[ie]=yy;
        if(ie<=3) gav+=3.0*yy;
        else if(ie>9) gav+=yy;
        else gav+=2.0*yy;
        std::printf(" GAMMA(%c,%c,%c,%c) = %14.7e\n", alab[ia][0], alab[ib][0], alab[ic][0], alab[id][0], gamma[ie]);
    }
    double gave=gav/15.0;
    std::printf("\n\n  AVERAGE GAMMA VALUE AT %10.5f = %14.7e a.u. = %14.7e ESU (X10-39)\n\n", omega, gave, gave/1.985425939776565);
}

void ngefis(std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& gd3,
            std::vector<std::vector<double>>& ud3, std::vector<std::vector<double>>& g1,
            std::vector<std::vector<double>>& u1, std::vector<std::vector<double>>& gs,
            std::vector<std::vector<double>>& usmd, std::vector<std::vector<double>>& eps,
            std::vector<std::vector<double>>& us){
    const int ida[16]={0,1,2,3,1,1,2,2,3,3,1,1,2,2,3,3};
    const int idb[16]={0,1,2,3,2,3,1,3,1,2,1,1,2,2,3,3};
    const int idc[16]={0,1,2,3,1,1,2,2,3,3,2,3,1,3,1,2};
    const int idd[16]={0,1,2,3,2,3,1,3,1,2,2,3,1,3,1,2};
    const int ip[4][4]={{0,0,0,0},{0,1,2,3},{0,2,4,5},{0,3,5,6}};
    const int ipair[4][4]={{0,0,0,0},{0,1,4,7},{0,2,5,8},{0,3,6,9}};
    const char* alab[4]={"","X","Y","Z"};
    double one=1.0;
    int msq=norbs*norbs;
    std::printf("\n\n GAMMA (DC-EFISHG) AT %10.5f EV.\n\n", omega);
    int jgarc=16, juarc=13, jurec=1, jgrec=4;
    int jg2rec=55, ju2rec=49, ju2mrc=67, jeprec=61;
    double gav=0.0;
    double gamma[16];
    for(int ie=1;ie<=15;++ie){
        int ia=ida[ie], ib=idb[ie], ic=idc[ie], id=idd[ie];
        int icd=ip[ic][id], ibd=ipair[ib][id], ibc=ipair[ib][ic];
        daread(x, msq, jgarc+ia);
        fhpatn(gd3, x, norbs, 2, one);
        daread(x, msq, juarc+ia);
        fhpatn(ud3, x, norbs, 2, (-one));
        double yy=0.0;
        for(int imove=1;imove<=3;++imove){
            int j2, j34=0, j3u=0, j3g=0, j3e=0, j3um=0;
            if(imove==1){ j2=ib; j34=icd; }
            else if(imove==2){ j2=ic+6; j3u=ibd+24; j3g=ibd+27; j3e=ibd+30; j3um=ibd+33; }
            else { j2=id+6; j3u=ibc+24; j3g=ibc+27; j3e=ibc+30; j3um=ibc+33; }
            daread(u1, msq, jurec+j2);
            daread(g1, msq, jgrec+j2);
            if(imove==1){
                daread(gs, msq, jg2rec+j34);
                daread(us, msq, ju2rec+j34);
                daread(usmd, msq, ju2mrc+j34);
                daread(eps, msq, jeprec+j34);
            } else {
                daread(gs, msq, jg2rec+j3g);
                daread(us, msq, ju2rec+j3u);
                daread(usmd, msq, ju2mrc+j3um);
                daread(eps, msq, jeprec+j3e);
            }
            yy+=trsub(ud3,g1,us,nclose,norbs,norbs);
            yy-=trsub(usmd,g1,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,g1,us,norbs,nclose,norbs);
            yy+=trsub(usmd,g1,ud3,norbs,nclose,norbs);
            yy+=trsub(ud3,gs,u1,nclose,norbs,norbs);
            yy+=trsub(u1,gs,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,eps,u1,norbs,nclose,norbs);
            yy-=trsub(u1,eps,ud3,norbs,nclose,norbs);
            yy+=trsub(u1,gd3,us,nclose,norbs,norbs);
            yy-=trsub(usmd,gd3,u1,nclose,norbs,norbs);
            yy-=trsub(u1,gd3,us,norbs,nclose,norbs);
            yy+=trsub(usmd,gd3,u1,norbs,nclose,norbs);
        }
        gamma[ie]=yy;
        if(ie<=3) gav+=3.0*yy;
        else if(ie>9) gav+=yy;
        else gav+=2.0*yy;
        std::printf(" GAMMA(%c,%c,%c,%c) = %14.7e\n", alab[ia][0], alab[ib][0], alab[ic][0], alab[id][0], gamma[ie]);
    }
    double gave=gav/15.0;
    std::printf("\n\n AVERAGE GAMMA VALUE AT %10.5f EV = %14.7e a.u. = %14.7e ESU (X10-39)\n\n", omega, gave, gave/1.985425939776565);
}

void ngamtg(std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& gd3,
            std::vector<std::vector<double>>& ud3, std::vector<std::vector<double>>& g1,
            std::vector<std::vector<double>>& u1, std::vector<std::vector<double>>& gs,
            std::vector<std::vector<double>>& usmd, std::vector<std::vector<double>>& eps,
            std::vector<std::vector<double>>& us){
    const int ida[10]={0,1,2,3,1,1,2,2,3,3};
    const int idb[10]={0,1,2,3,1,1,2,2,3,3};
    const int idc[10]={0,1,2,3,2,3,1,3,1,2};
    const int idd[10]={0,1,2,3,2,3,1,3,1,2};
    const int ipair[4][4]={{0,0,0,0},{0,1,2,3},{0,2,4,5},{0,3,5,6}};
    const char* alab[4]={"","X","Y","Z"};
    double one=1.0;
    int msq=norbs*norbs;
    std::printf("\n\n GAMMA (THIRD HARMONIC GENERATION) AT %10.5f EV.\n\n", omega);
    int jgarc=22, juarc=19, jurec=7, jgrec=10;
    int jg2rec=55, ju2rec=49, ju2mrc=67, jeprec=61;
    double gav=0.0;
    double gamma[10];
    for(int ie=1;ie<=9;++ie){
        int ia=ida[ie], ib=idb[ie], ic=idc[ie], id=idd[ie];
        int icd=ipair[ic][id], ibd=ipair[ib][id], ibc=ipair[ib][ic];
        daread(x, msq, jgarc+ia);
        fhpatn(gd3, x, norbs, 2, one);
        daread(x, msq, juarc+ia);
        fhpatn(ud3, x, norbs, 2, (-one));
        double yy=0.0;
        for(int imove=1;imove<=3;++imove){
            int j2, j34;
            if(imove==2){ j2=ic; j34=ibd; }
            else if(imove==3){ j2=id; j34=ibc; }
            else { j2=ib; j34=icd; }
            daread(u1, msq, jurec+j2);
            daread(g1, msq, jgrec+j2);
            daread(gs, msq, jg2rec+j34);
            daread(us, msq, ju2rec+j34);
            daread(usmd, msq, ju2mrc+j34);
            daread(eps, msq, jeprec+j34);
            yy+=trsub(ud3,g1,us,nclose,norbs,norbs);
            yy-=trsub(usmd,g1,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,g1,us,norbs,nclose,norbs);
            yy+=trsub(usmd,g1,ud3,norbs,nclose,norbs);
            yy+=trsub(ud3,gs,u1,nclose,norbs,norbs);
            yy+=trsub(u1,gs,ud3,nclose,norbs,norbs);
            yy-=trsub(ud3,eps,u1,norbs,nclose,norbs);
            yy-=trsub(u1,eps,ud3,norbs,nclose,norbs);
            yy+=trsub(u1,gd3,us,nclose,norbs,norbs);
            yy-=trsub(usmd,gd3,u1,nclose,norbs,norbs);
            yy-=trsub(u1,gd3,us,norbs,nclose,norbs);
            yy+=trsub(usmd,gd3,u1,norbs,nclose,norbs);
        }
        gamma[ie]=yy;
        gav+=yy;
        std::printf(" GAMMA(%c,%c,%c,%c) = %19.5f\n", alab[ia][0], alab[ib][0], alab[ic][0], alab[id][0], gamma[ie]);
    }
    double gave=gav/5.0;
    std::printf("\n\n AVERAGE GAMMA VALUE AT %10.5f = %19.5f a.u. = %19.5f ESU (X10-39)\n\n", omega, gave, gave/1.985425939776565);
}

void beopor(int iwflb, int maxitu, double btol, std::vector<std::vector<double>>& ua,
            std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& f,
            std::vector<std::vector<double>>& ga, std::vector<std::vector<double>>& gb,
            std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& h1,
            std::vector<std::vector<double>>& d, std::vector<std::vector<double>>& da,
            std::vector<std::vector<double>>& uab, std::vector<std::vector<double>>& uold1,
            std::vector<std::vector<double>>& g, std::vector<std::vector<double>>& x,
            std::vector<std::vector<double>>& c, std::vector<double>& w){
    const int ida[10]={0,1,1,1,2,2,2,3,3,3};
    const int idb[10]={0,1,2,3,1,2,3,1,2,3};
    const char* alab[4]={"","X","Y","Z"};
    double one=1.0;
    int maxsq=norbs*norbs;
    int iposu = (iwflb==2) ? 73 : 109;
    int iposg=iposu+9, ipose=iposg+9, iposum=ipose+9;
    if(iwflb==0) std::printf("\n +++++ BETA (STATIC) AT %15.5f EV.\n", omega);
    else if(iwflb==2) std::printf("\n +++++ BETA (ELECTROOPTIC POCKELS EFFECT) AT %15.5f EV.\n", omega);
    else std::printf("\n +++++ BETA (OPTICAL RECTIFICATION) AT %15.5f EV.\n", omega);
    double bavx=0.0, bavy=0.0, bavz=0.0;
    for(int id=1;id<=9;++id){
        double cmptim=second(1);
        int ia=ida[id], ib=idb[id];
        bool last=false;
        hmuf(h1, ia, coord, nfirst, nlast, nat, norbs, numat);
        zerom(uold1, norbs);
        zerom(uab, norbs);
        zerom(f, norbs);
        int jpu, jpg;
        if(iwflb==2 || iwflb==0){
            jpu=1+ia; daread(ua, maxsq, jpu);
            jpg=4+ia; daread(ga, maxsq, jpg);
        } else {
            jpu=7+ia; daread(ua, maxsq, jpu);
            jpg=10+ia; daread(ga, maxsq, jpg);
        }
        if(iwflb==3){
            jpu=7+ib; daread(x, maxsq, jpu); fhpatn(ub, x, norbs, 2, (-one));
            jpg=10+ib; daread(x, maxsq, jpg); fhpatn(gb, x, norbs, 2, one);
        } else if(iwflb==0){
            jpu=1+ib; daread(ub, maxsq, jpu);
            jpg=4+ib; daread(gb, maxsq, jpg);
        } else {
            jpu=7+ib; daread(ub, maxsq, jpu);
            jpg=10+ib; daread(gb, maxsq, jpg);
        }
        tf(ua, ga, ub, gb, t, norbs);
        bdenin(x, ua, ub, c, norbs, nclose);
        bdenup(x, uab, c, d, da, norbs, nclose);
        double betaw=aval(h1,d,norbs);
        ffreq2(f, d, w);
        ffreq1(f, d, da, da, norbs);
        zerom(da, norbs);
        hplusf(f, da, norbs);
        int icount=0;
        double diff=0.0, maxu=0.0;
        ++icount;
        if(icount>=maxitu) last=true;
        transf(f, g, c, norbs);
        bmakuf(ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, diff, iwflb, maxu, btol);
        bdenin(x, ua, ub, c, norbs, nclose);
        bdenup(x, uab, c, d, da, norbs, nclose);
        betaw=aval(h1,d,norbs);
        zerom(f, norbs);
        ffreq2(f, d, w);
        ffreq1(f, d, da, da, norbs);
        zerom(da, norbs);
        hplusf(f, da, norbs);
        while(!last){
            ++icount;
            if(icount>=maxitu) last=true;
            transf(f, g, c, norbs);
            bmakuf(ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, diff, iwflb, maxu, btol);
            bdenin(x, ua, ub, c, norbs, nclose);
            bdenup(x, uab, c, d, da, norbs, nclose);
            betaw=aval(h1,d,norbs);
            zerom(f, norbs);
            ffreq2(f, d, w);
            ffreq1(f, d, da, da, norbs);
            zerom(da, norbs);
            hplusf(f, da, norbs);
        }
        cmptim=second(1)-cmptim;
        std::printf("\n CONVERGED IN%4d ITERATIONS IN%10.2f SECONDS\n", icount, cmptim);
        std::printf(" MAXIMUM UAB ELEMENT =%15.5f,  MAXIMUM DIFFERENCE =%15.5f\n\n", maxu, diff);
        for(int ic=1;ic<=3;++ic){
            hmuf(h1, ic, coord, nfirst, nlast, nat, norbs, numat);
            betaw=aval(h1,d,norbs);
            std::printf("      BETA(%c,%c,%c) = %15.5f\n", alab[ic][0], alab[ia][0], alab[ib][0], betaw);
            if(id==1 && ic==1) bavx+=3.0*betaw;
            else if((id==5 || id==9) && ic==1) bavx+=betaw;
            else if((id==2 || id==4) && ic==2) bavx+=betaw;
            else if((id==3 || id==7) && ic==3) bavx+=betaw;
            if(id==5 && ic==2) bavy+=3.0*betaw;
            else if((id==2 || id==4) && ic==1) bavy+=betaw;
            else if((id==1 || id==9) && ic==2) bavy+=betaw;
            else if((id==6 || id==8) && ic==3) bavy+=betaw;
            if(id==9 && ic==3) bavz+=3.0*betaw;
            else if((id==3 || id==7) && ic==1) bavz+=betaw;
            else if((id==6 || id==8) && ic==2) bavz+=betaw;
            else if((id==1 || id==5) && ic==3) bavz+=betaw;
        }
        epsab(h1, eigs, g, ga, gb, ua, ub, uab, da, norbs, nclose, iwflb);
        dawrit(uab, maxsq, iposu+id);
        dawrit(g, maxsq, iposg+id);
        dawrit(h1, maxsq, ipose+id);
        dawrit(da, maxsq, iposum+id);
    }
    bavx/=5.0; bavy/=5.0; bavz/=5.0;
    double bvec=std::sqrt(bavx*bavx+bavy*bavy+bavz*bavz);
    std::printf("\n");
    if(bvec<1e10){
        std::printf(" AVERAGE BETA X  VALUE AT%10.5f EV = %11.4f\n", omega, bavx);
        std::printf(" AVERAGE BETA Y  VALUE AT%10.5f EV = %11.4f\n", omega, bavy);
        std::printf(" AVERAGE BETA Z  VALUE AT%10.5f EV = %11.4f\n\n", omega, bavz);
        std::printf(" AVERAGE BETA    VALUE AT%10.5f EV = %11.4f\n", omega, bvec);
    } else if(bvec<1e15){
        std::printf(" AVERAGE BETA X  VALUE AT%10.5f EV = %16.4f\n", omega, bavx);
        std::printf(" AVERAGE BETA Y  VALUE AT%10.5f EV = %16.4f\n", omega, bavy);
        std::printf(" AVERAGE BETA Z  VALUE AT%10.5f EV = %16.4f\n\n", omega, bavz);
        std::printf(" AVERAGE BETA    VALUE AT%10.5f EV = %16.4f\n", omega, bvec);
    } else {
        std::printf(" AVERAGE BETA X  VALUE AT%10.5f EV = %21.4f\n", omega, bavx);
        std::printf(" AVERAGE BETA Y  VALUE AT%10.5f EV = %21.4f\n", omega, bavy);
        std::printf(" AVERAGE BETA Z  VALUE AT%10.5f EV = %21.4f\n\n", omega, bavz);
        std::printf(" AVERAGE BETA    VALUE AT%10.5f EV = %21.4f\n", omega, bvec);
    }
    std::printf("\n");
}

void betaf(int iwflb, int maxitu, double btol, std::vector<std::vector<double>>& ua,
           std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& f,
           std::vector<std::vector<double>>& ga, std::vector<std::vector<double>>& gb,
           std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& h1,
           std::vector<std::vector<double>>& d, std::vector<std::vector<double>>& da,
           std::vector<std::vector<double>>& uab, std::vector<std::vector<double>>& uold1,
           std::vector<std::vector<double>>& g, std::vector<std::vector<double>>& x,
           std::vector<std::vector<double>>& c, std::vector<double>& w){
    const double au_to_esu=8.639418e-33;
    const int ida[7]={0,1,1,1,2,2,3};
    const int idb[7]={0,1,2,3,2,3,3};
    const char* alab[4]={"","X","Y","Z"};
    bool debug=keywrd.find(" BETAF ")!=std::string::npos;
    double one=1.0;
    int maxsq=norbs*norbs;
    int iposu=25+24*iwflb;
    int iposg=iposu+6, ipose=iposg+6, iposum=ipose+6;
    if(iwflb==0) std::printf("\n +++++ BETA (STATIC) AT %15.5f EV.\n", omega);
    else std::printf("\n +++++ BETA (SECOND HARMONIC GENERATION) AT %13.5f EV.\n", omega);
    double bavx=0.0, bavy=0.0, bavz=0.0;
    std::vector<std::vector<std::vector<double>>> allbet(4, std::vector<std::vector<double>>(4, std::vector<double>(4,0.0)));
    for(int id=1;id<=6;++id){
        double cmptim=second(1);
        int ia=ida[id], ib=idb[id];
        bool last=false;
        hmuf(h1, ia, coord, nfirst, nlast, nat, norbs, numat);
        zerom(uold1, norbs);
        zerom(uab, norbs);
        zerom(f, norbs);
        int jpu, jpg;
        if(iwflb==2 || iwflb==0){
            jpu=1+ia; daread(ua, maxsq, jpu);
            jpg=4+ia; daread(ga, maxsq, jpg);
        } else {
            jpu=7+ia; daread(ua, maxsq, jpu);
            jpg=10+ia; daread(ga, maxsq, jpg);
        }
        if(iwflb==3){
            jpu=7+ib; daread(x, maxsq, jpu); fhpatn(ub, x, norbs, 2, (-one));
            jpg=10+ib; daread(x, maxsq, jpg); fhpatn(gb, x, norbs, 2, one);
        } else if(iwflb==0){
            jpu=1+ib; daread(ub, maxsq, jpu);
            jpg=4+ib; daread(gb, maxsq, jpg);
        } else {
            jpu=7+ib; daread(ub, maxsq, jpu);
            jpg=10+ib; daread(gb, maxsq, jpg);
        }
        tf(ua, ga, ub, gb, t, norbs);
        bdenin(x, ua, ub, c, norbs, nclose);
        bdenup(x, uab, c, d, da, norbs, nclose);
        double betaw=aval(h1,d,norbs);
        ffreq2(f, d, w);
        ffreq1(f, d, da, da, norbs);
        zerom(da, norbs);
        hplusf(f, da, norbs);
        int icount=0;
        double diff=0.0, maxu=0.0;
        ++icount;
        if(icount>=maxitu) last=true;
        transf(f, g, c, norbs);
        bmakuf(ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, diff, iwflb, maxu, btol);
        bdenup(x, uab, c, d, da, norbs, nclose);
        betaw=aval(h1,d,norbs);
        zerom(f, norbs);
        ffreq2(f, d, w);
        ffreq1(f, d, da, da, norbs);
        zerom(da, norbs);
        hplusf(f, da, norbs);
        while(!last){
            ++icount;
            if(icount>=maxitu) last=true;
            transf(f, g, c, norbs);
            bmakuf(ua, ub, uab, t, uold1, g, eigs, last, norbs, nclose, diff, iwflb, maxu, btol);
            bdenup(x, uab, c, d, da, norbs, nclose);
            betaw=aval(h1,d,norbs);
            zerom(f, norbs);
            ffreq2(f, d, w);
            ffreq1(f, d, da, da, norbs);
            zerom(da, norbs);
            hplusf(f, da, norbs);
        }
        cmptim=second(1)-cmptim;
        if(debug){
            std::printf("\n CONVERGED IN%4d ITERATIONS IN%10.2f SECONDS\n", icount, cmptim);
            std::printf(" MAXIMUM UAB ELEMENT =%15.5f,  MAXIMUM DIFFERENCE =%15.5f\n\n", maxu, diff);
        }
        for(int ic=1;ic<=3;++ic){
            hmuf(h1, ic, coord, nfirst, nlast, nat, norbs, numat);
            betaw=aval(h1,d,norbs);
            allbet[ic][ia][ib]=betaw;
            if(debug) std::printf("      BETA(%c,%c,%c) = %15.5f\n", alab[ic][0], alab[ia][0], alab[ib][0], betaw);
            if(id==1 && ic==1) bavx+=3.0*betaw;
            else if(id==2 && ic==2) bavx+=2.0*betaw;
            else if(id==3 && ic==3) bavx+=2.0*betaw;
            else if((id==4 || id==6) && ic==1) bavx+=betaw;
            if(id==4 && ic==2) bavy+=3.0*betaw;
            else if(id==2 && ic==1) bavy+=2.0*betaw;
            else if(id==5 && ic==3) bavy+=2.0*betaw;
            else if((id==1 || id==6) && ic==2) bavy+=betaw;
            if(id==6 && ic==3) bavz+=3.0*betaw;
            else if(id==3 && ic==1) bavz+=2.0*betaw;
            else if(id==5 && ic==2) bavz+=2.0*betaw;
            else if((id==4 || id==1) && ic==3) bavz+=betaw;
        }
        epsab(h1, eigs, g, ga, gb, ua, ub, uab, da, norbs, nclose, iwflb);
        dawrit(uab, maxsq, iposu+id);
        dawrit(g, maxsq, iposg+id);
        dawrit(h1, maxsq, ipose+id);
        dawrit(da, maxsq, iposum+id);
    }
    std::printf("                           COMPONENTS OF BETA\n");
    bavx/=5.0; bavy/=5.0; bavz/=5.0;
    double bvec=std::sqrt(bavx*bavx+bavy*bavy+bavz*bavz);
    std::printf("           *XX         *XY         *YY         *XZ         *YZ         *ZZ\n");
    if(bvec<0.999e5){
        std::printf("   *=X %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n", allbet[1][1][1], allbet[1][1][2], allbet[1][2][2], allbet[1][1][3], allbet[1][2][3], allbet[1][3][3]);
        std::printf("   *=Y %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n", allbet[2][1][1], allbet[2][1][2], allbet[2][2][2], allbet[2][1][3], allbet[2][2][3], allbet[2][3][3]);
        std::printf("   *=Z %12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n", allbet[3][1][1], allbet[3][1][2], allbet[3][2][2], allbet[3][1][3], allbet[3][2][3], allbet[3][3][3]);
    } else {
        std::printf("   *=X %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n", allbet[1][1][1], allbet[1][1][2], allbet[1][2][2], allbet[1][1][3], allbet[1][2][3], allbet[1][3][3]);
        std::printf("   *=Y %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n", allbet[2][1][1], allbet[2][1][2], allbet[2][2][2], allbet[2][1][3], allbet[2][2][3], allbet[2][3][3]);
        std::printf("   *=Z %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n", allbet[3][1][1], allbet[3][1][2], allbet[3][2][2], allbet[3][1][3], allbet[3][2][3], allbet[3][3][3]);
    }
    std::printf("\n");
    if(bvec<1e10){
        std::printf(" AVERAGE BETA X (SHG) VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n", omega, bavx, bavx*au_to_esu);
        std::printf(" AVERAGE BETA Y (SHG) VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n", omega, bavy, bavy*au_to_esu);
        std::printf(" AVERAGE BETA Z (SHG) VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n\n", omega, bavz, bavz*au_to_esu);
        std::printf(" AVERAGE BETA   (SHG) VALUE AT%10.5f EV = %15.4f a.u. = %14.6e ESU\n", omega, bvec, bvec*au_to_esu);
    } else {
        std::printf(" AVERAGE BETA X (SHG) VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n", omega, bavx, bavx*au_to_esu);
        std::printf(" AVERAGE BETA Y (SHG) VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n", omega, bavy, bavy*au_to_esu);
        std::printf(" AVERAGE BETA Z (SHG) VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n\n", omega, bavz, bavz*au_to_esu);
        std::printf(" AVERAGE BETA   (SHG) VALUE AT%10.5f EV = %13.6e a.u. = %14.6e ESU\n", omega, bvec, bvec*au_to_esu);
    }
    std::printf("   (1 a.u. = 8.639418X10-33 esu)\n\n");
}

// polar(): TDHF polarizability driver (polar.F90, 4372 lines) -- retained stub.
// TODO: full translation; subroutines bdenin/bdenup/bmakuf/epsab/fhpatn/nonor/
// nonope/nonbet/ngoke/ngidri/ngefis/ngamtg/beopor/betaf are already ported below.
void polar() {
    std::fprintf(stderr, "[polar] retained stub: TDHF polarizability driver not yet ported; "
                         "subroutine set (bdenin..betaf) is available.\n");
}

// ---- ffreq1 / ffreq2 (polar.F90:2471-2644) ----
void ffreq1(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& ptot,
            std::vector<std::vector<double>>& pa, std::vector<std::vector<double>>& pb, int ndim) {
    // One-centre contributions to the Fock matrix (TDHF, polar.F90:2471).
    for (int ii = 1; ii <= numat; ++ii) {
        int ia = nfirst[ii];
        int ib = std::min(nlast[ii], ia + 3);
        int ni = nat[ii];
        double ptpop = 0.0, papop = 0.0;
        int idx = ib - ia + 2;
        if (idx == 2) continue;            // s-only atom: handled elsewhere
        if (idx >= 4) {                    // sp + : accumulate s population
            ptpop = ptot[ib][ib] + ptot[ib-1][ib-1] + ptot[ib-2][ib-2];
            papop = pa[ib][ib] + pa[ib-1][ib-1] + pa[ib-2][ib-2];
        }
        // F(S,S)
        f[ia][ia] += pb[ia][ia]*parameters_C::gss[ni] + ptpop*parameters_C::gsp[ni]
                   - papop*parameters_C::hsp[ni];
        if (ni >= 3) {
            int iplus = ia + 1;
            for (int j = iplus; j <= ib; ++j) {
                // F(P,P)
                f[j][j] += ptot[ia][ia]*parameters_C::gsp[ni] - pa[ia][ia]*parameters_C::hsp[ni]
                         + pb[j][j]*parameters_C::gpp[ni]
                         + (ptpop - ptot[j][j])*parameters_C::gp2[ni]
                         - 0.5*(papop - pa[j][j])*(parameters_C::gpp[ni] - parameters_C::gp2[ni]);
                // F(S,P)
                f[ia][j] += 2.0*ptot[ia][j]*parameters_C::hsp[ni]
                          - pa[ia][j]*(parameters_C::hsp[ni] + parameters_C::gsp[ni]);
                f[j][ia] += 2.0*ptot[j][ia]*parameters_C::hsp[ni]
                          - pa[j][ia]*(parameters_C::hsp[ni] + parameters_C::gsp[ni]);
            }
            // F(P,P*)
            int iminus = ib - 1;
            for (int j = iplus; j <= iminus; ++j) {
                int icc = j + 1;
                for (int k = icc; k <= ib; ++k) {
                    f[j][k] += ptot[j][k]*(parameters_C::gpp[ni] - parameters_C::gp2[ni])
                             - 0.5*pa[j][k]*(parameters_C::gpp[ni] + parameters_C::gp2[ni]);
                    f[k][j] += ptot[k][j]*(parameters_C::gpp[ni] - parameters_C::gp2[ni])
                             - 0.5*pa[k][j]*(parameters_C::gpp[ni] + parameters_C::gp2[ni]);
                }
            }
        }
    }
}

void ffreq2(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& ptot,
            std::vector<double>& w) {
    // TDHF two-centre two-electron repulsion part of the Fock matrix (polar.F90:2561).
    int kk = 0;
    for (int ii = 1; ii <= numat; ++ii) {
        int iim1 = ii - 1;
        int ia = nfirst[ii];
        int ib = nlast[ii];
        for (int jj = 1; jj <= iim1; ++jj) {
            int ja = nfirst[jj];
            int jb = nlast[jj];
            for (int i = ia; i <= ib; ++i) {
                for (int j = ia; j <= i; ++j) {
                    double fij = (i == j) ? 0.5 : 1.0;
                    for (int k = ja; k <= jb; ++k) {
                        for (int l = ja; l <= k; ++l) {
                            double fkl = (k == l) ? 0.5 : 1.0;
                            kk += 1;
                            double a = w[kk];
                            double aint = a*fkl*fij;
                            double pkl = ptot[k][l] + ptot[l][k];
                            f[i][j] += aint*pkl;
                            f[j][i] += aint*pkl;
                            double pij = ptot[i][j] + ptot[j][i];
                            f[k][l] += aint*pij;
                            f[l][k] += aint*pij;
                            aint *= 0.5;
                            f[i][l] -= aint*ptot[j][k];
                            f[l][i] -= aint*ptot[k][j];
                            f[k][j] -= aint*ptot[l][i];
                            f[j][k] -= aint*ptot[i][l];
                            f[i][k] -= aint*ptot[j][l];
                            f[k][i] -= aint*ptot[l][j];
                            f[j][l] -= aint*ptot[i][k];
                            f[l][j] -= aint*ptot[k][i];
                        }
                    }
                }
            }
        }
        kk += ((ib - ia + 1) * (ib - ia + 2) / 2) * ((ib - ia + 1) * (ib - ia + 2) / 2);
    }
}
