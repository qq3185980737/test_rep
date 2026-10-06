// mndod.cpp — MNDO-d (d-orbital) Hamiltonian routines (from mndod.F90).
#include "mndod.h"
#include "parameters_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "funcon_C.h"
#include <cmath>
#include <vector>
#include <cstring>
namespace mndod_C {
extern std::vector<std::vector<int>> indexd, indx, ind2;
extern std::vector<int> isym;
extern std::vector<std::vector<std::vector<double>>> ch;
extern std::vector<std::vector<int>> indpp, inddp, inddd;
}
namespace funcon_C { extern double ev; }
using namespace parameters_C;
using common_arrays_C::po;
using funcon_C::ev;
using funcon_C::a0;
using molkst_C::method_pm7;
void inighd(int ni);

void aijm(int ni);
void ddpo(int ni);
void scprm(int ni,double&,double&,double&,double&,double&,double&,double&,double&,double&,double&,double&,double&);
void eiscor(double,double,double,double,double,int);

void inid(){
    for(int ni=1;ni<=107;++ni){
        if(!dorbs[ni]) continue;
        aijm(ni);
        if(zdn[ni] > 1e-4) inighd(ni);
        ddpo(ni);
    }
    for(int i=1;i<=106;++i){
        if(natorb[i]<6 || main_group[i]){
            if(am[i]<1e-4) am[i]=1.0;
            po[1][i]=0.5/am[i];
            if(ad[i]>1e-5) po[2][i]=0.5/ad[i];
            if(aq[i]>1e-5) po[3][i]=0.5/aq[i];
            po[7][i]=po[1][i];
            ddp[2][i]=dd[i];
            ddp[3][i]=qq[i]*std::sqrt(2.0);
        }
        po[9][i]=po[1][i];
        if(pocord[i]>1e-5) po[9][i]=pocord[i];
    }
    po[2][1]=0.0; po[3][1]=0.0;
}

namespace mndod_C { extern std::vector<std::vector<double>> repd; }
void inighd(int ni){
    const double s3=1.7320508, s5=2.23606797, s15=3.87298334;
    if(!dorbs[ni]) return;
    double r066=0,r266=0,r466=0,r016=0,r244=0,r036=0,r236=0,r155=0,r355=0,r125=0,r234=0,r246=0;
    scprm(ni,r066,r266,r466,r016,r244,r036,r236,r155,r355,r125,r234,r246);
    if(f0sd[ni]>0.001) r016=f0sd[ni];
    if(g2sd[ni]>0.001) r244=g2sd[ni];
    eiscor(r016,r066,r244,r266,r466,ni);
    std::vector<double> R(53,0.0);
    R[1]=r016;
    R[2]=2.0/(3.0*s5)*r125;
    R[3]=1.0/s15*r125;
    R[4]=2.0/(5.0*s5)*r234;
    R[5]=r036+4.0/35.0*r236;
    R[6]=r036+2.0/35.0*r236;
    R[7]=r036-4.0/35.0*r236;
    R[8]=-1.0/(3.0*s5)*r125;
    R[9]=std::sqrt(3.0/125.0)*r234;
    R[10]=s3/35.0*r236;
    R[11]=3.0/35.0*r236;
    R[12]=-1.0/(5.0*s5)*r234;
    R[13]=r036-2.0/35.0*r236;
    R[14]=-2.0*s3/35.0*r236;
    R[15]=-R[3]; R[16]=-R[11]; R[17]=-R[9]; R[18]=-R[14];
    R[19]=1.0/5.0*r244;
    R[20]=2.0/(7.0*s5)*r246;
    R[21]=R[20]/2.0; R[22]=-R[20];
    R[23]=4.0/15.0*r155+27.0/245.0*r355;
    R[24]=2.0*s3/15.0*r155-9.0*s3/245.0*r355;
    R[25]=1.0/15.0*r155+18.0/245.0*r355;
    R[26]=(-s3/15.0*r155)+12.0*s3/245.0*r355;
    R[27]=(-s3/15.0*r155)-3.0*s3/245.0*r355;
    R[28]=-R[27];
    R[29]=r066+4.0/49.0*r266+4.0/49.0*r466;
    R[30]=r066+2.0/49.0*r266-24.0/441.0*r466;
    R[31]=r066-4.0/49.0*r266+6.0/441.0*r466;
    R[32]=std::sqrt(3.0/245.0)*r246;
    R[33]=1.0/5.0*r155+24.0/245.0*r355;
    R[34]=1.0/5.0*r155-6.0/245.0*r355;
    R[35]=3.0/49.0*r355;
    R[36]=1.0/49.0*r266+30.0/441.0*r466;
    R[37]=s3/49.0*r266-5.0*s3/441.0*r466;
    R[38]=r066-2.0/49.0*r266-4.0/441.0*r466;
    R[39]=(-2.0*s3/49.0*r266)+10.0*s3/441.0*r466;
    R[40]=-R[32]; R[41]=-R[34]; R[42]=-R[35]; R[43]=-R[37];
    R[44]=3.0/49.0*r266+20.0/441.0*r466;
    R[45]=-R[39];
    R[46]=1.0/5.0*r155-3.0/35.0*r355;
    R[47]=-R[46];
    R[48]=4.0/49.0*r266+15.0/441.0*r466;
    R[49]=3.0/49.0*r266-5.0/147.0*r466;
    R[50]=-R[49];
    R[51]=r066+4.0/49.0*r266-34.0/441.0*r466;
    R[52]=35.0/441.0*r466;
    for(int I=1;I<=52;++I) mndod_C::repd[I][ni]=R[I];
    f0dd[ni]=r066; f2dd[ni]=r266; f4dd[ni]=r466;
    f0sd[ni]=r016; g2sd[ni]=r244;
    f0pd[ni]=r036; f2pd[ni]=r236; g1pd[ni]=r155; g3pd[ni]=r355;
}

double poij(int l, double d, double fg){
    const double epsil=1e-8, g1=0.382, g2=0.618;
    if(l==0) return 0.5*ev/fg;
    double dsq=d*d, ev4=ev*0.25, ev8=ev/8.0;
    double a1=0.1, a2=5.0, f1=0,f2=0;
    auto bracket=[&](auto body){
        for(int i=0;i<100;++i){
            double delta=a2-a1; if(delta<epsil) break;
            double y1=a1+delta*g1, y2=a1+delta*g2;
            body(y1,y2);
            if(f1<f2) a2=y2; else a1=y1;
        }
    };
    if(l==1){
        bracket([&](double y1,double y2){
            f1=std::pow(ev4*(1.0/y1-1.0/std::sqrt(y1*y1+dsq))-fg,2);
            f2=std::pow(ev4*(1.0/y2-1.0/std::sqrt(y2*y2+dsq))-fg,2);
        });
    } else if(l==2){
        bracket([&](double y1,double y2){
            f1=std::pow(ev8*(1.0/y1-2.0/std::sqrt(y1*y1+dsq*0.5)+1.0/std::sqrt(y1*y1+dsq))-fg,2);
            f2=std::pow(ev8*(1.0/y2-2.0/std::sqrt(y2*y2+dsq*0.5)+1.0/std::sqrt(y2*y2+dsq))-fg,2);
        });
    } else {
        bracket([&](double,double){});
    }
    return (f1>=f2)?a2:a1;
}

extern void to_point(double rij, double& point, double& cnst);

void reppd(int ni, int nj, double rij, std::vector<double>& ri, double& gab){
    const double td=2.0, pp=0.5;
    int nri[23]={0,1,-1,1,1,-1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,1,1,1,1,1};
    double ev1=ev/2, ev2=ev1/2, ev3=ev2/2, ev4=ev3/2;
    ri.assign(23,0.0);
    double r=rij/a0;
    bool si=natorb[ni]>=3, sj=natorb[nj]>=3;
    std::vector<double> arg(73,0.0), sqr(73,0.0);
    double aee=po[9][ni]+po[9][nj]; aee*=aee;
    gab=ev/std::sqrt(r*r+aee);
    aee=pp/am[ni]+pp/am[nj]; aee*=aee;
    if(!si && !sj){ ri[1]=ev/std::sqrt(r*r+aee); }
    else if(si && !sj){
        double da=dd[ni], qa=qq[ni]*td;
        double ade=pp/ad[ni]+pp/am[nj]; ade*=ade;
        double aqe=pp/aq[ni]+pp/am[nj]; aqe*=aqe;
        double rsq=r*r, xxx;
        arg[1]=rsq+aee;
        xxx=r+da; arg[2]=xxx*xxx+ade;
        xxx=r-da; arg[3]=xxx*xxx+ade;
        xxx=r+qa; arg[4]=xxx*xxx+aqe;
        xxx=r-qa; arg[5]=xxx*xxx+aqe;
        arg[6]=rsq+aqe; arg[7]=arg[6]+qa*qa;
        for(int i=1;i<=7;++i) sqr[i]=std::sqrt(arg[i]);
        double ee=ev/sqr[1];
        ri[1]=ee; ri[2]=ev1/sqr[2]-ev1/sqr[3];
        ri[3]=ee+ev2/sqr[4]+ev2/sqr[5]-ev1/sqr[6];
        ri[4]=ee+ev1/sqr[7]-ev1/sqr[6];
    }
    else if(!si && sj){
        double db=dd[nj], qb=qq[nj]*td;
        double aed=pp/am[ni]+pp/ad[nj]; aed*=aed;
        double aeq=pp/am[ni]+pp/aq[nj]; aeq*=aeq;
        double rsq=r*r, xxx;
        arg[1]=rsq+aee;
        xxx=r-db; arg[2]=xxx*xxx+aed;
        xxx=r+db; arg[3]=xxx*xxx+aed;
        xxx=r-qb; arg[4]=xxx*xxx+aeq;
        xxx=r+qb; arg[5]=xxx*xxx+aeq;
        arg[6]=rsq+aeq; arg[7]=arg[6]+qb*qb;
        for(int i=1;i<=7;++i) sqr[i]=std::sqrt(arg[i]);
        double ee=ev/sqr[1];
        ri[1]=ee; ri[5]=ev1/sqr[2]-ev1/sqr[3];
        ri[11]=ee+ev2/sqr[4]+ev2/sqr[5]-ev1/sqr[6];
        ri[12]=ee+ev1/sqr[7]-ev1/sqr[6];
    }
    else {
        double da=dd[ni], db=dd[nj], qa=qq[ni]*td, qb=qq[nj]*td;
        double ade=pp/ad[ni]+pp/am[nj]; ade*=ade;
        double aqe=pp/aq[ni]+pp/am[nj]; aqe*=aqe;
        double aed=pp/am[ni]+pp/ad[nj]; aed*=aed;
        double aeq=pp/am[ni]+pp/aq[nj]; aeq*=aeq;
        double axx=pp/ad[ni]+pp/ad[nj]; axx*=axx;
        double adq=pp/ad[ni]+pp/aq[nj]; adq*=adq;
        double aqd=pp/aq[ni]+pp/ad[nj]; aqd*=aqd;
        double aqq=pp/aq[ni]+pp/aq[nj]; aqq*=aqq;
        double rsq=r*r, xxx;
        arg[1]=rsq+aee;
        xxx=r+da; arg[2]=xxx*xxx+ade; xxx=r-da; arg[3]=xxx*xxx+ade;
        xxx=r-qa; arg[4]=xxx*xxx+aqe; xxx=r+qa; arg[5]=xxx*xxx+aqe;
        arg[6]=rsq+aqe; arg[7]=arg[6]+qa*qa;
        xxx=r-db; arg[8]=xxx*xxx+aed; xxx=r+db; arg[9]=xxx*xxx+aed;
        xxx=r-qb; arg[10]=xxx*xxx+aeq; xxx=r+qb; arg[11]=xxx*xxx+aeq;
        arg[12]=rsq+aeq; arg[13]=arg[12]+qb*qb;
        xxx=da-db; arg[14]=rsq+axx+xxx*xxx;
        xxx=da+db; arg[15]=rsq+axx+xxx*xxx;
        xxx=r+da-db; arg[16]=xxx*xxx+axx; xxx=r-da+db; arg[17]=xxx*xxx+axx;
        xxx=r-da-db; arg[18]=xxx*xxx+axx; xxx=r+da+db; arg[19]=xxx*xxx+axx;
        xxx=r+da; arg[20]=xxx*xxx+adq; arg[21]=arg[20]+qb*qb;
        xxx=r-da; arg[22]=xxx*xxx+adq; arg[23]=arg[22]+qb*qb;
        xxx=r-db; arg[24]=xxx*xxx+aqd; arg[25]=arg[24]+qa*qa;
        xxx=r+db; arg[26]=xxx*xxx+aqd; arg[27]=arg[26]+qa*qa;
        xxx=r+da-qb; arg[28]=xxx*xxx+adq; xxx=r-da-qb; arg[29]=xxx*xxx+adq;
        xxx=r+da+qb; arg[30]=xxx*xxx+adq; xxx=r-da+qb; arg[31]=xxx*xxx+adq;
        xxx=r+qa-db; arg[32]=xxx*xxx+aqd; xxx=r+qa+db; arg[33]=xxx*xxx+aqd;
        xxx=r-qa-db; arg[34]=xxx*xxx+aqd; xxx=r-qa+db; arg[35]=xxx*xxx+aqd;
        arg[36]=rsq+aqq;
        xxx=qa-qb; arg[37]=arg[36]+xxx*xxx;
        xxx=qa+qb; arg[38]=arg[36]+xxx*xxx;
        arg[39]=arg[36]+qa*qa; arg[40]=arg[36]+qb*qb; arg[41]=arg[39]+qb*qb;
        xxx=r-qb; arg[42]=xxx*xxx+aqq; arg[43]=arg[42]+qa*qa;
        xxx=r+qb; arg[44]=xxx*xxx+aqq; arg[45]=arg[44]+qa*qa;
        xxx=r+qa; arg[46]=xxx*xxx+aqq; arg[47]=arg[46]+qb*qb;
        xxx=r-qa; arg[48]=xxx*xxx+aqq; arg[49]=arg[48]+qb*qb;
        xxx=r+qa-qb; arg[50]=xxx*xxx+aqq; xxx=r+qa+qb; arg[51]=xxx*xxx+aqq;
        xxx=r-qa-qb; arg[52]=xxx*xxx+aqq; xxx=r-qa+qb; arg[53]=xxx*xxx+aqq;
        qa=qq[ni]; qb=qq[nj];
        double yyy,zzz,www;
        xxx=da-qb; xxx*=xxx; yyy=r-qb; yyy*=yyy; zzz=da+qb; zzz*=zzz; www=r+qb; www*=www;
        arg[54]=xxx+yyy+adq; arg[55]=xxx+www+adq; arg[56]=zzz+yyy+adq; arg[57]=zzz+www+adq;
        xxx=qa-db; xxx*=xxx; yyy=qa+db; yyy*=yyy; zzz=r+qa; zzz*=zzz; www=r-qa; www*=www;
        arg[58]=zzz+xxx+aqd; arg[59]=www+xxx+aqd; arg[60]=zzz+yyy+aqd; arg[61]=www+yyy+aqd;
        xxx=qa-qb; xxx*=xxx;
        arg[62]=arg[36]+td*xxx;
        yyy=qa+qb; yyy*=yyy; arg[63]=arg[36]+td*yyy;
        arg[64]=arg[36]+td*(qa*qa+qb*qb);
        zzz=r+qa-qb; zzz*=zzz; arg[65]=zzz+xxx+aqq; arg[66]=zzz+yyy+aqq;
        zzz=r+qa+qb; zzz*=zzz; arg[67]=zzz+xxx+aqq; arg[68]=zzz+yyy+aqq;
        zzz=r-qa-qb; zzz*=zzz; arg[69]=zzz+xxx+aqq; arg[70]=zzz+yyy+aqq;
        zzz=r-qa+qb; zzz*=zzz; arg[71]=zzz+xxx+aqq; arg[72]=zzz+yyy+aqq;
        for(int i=1;i<=72;++i) sqr[i]=std::sqrt(arg[i]);
        double ee=ev/sqr[1];
        double dze=-ev1/sqr[2]+ev1/sqr[3];
        double qzze=ev2/sqr[4]+ev2/sqr[5]-ev1/sqr[6];
        double qxxe=ev1/sqr[7]-ev1/sqr[6];
        double edz=-ev1/sqr[8]+ev1/sqr[9];
        double eqzz=ev2/sqr[10]+ev2/sqr[11]-ev1/sqr[12];
        double eqxx=ev1/sqr[13]-ev1/sqr[12];
        double dxdx=ev1/sqr[14]-ev1/sqr[15];
        double dzdz=ev2/sqr[16]+ev2/sqr[17]-ev2/sqr[18]-ev2/sqr[19];
        double dzqxx=ev2/sqr[20]-ev2/sqr[21]-ev2/sqr[22]+ev2/sqr[23];
        double qxxdz=ev2/sqr[24]-ev2/sqr[25]-ev2/sqr[26]+ev2/sqr[27];
        double dzqzz=-ev3/sqr[28]+ev3/sqr[29]-ev3/sqr[30]+ev3/sqr[31]-ev2/sqr[22]+ev2/sqr[20];
        double qzzdz=-ev3/sqr[32]+ev3/sqr[33]-ev3/sqr[34]+ev3/sqr[35]+ev2/sqr[24]-ev2/sqr[26];
        double qxxqxx=ev3/sqr[37]+ev3/sqr[38]-ev2/sqr[39]-ev2/sqr[40]+ev2/sqr[36];
        double qxxqyy=ev2/sqr[41]-ev2/sqr[39]-ev2/sqr[40]+ev2/sqr[36];
        double qxxqzz=ev3/sqr[43]+ev3/sqr[45]-ev3/sqr[42]-ev3/sqr[44]-ev2/sqr[39]+ev2/sqr[36];
        double qzzqxx=ev3/sqr[47]+ev3/sqr[49]-ev3/sqr[46]-ev3/sqr[48]-ev2/sqr[40]+ev2/sqr[36];
        double qzzqzz=ev4/sqr[50]+ev4/sqr[51]+ev4/sqr[52]+ev4/sqr[53]-ev3/sqr[48]-ev3/sqr[46]-ev3/sqr[42]-ev3/sqr[44]+ev2/sqr[36];
        double dxqxz=-ev2/sqr[54]+ev2/sqr[55]+ev2/sqr[56]-ev2/sqr[57];
        double qxzdx=-ev2/sqr[58]+ev2/sqr[59]+ev2/sqr[60]-ev2/sqr[61];
        double qxzqxz=ev3/sqr[65]-ev3/sqr[67]-ev3/sqr[69]+ev3/sqr[71]-ev3/sqr[66]+ev3/sqr[68]+ev3/sqr[70]-ev3/sqr[72];
        ri[1]=ee; ri[2]=-dze; ri[3]=ee+qzze; ri[4]=ee+qxxe; ri[5]=-edz;
        ri[6]=dzdz; ri[7]=dxdx; ri[8]=-edz-qzzdz; ri[9]=-edz-qxxdz; ri[10]=-qxzdx;
        ri[11]=ee+eqzz; ri[12]=ee+eqxx; ri[13]=-dze-dzqzz; ri[14]=-dze-dzqxx; ri[15]=-dxqxz;
        ri[16]=ee+eqzz+qzze+qzzqzz; ri[17]=ee+eqzz+qxxe+qxxqzz;
        ri[18]=ee+eqxx+qzze+qzzqxx; ri[19]=ee+eqxx+qxxe+qxxqxx;
        ri[20]=qxzqxz; ri[21]=ee+eqxx+qxxe+qxxqyy;
        ri[22]=pp*(qxxqxx-qxxqyy);
    }
    if(method_pm7){
        double point=0,cnst=0; to_point(rij,point,cnst);
        for(int idx=1;idx<=22;++idx){
            bool blend=(idx==1||idx==3||idx==4||idx==11||idx==12||idx==16||idx==17||idx==18||idx==19||idx==21);
            ri[idx]=ri[idx]*cnst+(blend?(1.0-cnst)*point:0.0);
        }
        gab=gab*cnst+(1.0-cnst)*point;
    }
    for(int idx=1;idx<=22;++idx) ri[idx]*=nri[idx];
}

namespace mndod_C {
extern std::vector<std::vector<int>> indexd, ind2;
extern std::vector<int> isym;
}
extern double rijkl(int ni,int nj,int ij,int kl,int li,int lj,int lk,int ll,int ic,double r);
using mndod_C::indexd; using mndod_C::ind2; using mndod_C::isym;

void reppd2(int ni,int nj,double r,const std::vector<double>& ri,
            std::vector<double>& rep,double core[11][3]){
    int ipos[35]={0,1,5,11,12,12,2,6,13,14,14,3,8,16,18,18,7,15,10,20,4,9,17,19,21,7,15,10,20,22,4,9,17,21,19};
    int lorb[10]={0,0,1,1,1,2,2,2,2,2};
    rep.assign(492,0.0);
    for(int idx=1;idx<=34;++idx) rep[idx]=ri[ipos[idx]];
    if(dorbs[ni] || dorbs[nj]){
        int lasti = dorbs[ni]?9:(ni<3?1:4);
        int lastk = dorbs[nj]?9:(nj<3?1:4);
        int ij=0,kl=0;
        for(int i=1;i<=lasti;++i){
            int li=lorb[i];
            for(int j=1;j<=i;++j){
                bool coul=(i==j); int lj=lorb[j];
                ij=indexd[i][j];
                for(int k=1;k<=lastk;++k){
                    int lk=lorb[k];
                    for(int l=1;l<=k;++l){
                        bool coulomb=(coul && k==l); int ll=lorb[l];
                        kl=indexd[k][l];
                        int numb=ind2[ij][kl];
                        if(numb<=34) continue;
                        int nold=isym[numb];
                        if(nold>=35) rep[numb]=rep[nold];
                        else if(nold<=-35) rep[numb]=-rep[-nold];
                        else {
                            rep[numb]=rijkl(ni,nj,ij,kl,li,lj,lk,ll,0,r)*ev;
                            if(method_pm7){
                                double point=0,cnst=0; to_point(r*a0,point,cnst);
                                rep[numb]=rep[numb]*cnst+(coulomb?(1.0-cnst)*point:0.0);
                            }
                        }
                    }
                }
            }
        }
        for(int q=5;q<=10;++q){ core[q][1]=0.0; core[q][2]=0.0; }
        if(dorbs[nj]){
            kl=indexd[5][1]; core[5][2]=-rijkl(ni,nj,ij,kl,0,0,2,0,1,r)*ev*tore[ni];
            kl=indexd[5][2]; core[6][2]=-rijkl(ni,nj,ij,kl,0,0,2,1,1,r)*ev*tore[ni];
            kl=indexd[5][5]; core[7][2]=-rijkl(ni,nj,ij,kl,0,0,2,2,1,r)*ev*tore[ni];
            kl=indexd[6][3]; core[8][2]=-rijkl(ni,nj,ij,kl,0,0,2,1,1,r)*ev*tore[ni];
            kl=indexd[6][6]; core[9][2]=-rijkl(ni,nj,ij,kl,0,0,2,2,1,r)*ev*tore[ni];
            kl=indexd[8][8]; core[10][2]=-rijkl(ni,nj,ij,kl,0,0,2,2,1,r)*ev*tore[ni];
        }
        if(dorbs[ni]){
            kl=indexd[5][1]; core[5][1]=-rijkl(ni,nj,kl,ij,2,0,0,0,2,r)*ev*tore[nj];
            kl=indexd[5][2]; core[6][1]=-rijkl(ni,nj,kl,ij,2,1,0,0,2,r)*ev*tore[nj];
            kl=indexd[5][5]; core[7][1]=-rijkl(ni,nj,kl,ij,2,2,0,0,2,r)*ev*tore[nj];
            kl=indexd[6][3]; core[8][1]=-rijkl(ni,nj,kl,ij,2,1,0,0,2,r)*ev*tore[nj];
            kl=indexd[6][6]; core[9][1]=-rijkl(ni,nj,kl,ij,2,2,0,0,2,r)*ev*tore[nj];
            kl=indexd[8][8]; core[10][1]=-rijkl(ni,nj,kl,ij,2,2,0,0,2,r)*ev*tore[nj];
        }
    }
}

namespace mndod_C {
extern std::vector<std::vector<int>> indx;
extern std::vector<std::vector<std::vector<double>>> ch;
}
extern double charg(double r,int l1,int l2,int mm,double da,double db,double add);
using mndod_C::indx; using mndod_C::ch;

double rijkl(int ni,int nj,int ij,int kl,int li,int lj,int lk,int ll,int ic,double r){
    int l1min=std::abs(li-lj), l1max=li+lj;
    int lij=indx[li+1][lj+1];
    int l2min=std::abs(lk-ll), l2max=lk+ll;
    int lkl=indx[lk+1][ll+1];
    l1max=std::min(l1max,2); l1min=std::min(l1min,2);
    l2max=std::min(l2max,2); l2min=std::min(l2min,2);
    double sum=0.0, pij=0, dij=0, pkl=0, dkl=0;
    for(int l1=l1min;l1<=l1max;++l1){
        if(l1==0){
            if(lij==1){ pij=po[1][ni]; if(ic==1) pij=po[9][ni]; }
            else if(lij==3){ pij=po[7][ni]; }
            else if(lij==6){ pij=po[8][ni]; }
        } else {
            dij=ddp[lij][ni]; pij=po[lij][ni];
        }
        for(int l2=l2min;l2<=l2max;++l2){
            if(l2==0){
                if(lkl==1){ pkl=po[1][nj]; if(ic==2) pkl=po[9][nj]; }
                else if(lkl==3){ pkl=po[7][nj]; }
                else if(lkl==6){ pkl=po[8][nj]; }
            } else {
                dkl=ddp[lkl][nj]; pkl=po[lkl][nj];
            }
            double add=(pij+pkl)*(pij+pkl);
            int lmin=std::min(l1,l2);
            double s1=0.0;
            for(int m=-lmin;m<=lmin;++m){
                // F90 mndod.F90:1079 ccc = ch(ij,l1,m)*ch(kl,l2,m); ch dims (45,0:2,-2:2)
                // C++ ch third dim is m+2 (0..4 for m=-2..2).
                double ccc=ch[ij][l1][m+2]*ch[kl][l2][m+2];
                if(ccc==0.0) continue;
                int mm=std::abs(m);
                double cg=charg(r,l1,l2,mm,dij,dkl,add);
                s1+=cg*ccc;
            }
            sum+=s1;
        }
    }
    return sum;
}

namespace mndod_C {
extern thread_local double sp[4][4], sd[6][6]; extern thread_local double pp[7][4][4], dp[16][6][4], d_d[16][6][6];

extern thread_local double cored[11][3];
}
extern void rotmat(int nj,int ni,const std::vector<double>&,const std::vector<double>&,double&);
extern void spcore(int ni,int nj,double r,double cored[11][3]);
extern void tx(int ii,int kk,const std::vector<double>& rep,char* logv,double* v);
extern void w2mat(const std::vector<double>& ww,std::vector<double>& w,int& kr,int iw,int jw);
extern void elenuc(int,int,int,int,std::vector<double>&);
extern void ccrep(int ni,int nj,double& r,double gab,double& enuclr);
using mndod_C::sp; using mndod_C::sd; using mndod_C::pp; using mndod_C::dp; using mndod_C::d_d;
using mndod_C::cored;
using mndod_C::inddd;

void rotatd(int ni,int nj,const std::vector<double>& ci,const std::vector<double>& cj,
            std::vector<double>& w,int& kr,double& enuc){
    int met[46]={0,1,2,3,2,3,3,2,3,3,3,4,5,5,5,6,4,5,5,5,6,6,4,5,5,5,6,6,6,4,5,5,5,6,6,6,6,4,5,5,5,6,6,6,6,6}; // F90 data met ends with 5*6
    double r=0,gab=0;
    // Flat 1-D scratch matrices (row stride 46). rotatd is called ~650x per
    // dcart; per-call heap allocation of the 46x46 matrices was a measurable
    // hidden cost, and flat access removes the vector<vector> double
    // indirection in the hot four-level loop below.
    std::vector<double> ri(23,0), rep(492,0), ww(2026,0), en(172,0);
    char logv[46 * 46];
    double v[46 * 46];
    std::memset(logv, 0, sizeof(logv));
    std::memset(v, 0, sizeof(v));
    rotmat(nj,ni,ci,cj,r);
    reppd(ni,nj,r,ri,gab);
    double rij=r; r=r/a0;
    spcore(ni,nj,r,cored);
    double const_,point;
    if(method_pm7){ to_point(rij,point,const_); } else { const_=1.0; point=0.0; }
    reppd2(ni,nj,r,ri,rep,cored);
    point=-(ev/r)*tore[nj];
    for(int q=1;q<=10;++q){
        bool blend=(q==1||q==3||q==4||q==7||q==9||q==10);
        cored[q][1]=cored[q][1]*const_+(blend?(1.0-const_)*point:0.0);
    }
    point=-(ev/r)*tore[ni];
    for(int q=1;q<=10;++q){
        bool blend=(q==1||q==3||q==4||q==7||q==9||q==10);
        cored[q][2]=cored[q][2]*const_+(blend?(1.0-const_)*point:0.0);
    }
    int ii=natorb[ni], kk=natorb[nj];
    if (ni==1 && nj==2) {
        fprintf(stderr, "[RT] ni=%d nj=%d ci1=%+.6f r=%.4f ri1=%.4f rep1=%.4f rep3=%.4f\n",
                ni, nj, ci[1], r, ri[1], rep[1], rep[3]);
        fflush(stderr);
    }
    if(ii*kk>0){
        int limij=indx[ii][ii], limkl=indx[kk][kk];
        auto indw=[&](int i,int j,int kl_){ return (indx[i][j]-1)*limkl + kl_; };
        tx(ii,kk,rep,logv,v);
            for(int i1=1;i1<=ii;++i1) for(int j1=1;j1<=i1;++j1){
            int ij=indexd[i1][j1]; int jj=indx[i1][j1]; int mm=met[jj];
            for(int k=1;k<=kk;++k) for(int l=1;l<=k;++l){
                int kl=indx[k][l];
                if(!logv[ij*46 + kl]) continue;
                double wrepp=v[ij*46 + kl];                if(mm==1){ ww[indw(1,1,kl)]=wrepp; }
                else if(mm==2){ for(int i=1;i<=3;++i) ww[indw(i+1,1,kl)]+=sp[i1-1][i]*wrepp; }
                else if(mm==3){
                    for(int i=1;i<=3;++i){
                        double cc=pp[i][i1-1][j1-1];
                        ww[indw(i+1,i+1,kl)]+=cc*wrepp;
                        int iminus=i-1; if(iminus==0) continue;
                        for(int j=1;j<=iminus;++j){ cc=pp[1+i+j][i1-1][j1-1]; ww[indw(i+1,j+1,kl)]+=cc*wrepp; }
                    }
                }
                else if(mm==4){ for(int i=1;i<=5;++i) ww[indw(i+4,1,kl)]+=sd[i1-4][i]*wrepp; }
                else if(mm==5){
                    for(int i=1;i<=5;++i) for(int j=1;j<=3;++j){
                        int ij1=3*(i-1)+j;
                        ww[indw(i+4,j+1,kl)]+=dp[ij1][i1-4][j1-1]*wrepp;
                    }
                }
                else if(mm==6){
                    for(int i=1;i<=5;++i){
                        double cc=d_d[i][i1-4][j1-4];
                        ww[indw(i+4,i+4,kl)]+=cc*wrepp;
                        int iminus=i-1; if(iminus==0) continue;
                        for(int j=1;j<=iminus;++j){ int ij1=inddd[i][j]; cc=d_d[ij1][i1-4][j1-4]; ww[indw(i+4,j+4,kl)]+=cc*wrepp; }
                    }
                }
            }
        }
    }
    if(method_pm7){
        double sum=0;
        if(iod[ni]>0){
            int k = (natorb[nj]==9)?45:(natorb[nj]==4?10:1);
            if(natorb[nj]>1){
                sum=0;
                for(int i=5;i<=9;++i){ int j=k*((i*(i+1))/2-1); sum+=ww[j+3]+ww[j+6]+ww[j+10]; }
                sum=ww[1]-sum/15.0;
                for(int i=5;i<=9;++i){ int j=k*((i*(i+1))/2-1); for(int l=2;l<=4;++l) ww[(l*(l+1))/2+j]+=sum; }
                sum=cored[1][1]-(cored[3][1]+2.0*cored[4][1])/3.0;
                cored[3][1]+=sum; cored[4][1]+=sum;
            }
            sum=0; for(int i=5;i<=9;++i) sum+=ww[k*((i*(i+1))/2-1)+1];
            sum=ww[1]-sum/5.0;
            for(int i=5;i<=9;++i) ww[k*((i*(i+1))/2-1)+1]+=sum;
            sum=cored[1][1]-(cored[7][1]+2.0*cored[9][1]+2.0*cored[10][1])/5.0;
            cored[7][1]+=sum; cored[9][1]+=sum; cored[10][1]+=sum;
        }
        if(iod[nj]>0){
            sum=0;
            if(iod[ni]>0){
                for(int i=5;i<=9;++i){ int j=45*((i*(i+1))/2-1); sum+=ww[j+15]+ww[j+21]+ww[j+28]+ww[j+36]+ww[j+45]; }
                sum=ww[1]-sum/25.0;
                for(int i=5;i<=9;++i){ int j=45*((i*(i+1))/2-1); for(int k=5;k<=9;++k) ww[(k*(k+1))/2+j]+=sum; }
            }
            if(natorb[ni]>1){
                sum=0; for(int i=2;i<=4;++i){ int j=45*((i*(i+1))/2-1); sum+=ww[j+15]+ww[j+21]+ww[j+28]+ww[j+36]+ww[j+45]; }
                sum=ww[1]-sum/15.0;
                for(int i=2;i<=4;++i){ int j=45*((i*(i+1))/2-1); for(int k=5;k<=9;++k) ww[(k*(k+1))/2+j]+=sum; }
                sum=cored[1][2]-(cored[3][2]+2.0*cored[4][2])/3.0;
                cored[3][2]+=sum; cored[4][2]+=sum;
            }
            sum=ww[1]-(ww[15]+ww[21]+ww[28]+ww[36]+ww[45])/5.0;
            for(int k=5;k<=9;++k) ww[(k*(k+1))/2]+=sum;
            sum=cored[1][2]-(cored[7][2]+2.0*cored[9][2]+2.0*cored[10][2])/5.0;
            cored[7][2]+=sum; cored[9][2]+=sum; cored[10][2]+=sum;
        }
    }
    int iw=(ii*(ii+1))/2, jw=(kk*(kk+1))/2;
    if (ni==1 && nj==2) {
        fprintf(stderr, "[RT] ci1=%+.6f ww1=%.4f ww2=%.4f ww3=%.4f v47=%.4f v1=%.4f id11=%d ind21=%d limij=%d limkl=%d pp11=%.6f pp22=%.6f c31=%.6f c41=%.6f\n",
                ci[1], ww[1], ww[2], ww[3], v[1*46+1], v[1], indexd[1][1], ind2[1][1], iw, jw,
                pp[1][1][1], pp[1][2][2], pp[1][3][3], cored[3][1], cored[4][1]);
        fflush(stderr);
    }
    { static int prt=0; if (prt<3) { fprintf(stderr, "[RT2] ni=%d nj=%d natorbni=%d\n", ni, nj, natorb[ni]); fflush(stderr); prt++; } }
    w2mat(ww,w,kr,iw,jw);
    double sum=0;
    for(int i=1;i<=9;++i) for(int j=1;j<=9;++j){ int ii_=45*((i*(i+1))/2-1); int jj=(j*(j+1))/2; sum+=ww[ii_+jj]; }
    en.assign(172,0.0);
    int li=natorb[ni], lj=natorb[nj];
    elenuc(1,li,li+1,li+lj,en);
    if (ni==1 && nj==2) {
        fprintf(stderr, "[EN] ci1=%+.6f en0=%.4f en1=%.4f en2=%.4f en3=%.4f en4=%.4f en5=%.4f en6=%.4f\n",
                ci[1], en[0], en[1], en[2], en[3], en[4], en[5], en[6]);
        fflush(stderr);
    }
    ccrep(ni,nj,r,gab,enuc);
}

void rotmat(int nj,int ni,const std::vector<double>& coordi,const std::vector<double>& coordj,double& r){
    const double small=1e-7, pt5sq3=0.8660254037841;
    std::vector<std::vector<double>> p(4,std::vector<double>(4,0.0));
    std::vector<std::vector<double>> d(6,std::vector<double>(6,0.0));
    double x11=coordj[1]-coordi[1], x22=coordj[2]-coordi[2], x33=coordj[3]-coordi[3];
    double b=x11*x11+x22*x22;
    r=std::sqrt(b+x33*x33);
    double sqb=std::sqrt(b), sb=sqb/r;
    double ca=0,sa=0,cb=0;
    if(sb>small){ ca=x11/sqb; sa=x22/sqb; cb=x33/r; }
    else { sb=0;
        if(x33<0){ ca=-1; cb=-1; } else if(x33>0){ ca=1; cb=1; } else { ca=0; cb=0; }
    }
    p[1][1]=ca*sb; p[2][1]=ca*cb; p[3][1]=-sa;
    p[1][2]=sa*sb; p[2][2]=sa*cb; p[3][2]=ca;
    p[1][3]=cb;     p[2][3]=-sb;   p[3][3]=0;
    double c2a=0,c2b=0,s2a=0,s2b=0;
    if(dorbs[ni]||dorbs[nj]){
        c2a=2*ca*ca-1; c2b=2*cb*cb-1; s2a=2*sa*ca; s2b=2*sb*cb;
        d[1][1]=pt5sq3*c2a*sb*sb; d[2][1]=0.5*c2a*s2b; d[3][1]=-s2a*sb;
        d[4][1]=c2a*(cb*cb+0.5*sb*sb); d[5][1]=-s2a*cb;
        d[1][2]=pt5sq3*ca*s2b; d[2][2]=ca*c2b; d[3][2]=-sa*cb; d[4][2]=-0.5*ca*s2b; d[5][2]=sa*sb;
        d[1][3]=cb*cb-0.5*sb*sb; d[2][3]=-pt5sq3*s2b; d[3][3]=0; d[4][3]=pt5sq3*sb*sb; d[5][3]=0;
        d[1][4]=pt5sq3*sa*s2b; d[2][4]=sa*c2b; d[3][4]=ca*cb; d[4][4]=-0.5*sa*s2b; d[5][4]=-ca*sb;
        d[1][5]=pt5sq3*s2a*sb*sb; d[2][5]=0.5*s2a*s2b; d[3][5]=c2a*sb;
        d[4][5]=s2a*(cb*cb+0.5*sb*sb); d[5][5]=c2a*cb;
    }
    for(int k=1;k<=3;++k) for(int j=1;j<=3;++j) sp[k][j]=p[k][j];
    for(int k=1;k<=3;++k){
        pp[1][k][k]=p[k][1]*p[k][1];
        pp[2][k][k]=p[k][2]*p[k][2];
        pp[3][k][k]=p[k][3]*p[k][3];
        pp[4][k][k]=p[k][1]*p[k][2];
        pp[5][k][k]=p[k][1]*p[k][3];
        pp[6][k][k]=p[k][2]*p[k][3];
        if(k==1) continue;
        for(int j=1;j<=k-1;++j){
            pp[1][k][j]=2*p[k][1]*p[j][1];
            pp[2][k][j]=2*p[k][2]*p[j][2];
            pp[3][k][j]=2*p[k][3]*p[j][3];
            pp[4][k][j]=p[k][1]*p[j][2]+p[k][2]*p[j][1];
            pp[5][k][j]=p[k][1]*p[j][3]+p[k][3]*p[j][1];
            pp[6][k][j]=p[k][2]*p[j][3]+p[k][3]*p[j][2];
        }
    }
    if(dorbs[ni]||dorbs[nj]){
        for(int a=1;a<=5;++a) for(int b=1;b<=5;++b) sd[a][b]=d[a][b];
        for(int k=1;k<=5;++k){
            for(int t=1;t<=3;++t){
                dp[1][k][t]=d[k][1]*p[t][1]; dp[2][k][t]=d[k][1]*p[t][2]; dp[3][k][t]=d[k][1]*p[t][3];
                dp[4][k][t]=d[k][2]*p[t][1]; dp[5][k][t]=d[k][2]*p[t][2]; dp[6][k][t]=d[k][2]*p[t][3];
                dp[7][k][t]=d[k][3]*p[t][1]; dp[8][k][t]=d[k][3]*p[t][2]; dp[9][k][t]=d[k][3]*p[t][3];
                dp[10][k][t]=d[k][4]*p[t][1]; dp[11][k][t]=d[k][4]*p[t][2]; dp[12][k][t]=d[k][4]*p[t][3];
                dp[13][k][t]=d[k][5]*p[t][1]; dp[14][k][t]=d[k][5]*p[t][2]; dp[15][k][t]=d[k][5]*p[t][3];
            }
        }
        for(int k=1;k<=5;++k){
            d_d[1][k][k]=d[k][1]*d[k][1]; d_d[2][k][k]=d[k][2]*d[k][2];
            d_d[3][k][k]=d[k][3]*d[k][3]; d_d[4][k][k]=d[k][4]*d[k][4];
            d_d[5][k][k]=d[k][5]*d[k][5];
            d_d[6][k][k]=d[k][1]*d[k][2]; d_d[7][k][k]=d[k][1]*d[k][3]; d_d[8][k][k]=d[k][2]*d[k][3];
            d_d[9][k][k]=d[k][1]*d[k][4]; d_d[10][k][k]=d[k][2]*d[k][4]; d_d[11][k][k]=d[k][3]*d[k][4];
            d_d[12][k][k]=d[k][1]*d[k][5]; d_d[13][k][k]=d[k][2]*d[k][5];
            d_d[14][k][k]=d[k][3]*d[k][5]; d_d[15][k][k]=d[k][4]*d[k][5];
            if(k==1) continue;
            for(int j=1;j<=k-1;++j){
                d_d[1][k][j]=2*d[k][1]*d[j][1]; d_d[2][k][j]=2*d[k][2]*d[j][2];
                d_d[3][k][j]=2*d[k][3]*d[j][3]; d_d[4][k][j]=2*d[k][4]*d[j][4];
                d_d[5][k][j]=2*d[k][5]*d[j][5];
                d_d[6][k][j]=d[k][1]*d[j][2]+d[k][2]*d[j][1];
                d_d[7][k][j]=d[k][1]*d[j][3]+d[k][3]*d[j][1];
                d_d[8][k][j]=d[k][2]*d[j][3]+d[k][3]*d[j][2];
                d_d[9][k][j]=d[k][1]*d[j][4]+d[k][4]*d[j][1];
                d_d[10][k][j]=d[k][2]*d[j][4]+d[k][4]*d[j][2];
                d_d[11][k][j]=d[k][3]*d[j][4]+d[k][4]*d[j][3];
                d_d[12][k][j]=d[k][1]*d[j][5]+d[k][5]*d[j][1];
                d_d[13][k][j]=d[k][2]*d[j][5]+d[k][5]*d[j][2];
                d_d[14][k][j]=d[k][3]*d[j][5]+d[k][5]*d[j][3];
                d_d[15][k][j]=d[k][4]*d[j][5]+d[k][5]*d[j][4];
            }
        }
    }
}

void spcore(int ni,int nj,double r,double core[11][3]){
    double pxy[8]={0,1.0,-0.5,-0.5,0.5,0.25,0.25,0.5};
    for(int q=1;q<=4;++q){ core[q][1]=0; core[q][2]=0; }
    double r2=r*r;
    double aci=po[9][ni], acj=po[9][nj];
    double ssi=(aci+po[1][nj])*(aci+po[1][nj]);
    double ssj=(acj+po[1][ni])*(acj+po[1][ni]);
    core[1][1]=-tore[nj]*ev/std::sqrt(r2+ssj);
    core[1][2]=-tore[ni]*ev/std::sqrt(r2+ssi);
    if(ni>=3||nj>=3){
        if(ni>=3){
            double ppj=(acj+po[7][ni])*(acj+po[7][ni]);
            double da=ddp[2][ni]; double qa=ddp[3][ni]/std::sqrt(2.0);
            double twoqa=qa+qa;
            double adj=(po[2][ni]+acj)*(po[2][ni]+acj);
            double aqj=(po[3][ni]+acj)*(po[3][ni]+acj);
            double xj[8];
            xj[1]=r2+ppj; xj[2]=r2+aqj;
            xj[3]=(r+da)*(r+da)+adj; xj[4]=(r-da)*(r-da)+adj;
            xj[5]=(r-twoqa)*(r-twoqa)+aqj; xj[6]=(r+twoqa)*(r+twoqa)+aqj;
            xj[7]=r2+twoqa*twoqa+aqj;
            for(int i=1;i<=7;++i) xj[i]=pxy[i]/std::sqrt(xj[i]);
            double aj2=(xj[3]+xj[4])*ev;
            double aj3=(xj[1]+xj[2]+xj[5]+xj[6])*ev;
            double aj4=(xj[1]+xj[2]+xj[7])*ev;
            core[2][1]=-tore[nj]*aj2; core[3][1]=-tore[nj]*aj3; core[4][1]=-tore[nj]*aj4;
        }
        if(nj>=3){
            double ppi=(aci+po[7][nj])*(aci+po[7][nj]);
            double db=ddp[2][nj]; double qb=ddp[3][nj]/std::sqrt(2.0);
            double adi=(po[2][nj]+aci)*(po[2][nj]+aci);
            double aqi=(po[3][nj]+aci)*(po[3][nj]+aci);
            double twoqb=qb+qb;
            double xi[8];
            xi[1]=r2+ppi; xi[2]=r2+aqi;
            xi[3]=(r+db)*(r+db)+adi; xi[4]=(r-db)*(r-db)+adi;
            xi[5]=(r-twoqb)*(r-twoqb)+aqi; xi[6]=(r+twoqb)*(r+twoqb)+aqi;
            xi[7]=r2+twoqb*twoqb+aqi;
            for(int i=1;i<=7;++i) xi[i]=pxy[i]/std::sqrt(xi[i]);
            double ai2=-(xi[3]+xi[4])*ev;
            double ai3=(xi[1]+xi[2]+xi[5]+xi[6])*ev;
            double ai4=(xi[1]+xi[2]+xi[7])*ev;
            core[2][2]=-tore[ni]*ai2; core[3][2]=-tore[ni]*ai3; core[4][2]=-tore[ni]*ai4;
        }
    }
}

void tx(int ii,int kk,const std::vector<double>& rep,char* logv,double* v){
    int met[46]={0,1,2,3,2,3,3,2,3,3,3,4,5,5,5,6,4,5,5,5,6,6,4,5,5,5,6,6,6,4,5,5,5,6,6,6,6,4,5,5,5,6,6,6,6,6}; // F90 data met ends with 5*6
    int limkl=indx[kk][kk];
    for(int a=1;a<=45;++a) for(int b=1;b<=limkl;++b){ logv[a*46+b]=false; v[a*46+b]=0.0; }
    for(int i1=1;i1<=ii;++i1) for(int j1=1;j1<=i1;++j1){
        int ij=indexd[i1][j1];
        for(int k1=1;k1<=kk;++k1) for(int l1=1;l1<=k1;++l1){
            int kl=indexd[k1][l1]; int nd=ind2[ij][kl];
            if(nd==0) continue;
            double wrepp=rep[nd];
            int ll=indx[k1][l1]; int mm=met[ll];
            if(mm==1){ v[ij*46 + 1]=wrepp; }
            else if(mm==2){ int k=k1-1;
                v[ij*46 + 2]+=sp[k][1]*wrepp; v[ij*46 + 4]+=sp[k][2]*wrepp; v[ij*46 + 7]+=sp[k][3]*wrepp; }
            else if(mm==3){ int k=k1-1,l=l1-1;
                v[ij*46 + 3]+=pp[1][k][l]*wrepp; v[ij*46 + 6]+=pp[2][k][l]*wrepp; v[ij*46 + 10]+=pp[3][k][l]*wrepp;
                v[ij*46 + 5]+=pp[4][k][l]*wrepp; v[ij*46 + 8]+=pp[5][k][l]*wrepp; v[ij*46 + 9]+=pp[6][k][l]*wrepp;
                continue; }
            else if(mm==4){ int k=k1-4;
                v[ij*46 + 11]+=sd[k][1]*wrepp; v[ij*46 + 16]+=sd[k][2]*wrepp; v[ij*46 + 22]+=sd[k][3]*wrepp;
                v[ij*46 + 29]+=sd[k][4]*wrepp; v[ij*46 + 37]+=sd[k][5]*wrepp; }
            else if(mm==5){ int k=k1-4,l=l1-1;
                v[ij*46 + 12]+=dp[1][k][l]*wrepp; v[ij*46 + 13]+=dp[2][k][l]*wrepp; v[ij*46 + 14]+=dp[3][k][l]*wrepp;
                v[ij*46 + 17]+=dp[4][k][l]*wrepp; v[ij*46 + 18]+=dp[5][k][l]*wrepp; v[ij*46 + 19]+=dp[6][k][l]*wrepp;
                v[ij*46 + 23]+=dp[7][k][l]*wrepp; v[ij*46 + 24]+=dp[8][k][l]*wrepp; v[ij*46 + 25]+=dp[9][k][l]*wrepp;
                v[ij*46 + 30]+=dp[10][k][l]*wrepp; v[ij*46 + 31]+=dp[11][k][l]*wrepp; v[ij*46 + 32]+=dp[12][k][l]*wrepp;
                v[ij*46 + 38]+=dp[13][k][l]*wrepp; v[ij*46 + 39]+=dp[14][k][l]*wrepp; v[ij*46 + 40]+=dp[15][k][l]*wrepp; }
            else if(mm==6){ int k=k1-4,l=l1-4;
                v[ij*46 + 15]+=d_d[1][k][l]*wrepp; v[ij*46 + 21]+=d_d[2][k][l]*wrepp; v[ij*46 + 28]+=d_d[3][k][l]*wrepp;
                v[ij*46 + 36]+=d_d[4][k][l]*wrepp; v[ij*46 + 45]+=d_d[5][k][l]*wrepp;
                v[ij*46 + 20]+=d_d[6][k][l]*wrepp; v[ij*46 + 26]+=d_d[7][k][l]*wrepp; v[ij*46 + 27]+=d_d[8][k][l]*wrepp;
                v[ij*46 + 33]+=d_d[9][k][l]*wrepp; v[ij*46 + 34]+=d_d[10][k][l]*wrepp; v[ij*46 + 35]+=d_d[11][k][l]*wrepp;
                v[ij*46 + 41]+=d_d[12][k][l]*wrepp; v[ij*46 + 42]+=d_d[13][k][l]*wrepp; v[ij*46 + 43]+=d_d[14][k][l]*wrepp;
                v[ij*46 + 44]+=d_d[15][k][l]*wrepp; }
        }
        for(int b=1;b<=limkl;++b) if(v[ij*46 + b]!=0.0) logv[ij*46 + b]=true;
    }
}

void w2mat(const std::vector<double>& ww,std::vector<double>& w,int& kr,int limij,int limkl){
    // Official w2mat (mndod.F90 L2017-2052): ww declared (limkl,limij), column-major,
    // so ww(kl,ij) linearizes to (ij-1)*limkl+kl == l as l increments (ij outer, kl inner).
    // => w(l)=ww(l): straight copy, NO transpose.
    int l=0;
    for(int ij=1;ij<=limij;++ij) for(int kl=1;kl<=limkl;++kl){ ++l; w[l]=ww[l]; }
    kr+=l;
}

namespace mndod_C {
extern std::vector<int> intij,intkl,intrep;
extern std::vector<std::vector<double>> repd;
}
using mndod_C::intij; using mndod_C::intkl; using mndod_C::intrep; using mndod_C::repd;

void wstore(int ni,int ilim,int& kr,std::vector<std::vector<double>>& w){
    for(int i=1;i<=ilim;++i) for(int j=1;j<=ilim;++j) w[i][j]=0.0;
    int ip=1; w[ip][ip]=gss[ni];
    if(natorb[ni]>2){
        int ipx=ip+2, ipy=ip+5, ipz=ip+9;
        w[ipx][ip]=gsp[ni]; w[ipy][ip]=gsp[ni]; w[ipz][ip]=gsp[ni];
        w[ip][ipx]=gsp[ni]; w[ip][ipy]=gsp[ni]; w[ip][ipz]=gsp[ni];
        w[ipx][ipx]=gpp[ni]; w[ipy][ipy]=gpp[ni]; w[ipz][ipz]=gpp[ni];
        w[ipy][ipx]=gp2[ni]; w[ipz][ipx]=gp2[ni]; w[ipz][ipy]=gp2[ni];
        w[ipx][ipy]=gp2[ni]; w[ipx][ipz]=gp2[ni]; w[ipy][ipz]=gp2[ni];
        w[ip+1][ip+1]=hsp[ni]; w[ip+3][ip+3]=hsp[ni]; w[ip+6][ip+6]=hsp[ni];
        w[ip+4][ip+4]=0.5*(gpp[ni]-gp2[ni]); w[ip+7][ip+7]=0.5*(gpp[ni]-gp2[ni]);
        w[ip+8][ip+8]=0.5*(gpp[ni]-gp2[ni]);
        if(ilim>10){
            int ij0=ip-1;
            for(int i=1;i<=243;++i){
                int ij=intij[i], kl=intkl[i], Int=intrep[i];
                w[ij+ij0][kl+ij0]=repd[Int][ni];
            }
        }
    }
    kr+=ilim*ilim;
}

double charg(double r,int l1,int l2,int m,double da,double db,double add){
    double res=0.0;
    auto S=[&](double x)->double { return 1.0/std::sqrt(x*x+add); };
    if(l1==0&&l2==0){ res=S(r); }
    else if(l1==1&&l2==0){ res=(-S(r+da)+S(r-da))/2.0; }
    else if(l1==0&&l2==1){ res=(S(r+db)-S(r-db))/2.0; }
    else if(l1==1&&l2==1&&m==0){
        double dzdz=S(r+da-db)+S(r-da+db)-S(r-da-db)-S(r+da+db);
        res=dzdz/4.0;
    }
    else if(l1==1&&l2==1&&m==1){
        double dxdx=2.0*S(std::sqrt((da-db)*(da-db)+r*r))-2.0*S(std::sqrt((da+db)*(da+db)+r*r));
        res=dxdx*0.25;
    }
    else if(l1==0&&l2==2){
        double qqzz=S(r-db)-2.0*S(std::sqrt(r*r+db*db))+S(r+db);
        res=qqzz/4.0;
    }
    else if(l1==2&&l2==0){
        double qzzq=S(r-da)-2.0*S(std::sqrt(r*r+da*da))+S(r+da);
        res=qzzq/4.0;
    }
    else if(l1==1&&l2==2&&m==0){
        double dzqzz=S(r-da-db)-2.0*S(std::sqrt((r-da)*(r-da)+db*db))+S(r+db-da)
            -S(r-db+da)+2.0*S(std::sqrt((r+da)*(r+da)+db*db))-S(r+da+db);
        res=dzqzz/8.0;
    }
    else if(l1==2&&l2==1&&m==0){
        double qzzdz=-S(r-da-db)+2.0*S(std::sqrt((r-db)*(r-db)+da*da))-S(r+da-db)
            +S(r-da+db)-2.0*S(std::sqrt((r+db)*(r+db)+da*da))+S(r+da+db);
        res=qzzdz/8.0;
    }
    else if(l1==2&&l2==2&&m==0){
        double zzzz=S(r-da-db)+S(r+da+db)+S(r-da+db)+S(r+da-db)
            -2.0*S(std::sqrt((r-da)*(r-da)+db*db))-2.0*S(std::sqrt((r-db)*(r-db)+da*da))
            -2.0*S(std::sqrt((r+da)*(r+da)+db*db))-2.0*S(std::sqrt((r+db)*(r+db)+da*da))
            +2.0*S(std::sqrt(r*r+(da-db)*(da-db)))+2.0*S(std::sqrt(r*r+(da+db)*(da+db)));
        double xyxy=4.0*S(std::sqrt(r*r+(da-db)*(da-db)))+4.0*S(std::sqrt(r*r+(da+db)*(da+db)))
            -8.0*S(std::sqrt(r*r+da*da+db*db));
        res=zzzz/16.0-xyxy/64.0;
    }
    else if(l1==1&&l2==2&&m==1){
        double ab=db/std::sqrt(2.0);
        double dxqxz=-2.0*S(std::sqrt((r-ab)*(r-ab)+(da-ab)*(da-ab)))+2.0*S(std::sqrt((r+ab)*(r+ab)+(da-ab)*(da-ab)))
            +2.0*S(std::sqrt((r-ab)*(r-ab)+(da+ab)*(da+ab)))-2.0*S(std::sqrt((r+ab)*(r+ab)+(da+ab)*(da+ab)));
        res=dxqxz/8.0;
    }
    else if(l1==2&&l2==1&&m==1){
        double aa=da/std::sqrt(2.0);
        double qxzdx=-2.0*S(std::sqrt((r+aa)*(r+aa)+(aa-db)*(aa-db)))+2.0*S(std::sqrt((r-aa)*(r-aa)+(aa-db)*(aa-db)))
            +2.0*S(std::sqrt((r+aa)*(r+aa)+(aa+db)*(aa+db)))-2.0*S(std::sqrt((r-aa)*(r-aa)+(aa+db)*(aa+db)));
        res=qxzdx/8.0;
    }
    else if(l1==2&&l2==2&&m==1){
        double aa=da/std::sqrt(2.0), ab=db/std::sqrt(2.0);
        auto S2=[&](double x,double y)->double { return 1.0/std::sqrt(x*x+y*y+add); };
        double qxzqxz=2.0*S2(r+aa-ab,aa-ab)-2.0*S2(r+aa+ab,aa-ab)-2.0*S2(r-aa-ab,aa-ab)+2.0*S2(r-aa+ab,aa-ab)
            -2.0*S2(r+aa-ab,aa+ab)+2.0*S2(r+aa+ab,aa+ab)+2.0*S2(r-aa-ab,aa+ab)-2.0*S2(r-aa+ab,aa+ab);
        res=qxzqxz/16.0;
    }
    else if(l1==2&&l2==2&&m==2){
        double xyxy=4.0*S(std::sqrt(r*r+(da-db)*(da-db)))+4.0*S(std::sqrt(r*r+(da+db)*(da+db)))
            -8.0*S(std::sqrt(r*r+da*da+db*db));
        res=xyxy/16.0;
    }
    return res;
}

namespace mndod_C {
extern std::vector<std::vector<int>> indpp, inddp, inddd;
}


using mndod_C::indpp; using mndod_C::inddp; using mndod_C::inddd;
void elenuc(int ia,int ib,int ja,int jb,std::vector<double>& h){
    int k=ia,l=ib,n=1;
    while(true){
        for(int i=k;i<=l;++i){
            int ind1=i-k;
            if(ind1==0){
                for(int j=k;j<=i;++j){
                    int ind2=j-k; int m=(i*(i-1))/2+j;
                    if(ind2==0) h[m]+=cored[1][n];
                    else if(ind2<4){ int ipp=indpp[ind1][ind2];
                        h[m]+=cored[3][n]*pp[ipp][1][1]+cored[4][n]*(pp[ipp][2][2]+pp[ipp][3][3]); }
                    else { int idd=inddd[ind1-3][ind2-3];
                        h[m]+=cored[7][n]*d_d[idd][1][1]+cored[9][n]*(d_d[idd][2][2]+d_d[idd][3][3])
                            +cored[10][n]*(d_d[idd][4][4]+d_d[idd][5][5]); }
                }
            } else if(ind1<4){
                for(int j=k;j<=i;++j){
                    int ind2=j-k; int m=(i*(i-1))/2+j;
                    if(ind2==0) h[m]+=sp[1][ind1]*cored[2][n];
                    else if(ind2<4){ int ipp=indpp[ind1][ind2];
                        h[m]+=cored[3][n]*pp[ipp][1][1]+cored[4][n]*(pp[ipp][2][2]+pp[ipp][3][3]); }
                    else { int idd=inddd[ind1-3][ind2-3];
                        h[m]+=cored[7][n]*d_d[idd][1][1]+cored[9][n]*(d_d[idd][2][2]+d_d[idd][3][3])
                            +cored[10][n]*(d_d[idd][4][4]+d_d[idd][5][5]); }
                }
            } else {
                for(int j=k;j<=i;++j){
                    int ind2=j-k; int m=(i*(i-1))/2+j;
                    if(ind2==0) h[m]+=sd[1][ind1-3]*cored[5][n];
                    else if(ind2<4){ int idp=inddp[ind1-3][ind2];
                        h[m]+=cored[6][n]*dp[idp][1][1]+cored[8][n]*(dp[idp][2][2]+dp[idp][3][3]); }
                    else { int idd=inddd[ind1-3][ind2-3];
                        h[m]+=cored[7][n]*d_d[idd][1][1]+cored[9][n]*(d_d[idd][2][2]+d_d[idd][3][3])
                            +cored[10][n]*(d_d[idd][4][4]+d_d[idd][5][5]); }
                }
            }
        }
        if(n==2) break;
        k=ja; l=jb; n=2;
    }
}

namespace mndod_C {
extern std::vector<int> iii, iiid;
}
extern double rsc(int,int,double,int,double,int,double,int,double);
using mndod_C::iii; using mndod_C::iiid;

void scprm(int ni,double& r066,double& r266,double& r466,double& r016,double& r244,
           double& r036,double& r236,double& r155,double& r355,double& r125,double& r234,double& r246){
    int ns=iii[ni], nd=iiid[ni];
    double es=zsn[ni], ep=zpn[ni], ed=zdn[ni];
    r016=rsc(0,ns,es,ns,es,nd,ed,nd,ed);
    r036=rsc(0,ns,ep,ns,ep,nd,ed,nd,ed);
    r066=rsc(0,nd,ed,nd,ed,nd,ed,nd,ed);
    r155=rsc(1,ns,ep,nd,ed,ns,ep,nd,ed);
    r125=rsc(1,ns,es,ns,ep,ns,ep,nd,ed);
    r244=rsc(2,ns,es,nd,ed,ns,es,nd,ed);
    r236=rsc(2,ns,ep,ns,ep,nd,ed,nd,ed);
    r266=rsc(2,nd,ed,nd,ed,nd,ed,nd,ed);
    r234=rsc(2,ns,ep,ns,ep,ns,es,nd,ed);
    r246=rsc(2,ns,es,nd,ed,nd,ed,nd,ed);
    r355=rsc(3,ns,ep,nd,ed,ns,ep,nd,ed);
    r466=rsc(4,nd,ed,nd,ed,nd,ed,nd,ed);
}

void eiscor(double r016,double r066,double r244,double r266,double r466,int ni){
    int ir016[101],ir066[101],ir244[101],ir266[101],ir466[101];
    std::fill(ir016,ir016+101,0); std::fill(ir066,ir066+101,0);
    std::fill(ir244,ir244+101,0); std::fill(ir266,ir266+101,0); std::fill(ir466,ir466+101,0);
    int a16[9]={0,2,4,6,5,10,12,14,16};
    int a66[9]={0,0,1,3,10,10,15,21,28};
    int a44[9]={0,1,2,3,5,5,6,7,8};
    int a66_[9]={0,0,8,15,35,35,35,43,50};
    int a46[9]={0,0,1,8,35,35,35,36,43};
    for(int i=1;i<=8;++i){ ir016[20+i]=a16[i]; ir066[20+i]=a66[i]; ir244[20+i]=a44[i]; ir266[20+i]=a66_[i]; ir466[20+i]=a46[i]; }
    ir016[29]=10; ir066[29]=45; ir244[29]=5; ir266[29]=70; ir466[29]=70;
    int b16[9]={0,2,4,4,5,10,7,8,0};
    int b66[9]={0,0,1,6,10,10,21,28,45};
    int b44[9]={0,1,2,4,5,5,5,5,0};
    int b66_[9]={0,0,8,21,35,35,43,50,70};
    int b46[9]={0,0,1,21,35,35,36,43,70};
    for(int i=1;i<=8;++i){ ir016[38+i]=b16[i]; ir066[38+i]=b66[i]; ir244[38+i]=b44[i]; ir266[38+i]=b66_[i]; ir466[38+i]=b46[i]; }
    ir016[47]=10; ir066[47]=45; ir244[47]=5; ir266[47]=70; ir466[47]=70;
    ir016[57]=2; ir066[57]=0; ir244[57]=1; ir266[57]=0; ir466[57]=0;
    ir016[71]=2; ir066[71]=0; ir244[71]=1; ir266[71]=0; ir466[71]=0;
    int c16[9]={0,4,6,5,10,12,14,9,10};
    int c66[9]={0,1,3,10,10,15,21,36,45};
    int c44[9]={0,2,3,5,5,6,7,5,5};
    int c66_[9]={0,8,15,35,35,35,43,56,70};
    int c46[9]={0,1,8,35,35,35,36,56,70};
    for(int i=1;i<=8;++i){ ir016[71+i]=c16[i]; ir066[71+i]=c66[i]; ir244[71+i]=c44[i]; ir266[71+i]=c66_[i]; ir466[71+i]=c46[i]; }
    ir016[80]=0; ir066[80]=0; ir244[80]=0; ir266[80]=0; ir466[80]=0;
    eisol[ni]=eisol[ni]+ir016[ni]*r016+ir066[ni]*r066-ir244[ni]*r244/5.0
        -ir266[ni]*r266/49.0-ir466[ni]*r466/49.0;
}

namespace mndod_C {
std::vector<double> fx(31,0.0);
std::vector<std::vector<double>> b(31,std::vector<double>(31,0.0));
}
using mndod_C::fx; using mndod_C::b;

void fbx(){
    fx[1]=1.0;
    for(int i=2;i<=30;++i) fx[i]=fx[i-1]*(i-1);
    for(int j=1;j<=30;++j) b[j][1]=1.0;
    for(int i=2;i<=30;++i) for(int j=2;j<=i;++j) b[i][j]=b[i-1][j-1]+b[i-1][j];
}

double rsc(int k,int na,double ea,int nb,double eb,int nc,double ec,int nd,double ed){
    double aea=std::log(ea), aeb=std::log(eb), aec=std::log(ec), aed=std::log(ed);
    int nab=na+nb, ncd=nc+nd;
    double ecd=ec+ed, eab=ea+eb, e=ecd+eab;
    int n=nab+ncd;
    double ae=std::log(e), a2=std::log(2.0), acd=std::log(ecd), aab=std::log(eab);
    double ff=fx[n]/std::sqrt(fx[2*na+1]*fx[2*nb+1]*fx[2*nc+1]*fx[2*nd+1]);
    double c=ev*ff*std::exp(na*aea+nb*aeb+nc*aec+nd*aed+0.5*(aea+aeb+aec+aed)+a2*(n+2)-ae*n);
    double s0=1.0/e, s1=0.0, s2=0.0;
    int m=ncd-k;
    for(int i=1;i<=m;++i){
        s0=s0*e/ecd;
        s1+=s0*(b[ncd-k][i]-b[ncd+k+1][i])/b[n][i];
    }
    int m1=m+1, m2=ncd+k+1;
    for(int i=m1;i<=m2;++i){
        s0=s0*e/ecd;
        s2+=s0*b[m2][i]/b[n][i];
    }
    double s3=std::exp(ae*n-acd*m2-aab*(nab-k))/b[n][m2];
    return c*(s1-s2+s3);
}
namespace mndod_C {
std::vector<std::vector<int>> indexd(10,std::vector<int>(10,0));
std::vector<std::vector<int>> indx(10,std::vector<int>(10,0));
std::vector<std::vector<int>> ind2(46,std::vector<int>(46,0));
std::vector<int> isym(492,0);
std::vector<std::vector<std::vector<double>>> ch(46,std::vector<std::vector<double>>(3,std::vector<double>(5,0.0)));
std::vector<std::vector<int>> indpp(4,std::vector<int>(4,0));
std::vector<std::vector<int>> inddp(6,std::vector<int>(4,0));
std::vector<std::vector<int>> inddd(6,std::vector<int>(6,0));
}
void fordd(){
    for(int i=1;i<=9;++i) for(int j=1;j<=i;++j){
        indexd[i][j]=-(j*(j-1))/2+i+9*(j-1);
        indx[i][j]=(i*(i-1))/2+j;
        indexd[j][i]=indexd[i][j];
        indx[j][i]=indx[i][j];
    }
    ind2[1][1]=1;
    ind2[1][2]=2;
    ind2[1][10]=3;
    ind2[1][18]=4;
    ind2[1][25]=5;
    ind2[2][1]=6;
    ind2[2][2]=7;
    ind2[2][10]=8;
    ind2[2][18]=9;
    ind2[2][25]=10;
    ind2[10][1]=11;
    ind2[10][2]=12;
    ind2[10][10]=13;
    ind2[10][18]=14;
    ind2[10][25]=15;
    ind2[3][3]=16;
    ind2[3][11]=17;
    ind2[11][3]=18;
    ind2[11][11]=19;
    ind2[18][1]=20;
    ind2[18][2]=21;
    ind2[18][10]=22;
    ind2[18][18]=23;
    ind2[18][25]=24;
    ind2[4][4]=25;
    ind2[4][12]=26;
    ind2[12][4]=27;
    ind2[12][12]=28;
    ind2[19][19]=29;
    ind2[25][1]=30;
    ind2[25][2]=31;
    ind2[25][10]=32;
    ind2[25][18]=33;
    ind2[25][25]=34;
    ind2[1][5]=35;
    ind2[1][13]=36;
    ind2[1][31]=37;
    ind2[1][21]=38;
    ind2[1][36]=39;
    ind2[1][28]=40;
    ind2[1][40]=41;
    ind2[1][43]=42;
    ind2[1][45]=43;
    ind2[2][5]=44;
    ind2[2][13]=45;
    ind2[2][31]=46;
    ind2[2][21]=47;
    ind2[2][36]=48;
    ind2[2][28]=49;
    ind2[2][40]=50;
    ind2[2][43]=51;
    ind2[2][45]=52;
    ind2[10][5]=53;
    ind2[10][13]=54;
    ind2[10][31]=55;
    ind2[10][21]=56;
    ind2[10][36]=57;
    ind2[10][28]=58;
    ind2[10][40]=59;
    ind2[10][43]=60;
    ind2[10][45]=61;
    ind2[3][20]=62;
    ind2[3][6]=63;
    ind2[3][14]=64;
    ind2[3][32]=65;
    ind2[3][23]=66;
    ind2[3][38]=67;
    ind2[3][30]=68;
    ind2[3][42]=69;
    ind2[11][20]=70;
    ind2[11][6]=71;
    ind2[11][14]=72;
    ind2[11][32]=73;
    ind2[11][23]=74;
    ind2[11][38]=75;
    ind2[11][30]=76;
    ind2[11][42]=77;
    ind2[18][5]=78;
    ind2[18][13]=79;
    ind2[18][31]=80;
    ind2[18][21]=81;
    ind2[18][36]=82;
    ind2[18][28]=83;
    ind2[18][40]=84;
    ind2[18][8]=85;
    ind2[18][16]=86;
    ind2[18][34]=87;
    ind2[18][43]=88;
    ind2[18][45]=89;
    ind2[4][26]=90;
    ind2[4][7]=91;
    ind2[4][15]=92;
    ind2[4][33]=93;
    ind2[4][29]=94;
    ind2[4][41]=95;
    ind2[4][24]=96;
    ind2[4][39]=97;
    ind2[12][26]=98;
    ind2[12][7]=99;
    ind2[12][15]=100;
    ind2[12][33]=101;
    ind2[12][29]=102;
    ind2[12][41]=103;
    ind2[12][24]=104;
    ind2[12][39]=105;
    ind2[19][27]=106;
    ind2[19][22]=107;
    ind2[19][37]=108;
    ind2[19][9]=109;
    ind2[19][17]=110;
    ind2[19][35]=111;
    ind2[25][5]=112;
    ind2[25][13]=113;
    ind2[25][31]=114;
    ind2[25][21]=115;
    ind2[25][36]=116;
    ind2[25][28]=117;
    ind2[25][40]=118;
    ind2[25][8]=119;
    ind2[25][16]=120;
    ind2[25][34]=121;
    ind2[25][43]=122;
    ind2[25][45]=123;
    ind2[5][1]=124;
    ind2[5][2]=125;
    ind2[5][10]=126;
    ind2[5][18]=127;
    ind2[5][25]=128;
    ind2[5][5]=129;
    ind2[5][13]=130;
    ind2[5][31]=131;
    ind2[5][21]=132;
    ind2[5][36]=133;
    ind2[5][28]=134;
    ind2[5][40]=135;
    ind2[5][43]=136;
    ind2[5][45]=137;
    ind2[13][1]=138;
    ind2[13][2]=139;
    ind2[13][10]=140;
    ind2[13][18]=141;
    ind2[13][25]=142;
    ind2[13][5]=143;
    ind2[13][13]=144;
    ind2[13][31]=145;
    ind2[13][21]=146;
    ind2[13][36]=147;
    ind2[13][28]=148;
    ind2[13][40]=149;
    ind2[13][43]=150;
    ind2[13][45]=151;
    ind2[20][3]=152;
    ind2[20][11]=153;
    ind2[20][20]=154;
    ind2[20][6]=155;
    ind2[20][14]=156;
    ind2[20][32]=157;
    ind2[20][23]=158;
    ind2[20][38]=159;
    ind2[20][30]=160;
    ind2[20][42]=161;
    ind2[26][4]=162;
    ind2[26][12]=163;
    ind2[26][26]=164;
    ind2[26][7]=165;
    ind2[26][15]=166;
    ind2[26][33]=167;
    ind2[26][29]=168;
    ind2[26][41]=169;
    ind2[26][24]=170;
    ind2[26][39]=171;
    ind2[31][1]=172;
    ind2[31][2]=173;
    ind2[31][10]=174;
    ind2[31][18]=175;
    ind2[31][25]=176;
    ind2[31][5]=177;
    ind2[31][13]=178;
    ind2[31][31]=179;
    ind2[31][21]=180;
    ind2[31][36]=181;
    ind2[31][28]=182;
    ind2[31][40]=183;
    ind2[31][43]=184;
    ind2[31][45]=185;
    ind2[6][3]=186;
    ind2[6][11]=187;
    ind2[6][20]=188;
    ind2[6][6]=189;
    ind2[6][14]=190;
    ind2[6][32]=191;
    ind2[6][23]=192;
    ind2[6][38]=193;
    ind2[6][30]=194;
    ind2[6][42]=195;
    ind2[14][3]=196;
    ind2[14][11]=197;
    ind2[14][20]=198;
    ind2[14][6]=199;
    ind2[14][14]=200;
    ind2[14][32]=201;
    ind2[14][23]=202;
    ind2[14][38]=203;
    ind2[14][30]=204;
    ind2[14][42]=205;
    ind2[21][1]=206;
    ind2[21][2]=207;
    ind2[21][10]=208;
    ind2[21][18]=209;
    ind2[21][25]=210;
    ind2[21][5]=211;
    ind2[21][13]=212;
    ind2[21][31]=213;
    ind2[21][21]=214;
    ind2[21][36]=215;
    ind2[21][28]=216;
    ind2[21][40]=217;
    ind2[21][8]=218;
    ind2[21][16]=219;
    ind2[21][34]=220;
    ind2[21][43]=221;
    ind2[21][45]=222;
    ind2[27][19]=223;
    ind2[27][27]=224;
    ind2[27][22]=225;
    ind2[27][37]=226;
    ind2[27][9]=227;
    ind2[27][17]=228;
    ind2[27][35]=229;
    ind2[32][3]=230;
    ind2[32][11]=231;
    ind2[32][20]=232;
    ind2[32][6]=233;
    ind2[32][14]=234;
    ind2[32][32]=235;
    ind2[32][23]=236;
    ind2[32][38]=237;
    ind2[32][30]=238;
    ind2[32][42]=239;
    ind2[36][1]=240;
    ind2[36][2]=241;
    ind2[36][10]=242;
    ind2[36][18]=243;
    ind2[36][25]=244;
    ind2[36][5]=245;
    ind2[36][13]=246;
    ind2[36][31]=247;
    ind2[36][21]=248;
    ind2[36][36]=249;
    ind2[36][28]=250;
    ind2[36][40]=251;
    ind2[36][8]=252;
    ind2[36][16]=253;
    ind2[36][34]=254;
    ind2[36][43]=255;
    ind2[36][45]=256;
    ind2[7][4]=257;
    ind2[7][12]=258;
    ind2[7][26]=259;
    ind2[7][7]=260;
    ind2[7][15]=261;
    ind2[7][33]=262;
    ind2[7][29]=263;
    ind2[7][41]=264;
    ind2[7][24]=265;
    ind2[7][39]=266;
    ind2[15][4]=267;
    ind2[15][12]=268;
    ind2[15][26]=269;
    ind2[15][7]=270;
    ind2[15][15]=271;
    ind2[15][33]=272;
    ind2[15][29]=273;
    ind2[15][41]=274;
    ind2[15][24]=275;
    ind2[15][39]=276;
    ind2[22][19]=277;
    ind2[22][27]=278;
    ind2[22][22]=279;
    ind2[22][37]=280;
    ind2[22][9]=281;
    ind2[22][17]=282;
    ind2[22][35]=283;
    ind2[28][1]=284;
    ind2[28][2]=285;
    ind2[28][10]=286;
    ind2[28][18]=287;
    ind2[28][25]=288;
    ind2[28][5]=289;
    ind2[28][13]=290;
    ind2[28][31]=291;
    ind2[28][21]=292;
    ind2[28][36]=293;
    ind2[28][28]=294;
    ind2[28][40]=295;
    ind2[28][8]=296;
    ind2[28][16]=297;
    ind2[28][34]=298;
    ind2[28][43]=299;
    ind2[28][45]=300;
    ind2[33][4]=301;
    ind2[33][12]=302;
    ind2[33][26]=303;
    ind2[33][7]=304;
    ind2[33][15]=305;
    ind2[33][33]=306;
    ind2[33][29]=307;
    ind2[33][41]=308;
    ind2[33][24]=309;
    ind2[33][39]=310;
    ind2[37][19]=311;
    ind2[37][27]=312;
    ind2[37][22]=313;
    ind2[37][37]=314;
    ind2[37][9]=315;
    ind2[37][17]=316;
    ind2[37][35]=317;
    ind2[40][1]=318;
    ind2[40][2]=319;
    ind2[40][10]=320;
    ind2[40][18]=321;
    ind2[40][25]=322;
    ind2[40][5]=323;
    ind2[40][13]=324;
    ind2[40][31]=325;
    ind2[40][21]=326;
    ind2[40][36]=327;
    ind2[40][28]=328;
    ind2[40][40]=329;
    ind2[40][8]=330;
    ind2[40][16]=331;
    ind2[40][34]=332;
    ind2[40][43]=333;
    ind2[40][45]=334;
    ind2[8][18]=335;
    ind2[8][25]=336;
    ind2[8][21]=337;
    ind2[8][36]=338;
    ind2[8][28]=339;
    ind2[8][40]=340;
    ind2[8][8]=341;
    ind2[8][16]=342;
    ind2[8][34]=343;
    ind2[16][18]=344;
    ind2[16][25]=345;
    ind2[16][21]=346;
    ind2[16][36]=347;
    ind2[16][28]=348;
    ind2[16][40]=349;
    ind2[16][8]=350;
    ind2[16][16]=351;
    ind2[16][34]=352;
    ind2[23][3]=353;
    ind2[23][11]=354;
    ind2[23][20]=355;
    ind2[23][6]=356;
    ind2[23][14]=357;
    ind2[23][32]=358;
    ind2[23][23]=359;
    ind2[23][38]=360;
    ind2[23][30]=361;
    ind2[23][42]=362;
    ind2[29][4]=363;
    ind2[29][12]=364;
    ind2[29][26]=365;
    ind2[29][7]=366;
    ind2[29][15]=367;
    ind2[29][33]=368;
    ind2[29][29]=369;
    ind2[29][41]=370;
    ind2[29][24]=371;
    ind2[29][39]=372;
    ind2[34][18]=373;
    ind2[34][25]=374;
    ind2[34][21]=375;
    ind2[34][36]=376;
    ind2[34][28]=377;
    ind2[34][40]=378;
    ind2[34][8]=379;
    ind2[34][16]=380;
    ind2[34][34]=381;
    ind2[38][3]=382;
    ind2[38][11]=383;
    ind2[38][20]=384;
    ind2[38][6]=385;
    ind2[38][14]=386;
    ind2[38][32]=387;
    ind2[38][23]=388;
    ind2[38][38]=389;
    ind2[38][30]=390;
    ind2[38][42]=391;
    ind2[41][4]=392;
    ind2[41][12]=393;
    ind2[41][26]=394;
    ind2[41][7]=395;
    ind2[41][15]=396;
    ind2[41][33]=397;
    ind2[41][29]=398;
    ind2[41][41]=399;
    ind2[41][24]=400;
    ind2[41][39]=401;
    ind2[43][1]=402;
    ind2[43][2]=403;
    ind2[43][10]=404;
    ind2[43][18]=405;
    ind2[43][25]=406;
    ind2[43][5]=407;
    ind2[43][13]=408;
    ind2[43][31]=409;
    ind2[43][21]=410;
    ind2[43][36]=411;
    ind2[43][28]=412;
    ind2[43][40]=413;
    ind2[43][43]=414;
    ind2[43][45]=415;
    ind2[9][19]=416;
    ind2[9][27]=417;
    ind2[9][22]=418;
    ind2[9][37]=419;
    ind2[9][9]=420;
    ind2[9][17]=421;
    ind2[9][35]=422;
    ind2[17][19]=423;
    ind2[17][27]=424;
    ind2[17][22]=425;
    ind2[17][37]=426;
    ind2[17][9]=427;
    ind2[17][17]=428;
    ind2[17][35]=429;
    ind2[24][4]=430;
    ind2[24][12]=431;
    ind2[24][26]=432;
    ind2[24][7]=433;
    ind2[24][15]=434;
    ind2[24][33]=435;
    ind2[24][29]=436;
    ind2[24][41]=437;
    ind2[24][24]=438;
    ind2[24][39]=439;
    ind2[30][3]=440;
    ind2[30][11]=441;
    ind2[30][20]=442;
    ind2[30][6]=443;
    ind2[30][14]=444;
    ind2[30][32]=445;
    ind2[30][23]=446;
    ind2[30][38]=447;
    ind2[30][30]=448;
    ind2[30][42]=449;
    ind2[35][19]=450;
    ind2[35][27]=451;
    ind2[35][22]=452;
    ind2[35][37]=453;
    ind2[35][9]=454;
    ind2[35][17]=455;
    ind2[35][35]=456;
    ind2[39][4]=457;
    ind2[39][12]=458;
    ind2[39][26]=459;
    ind2[39][7]=460;
    ind2[39][15]=461;
    ind2[39][33]=462;
    ind2[39][29]=463;
    ind2[39][41]=464;
    ind2[39][24]=465;
    ind2[39][39]=466;
    ind2[42][3]=467;
    ind2[42][11]=468;
    ind2[42][20]=469;
    ind2[42][6]=470;
    ind2[42][14]=471;
    ind2[42][32]=472;
    ind2[42][23]=473;
    ind2[42][38]=474;
    ind2[42][30]=475;
    ind2[42][42]=476;
    ind2[44][44]=477;
    ind2[45][1]=478;
    ind2[45][2]=479;
    ind2[45][10]=480;
    ind2[45][18]=481;
    ind2[45][25]=482;
    ind2[45][5]=483;
    ind2[45][13]=484;
    ind2[45][31]=485;
    ind2[45][21]=486;
    ind2[45][36]=487;
    ind2[45][28]=488;
    ind2[45][40]=489;
    ind2[45][43]=490;
    ind2[45][45]=491;
    isym[40]=38;
    isym[41]=39;
    isym[43]=42;
    isym[49]=47;
    isym[50]=48;
    isym[52]=51;
    isym[58]=56;
    isym[59]=57;
    isym[61]=60;
    isym[68]=66;
    isym[69]=67;
    isym[76]=74;
    isym[77]=75;
    isym[89]=88;
    isym[90]=62;
    isym[91]=63;
    isym[92]=64;
    isym[93]=65;
    isym[94]=-66;
    isym[95]=-67;
    isym[96]=66;
    isym[97]=67;
    isym[98]=70;
    isym[99]=71;
    isym[100]=72;
    isym[101]=73;
    isym[102]=-74;
    isym[103]=-75;
    isym[104]=74;
    isym[105]=75;
    isym[106]=86;
    isym[107]=86;
    isym[109]=85;
    isym[110]=86;
    isym[111]=87;
    isym[112]=78;
    isym[113]=79;
    isym[114]=80;
    isym[115]=83;
    isym[116]=84;
    isym[117]=81;
    isym[118]=82;
    isym[119]=-85;
    isym[120]=-86;
    isym[121]=-87;
    isym[122]=88;
    isym[123]=88;
    isym[128]=127;
    isym[134]=132;
    isym[135]=133;
    isym[137]=136;
    isym[142]=141;
    isym[148]=146;
    isym[149]=147;
    isym[151]=150;
    isym[160]=158;
    isym[161]=159;
    isym[162]=152;
    isym[163]=153;
    isym[164]=154;
    isym[165]=155;
    isym[166]=156;
    isym[167]=157;
    isym[168]=-158;
    isym[169]=-159;
    isym[170]=158;
    isym[171]=159;
    isym[176]=175;
    isym[182]=180;
    isym[183]=181;
    isym[185]=184;
    isym[194]=192;
    isym[195]=193;
    isym[204]=202;
    isym[205]=203;
    isym[222]=221;
    isym[224]=219;
    isym[225]=219;
    isym[227]=218;
    isym[228]=219;
    isym[229]=220;
    isym[238]=236;
    isym[239]=237;
    isym[256]=255;
    isym[257]=186;
    isym[258]=187;
    isym[259]=188;
    isym[260]=189;
    isym[261]=190;
    isym[262]=191;
    isym[263]=-192;
    isym[264]=-193;
    isym[265]=192;
    isym[266]=193;
    isym[267]=196;
    isym[268]=197;
    isym[269]=198;
    isym[270]=199;
    isym[271]=200;
    isym[272]=201;
    isym[273]=-202;
    isym[274]=-203;
    isym[275]=202;
    isym[276]=203;
    isym[277]=223;
    isym[278]=219;
    isym[279]=219;
    isym[280]=226;
    isym[281]=218;
    isym[282]=219;
    isym[283]=220;
    isym[284]=206;
    isym[285]=207;
    isym[286]=208;
    isym[287]=210;
    isym[288]=209;
    isym[289]=211;
    isym[290]=212;
    isym[291]=213;
    isym[292]=216;
    isym[293]=217;
    isym[294]=214;
    isym[295]=215;
    isym[296]=-218;
    isym[297]=-219;
    isym[298]=-220;
    isym[299]=221;
    isym[300]=221;
    isym[301]=230;
    isym[302]=231;
    isym[303]=232;
    isym[304]=233;
    isym[305]=234;
    isym[306]=235;
    isym[307]=-236;
    isym[308]=-237;
    isym[309]=236;
    isym[310]=237;
    isym[312]=253;
    isym[313]=253;
    isym[315]=252;
    isym[316]=253;
    isym[317]=254;
    isym[318]=240;
    isym[319]=241;
    isym[320]=242;
    isym[321]=244;
    isym[322]=243;
    isym[323]=245;
    isym[324]=246;
    isym[325]=247;
    isym[326]=250;
    isym[327]=251;
    isym[328]=248;
    isym[329]=249;
    isym[330]=-252;
    isym[331]=-253;
    isym[332]=-254;
    isym[333]=255;
    isym[334]=255;
    isym[336]=-335;
    isym[339]=-337;
    isym[340]=-338;
    isym[342]=337;
    isym[344]=223;
    isym[345]=-223;
    isym[346]=219;
    isym[347]=226;
    isym[348]=-219;
    isym[349]=-226;
    isym[350]=218;
    isym[351]=219;
    isym[352]=220;
    isym[363]=-353;
    isym[364]=-354;
    isym[365]=-355;
    isym[366]=-356;
    isym[367]=-357;
    isym[368]=-358;
    isym[369]=359;
    isym[370]=360;
    isym[371]=-361;
    isym[372]=-362;
    isym[374]=-373;
    isym[377]=-375;
    isym[378]=-376;
    isym[380]=375;
    isym[392]=-382;
    isym[393]=-383;
    isym[394]=-384;
    isym[395]=-385;
    isym[396]=-386;
    isym[397]=-387;
    isym[398]=388;
    isym[399]=389;
    isym[400]=-390;
    isym[401]=-391;
    isym[406]=405;
    isym[412]=410;
    isym[413]=411;
    isym[416]=335;
    isym[417]=337;
    isym[418]=337;
    isym[419]=338;
    isym[420]=341;
    isym[421]=337;
    isym[422]=343;
    isym[423]=223;
    isym[424]=219;
    isym[425]=219;
    isym[426]=226;
    isym[427]=218;
    isym[428]=219;
    isym[429]=220;
    isym[430]=353;
    isym[431]=354;
    isym[432]=355;
    isym[433]=356;
    isym[434]=357;
    isym[435]=358;
    isym[436]=-361;
    isym[437]=-362;
    isym[438]=359;
    isym[439]=360;
    isym[440]=353;
    isym[441]=354;
    isym[442]=355;
    isym[443]=356;
    isym[444]=357;
    isym[445]=358;
    isym[446]=361;
    isym[447]=362;
    isym[448]=359;
    isym[449]=360;
    isym[450]=373;
    isym[451]=375;
    isym[452]=375;
    isym[453]=376;
    isym[454]=379;
    isym[455]=375;
    isym[456]=381;
    isym[457]=382;
    isym[458]=383;
    isym[459]=384;
    isym[460]=385;
    isym[461]=386;
    isym[462]=387;
    isym[463]=-390;
    isym[464]=-391;
    isym[465]=388;
    isym[466]=389;
    isym[467]=382;
    isym[468]=383;
    isym[469]=384;
    isym[470]=385;
    isym[471]=386;
    isym[472]=387;
    isym[473]=390;
    isym[474]=391;
    isym[475]=388;
    isym[476]=389;
    isym[478]=402;
    isym[479]=403;
    isym[480]=404;
    isym[481]=405;
    isym[482]=405;
    isym[483]=407;
    isym[484]=408;
    isym[485]=409;
    isym[486]=410;
    isym[487]=411;
    isym[488]=410;
    isym[489]=411;
    isym[490]=415;
    isym[491]=414;
    ch[1][0][2]=1.E0;
    ch[2][1][2]=1.E0;
    ch[3][1][3]=1.E0;
    ch[4][1][1]=1.E0;
    ch[5][2][2]=1.15470054E0;
    ch[6][2][3]=1.E0;
    ch[7][2][1]=1.E0;
    ch[8][2][4]=1.E0;
    ch[9][2][0]=1.E0;
    ch[10][0][2]=1.E0;
    ch[10][2][2]=1.33333333E0;
    ch[11][2][3]=1.E0;
    ch[12][2][1]=1.E0;
    ch[13][1][2]=1.15470054E0;
    ch[14][1][3]=1.E0;
    ch[15][1][1]=1.E0;
    ch[18][0][2]=1.E0;
    ch[18][2][2]=-.66666667E0;
    ch[18][2][4]=1.E0;
    ch[19][2][0]=1.E0;
    ch[20][1][3]=-.57735027E0;
    ch[21][1][2]=1.E0;
    ch[23][1][3]=1.E0;
    ch[24][1][1]=1.E0;
    ch[25][0][2]=1.E0;
    ch[25][2][2]=-.66666667E0;
    ch[25][2][4]=-1.E0;
    ch[26][1][1]=-.57735027E0;
    ch[28][1][2]=1.E0;
    ch[29][1][1]=-1.E0;
    ch[30][1][3]=1.E0;
    ch[31][0][2]=1.E0;
    ch[31][2][2]=1.33333333E0;
    ch[32][2][3]=.57735027E0;
    ch[33][2][1]=.57735027E0;
    ch[34][2][4]=-1.15470054E0;
    ch[35][2][0]=-1.15470054E0;
    ch[36][0][2]=1.E0;
    ch[36][2][2]=.66666667E0;
    ch[36][2][4]=1.E0;
    ch[37][2][0]=1.E0;
    ch[38][2][3]=1.E0;
    ch[39][2][1]=1.E0;
    ch[40][0][2]=1.E0;
    ch[40][2][2]=.66666667E0;
    ch[40][2][4]=-1.E0;
    ch[41][2][1]=-1.E0;
    ch[42][2][3]=1.E0;
    ch[43][0][2]=1.E0;
    ch[43][2][2]=-1.33333333E0;
    ch[45][0][2]=1.E0;
    ch[45][2][2]=-1.33333333E0;
    indpp[1][1]=1;
    indpp[2][1]=4;
    indpp[3][1]=5;
    indpp[1][2]=4;
    indpp[2][2]=2;
    indpp[3][2]=6;
    indpp[1][3]=5;
    indpp[2][3]=6;
    indpp[3][3]=3;
    inddp[1][1]=1;
    inddp[2][1]=4;
    inddp[3][1]=7;
    inddp[4][1]=10;
    inddp[5][1]=13;
    inddp[1][2]=2;
    inddp[2][2]=5;
    inddp[3][2]=8;
    inddp[4][2]=11;
    inddp[5][2]=14;
    inddp[1][3]=3;
    inddp[2][3]=6;
    inddp[3][3]=9;
    inddp[4][3]=12;
    inddp[5][3]=15;
    inddd[1][1]=1;
    inddd[2][1]=6;
    inddd[3][1]=7;
    inddd[4][1]=9;
    inddd[5][1]=12;
    inddd[1][2]=6;
    inddd[2][2]=2;
    inddd[3][2]=8;
    inddd[4][2]=10;
    inddd[5][2]=13;
    inddd[1][3]=7;
    inddd[2][3]=8;
    inddd[3][3]=3;
    inddd[4][3]=11;
    inddd[5][3]=14;
    inddd[1][4]=9;
    inddd[2][4]=10;
    inddd[3][4]=11;
    inddd[4][4]=4;
    inddd[5][4]=15;
    inddd[1][5]=12;
    inddd[2][5]=13;
    inddd[3][5]=14;
    inddd[4][5]=15;
    inddd[5][5]=5;
}


// ---- aijm / ddpo (mndod.F90 2125 / 2320) ----
// aij(6,107): [type][atom];  ddp/po: flat l*107+atom (1-based).
namespace mndod_C {
extern std::vector<std::vector<double>> aij;
extern std::vector<int> iii, iiid;
extern std::vector<double> fx;
}
using mndod_C::aij; using mndod_C::iii; using mndod_C::iiid; using mndod_C::fx; using mndod_C::repd;

static double aijl(double z1, double z2, int n1, int n2, int l) {
    return fx[n1+n2+l+1] / std::sqrt(fx[2*n1+1]*fx[2*n2+1])
        * std::pow(2.0*z1/(z1+z2), n1) * std::sqrt(2.0*z1/(z1+z2))
        * std::pow(2.0*z2/(z1+z2), n2) * std::sqrt(2.0*z2/(z1+z2))
        * std::pow(2.0, l) / std::pow(z1+z2, l);
}

void aijm(int ni) {
    double z1 = zs[ni], z2 = zp[ni], z3 = zd[ni];
    int nsp = iii[ni];
    if (ni < 3) return;
    double zz = z1 * z2;
    if (zz < 0.01) return;
    aij[2][ni] = aijl(z1, z2, nsp, nsp, 1);
    aij[3][ni] = aijl(z2, z2, nsp, nsp, 2);
    if (dorbs[ni]) {
        int nd = iiid[ni];
        aij[4][ni] = aijl(z1, z3, nsp, nd, 2);
        aij[5][ni] = aijl(z2, z3, nsp, nd, 1);
        aij[6][ni] = aijl(z3, z3, nd, nd, 2);
    }
}

void ddpo(int ni) {
    double fg = gss[ni];
    po[1][ni] = poij(0, 1.0, fg);
    if (ni >= 3) {
            double d = aij[2][ni] / std::sqrt(12.0);
        fg = hsp[ni];
        ddp[2][ni] = d;
        po[2][ni] = poij(1, d, fg);
        po[7][ni] = po[1][ni];
        d = std::sqrt(aij[3][ni]*0.1);
        fg = 0.5*(gpp[ni]-gp2[ni]);
        ddp[3][ni] = d;
        po[3][ni] = poij(2, d, fg);
        if (dorbs[ni]) {
            double da = std::sqrt(1.0/60.0);
            d = std::sqrt(aij[4][ni]*da);
            fg = repd[19][ni];
            ddp[4][ni] = d;
            po[4][ni] = poij(2, d, fg);
            d = aij[5][ni]/std::sqrt(20.0);
            fg = repd[23][ni] - 1.8*repd[35][ni];
            ddp[5][ni] = d;
            po[5][ni] = poij(1, d, fg);
            fg = 0.2*(repd[29][ni]+2.0*repd[30][ni]+2.0*repd[31][ni]);
            po[8][ni] = poij(0, 1.0, fg);
            d = std::sqrt(aij[6][ni]/14.0);
            fg = repd[44][ni] - (20.0/35.0)*repd[52][ni];
            ddp[6][ni] = d;
            po[6][ni] = poij(2, d, fg);
        }
    }
}
