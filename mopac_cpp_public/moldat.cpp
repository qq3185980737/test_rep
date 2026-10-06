// moldat.cpp — molecular data / electron-counting driver (from moldat.F90).
#include "moldat.h"
#include <cmath>
#include <string>
#include <vector>

#include "big_swap.h"
#include "common_arrays_C.h"
#include "empiri.h"
#include "gmetry.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "parameters_C.h"
#include "reada.h"
#include "refer.h"
#include "set_up_dentate.h"
#include "symtrz.h"
#include "to_screen.h"
#include "vecprt.h"
#include "web_message.h"

namespace moldat_local {
// N_3_present / Si_O_H_present come from molkst_C (F90: use molkst_C)
int old_chrge=0;
int ndorbs=0, nmos=0;
bool method_pm6=false;
std::vector<std::string> txtatm;
std::vector<double> tvec;
int ilog=0;
}
using namespace molkst_C;
using namespace common_arrays_C;
using namespace moldat_local;
using namespace parameters_C;
// msdel comes from molkst_C
extern void setcup();
extern void check_CVS(bool);
extern void check_h(int&);
extern void setup_nhco(int&);
void moldat(int mode){
    int kharge=0,i=0,ndorbs_=0,ia=1,ib=0,nheavy=0,nnull=0,ii=0,k=0,k1=0,j=0,
        ielec=0,ndoubl=0,ne=0,nupp=0,ndown=0,l=0,iminr=0,jminr=0,icount=0,
        ireal=0,jreal=0,ni=0,n1=0,n4=0,n9=0;
    double elecs=0, yy=0, w=0, sum=0, rmin=100.0;
    bool debug=false,exci=false,sing=false,doub=false,trip=false,quar=false,quin=false,
         sext=false,sept=false,octe=false,none=false,birad=false,halfe=false,odd=false;
    std::string num1;
    std::vector<double> rxyz, c(3*(numat+1)+4, 0.0);   // dummy c passed to symtrz; sized safely
    debug = keywrd.find("MOLDAT")!=std::string::npos;
    i = keywrd.find(" CHARGE=");
    if(i!=-1){ kharge=(int)reada(keywrd,i); old_chrge=kharge; } else kharge=0;
    elecs=-kharge; ndorbs_=0;
    if(uss[1]>-1.0){ mopend("THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM"); return; }
    for(i=1;i<=100;++i){
        dorbs[i]=(zd[i]>1e-8);
        if(dorbs[i]) natorb[i]=9;
        else if(zp[i]>1e-20) natorb[i]=4;
        else if(zs[i]>1e-20) natorb[i]=1;
        else natorb[i]=0;
    }
    numat=0; ia=1; ib=0; nheavy=0; nnull=0;
    for(ii=1;ii<=natoms;++ii){
        if(labels[ii]!=99 && labels[ii]!=107){
            ++numat; nat[numat]=labels[ii]; nfirst[numat]=ia; ni=nat[numat];
            elecs+=tore[ni]; ib=ia+natorb[ni]-1;
            if(natorb[ni]==9) ndorbs_+=5;
            nlast[numat]=ib; uspd[ia]=uss[ni];
            if(ia!=ib){
                k=ia+1; k1=ia+3;
                for(j=k;j<=k1;++j) uspd[j]=upp[ni];
                if(ib>ia) ++nheavy; else ++nnull;
                if(k1!=ib){ k=k1+1; for(int kk=k;kk<=ib;++kk) uspd[kk]=udd[ni]; }
            }
        }
        ia=ib+1;
    }
    if(numat==1 && keywrd.find("FORCE")!=std::string::npos){ mopend("A SINGLE ATOM HAS NO VIBRATIONAL MODES"); return; }
    if(mode!=1) refer();
    gmetry(geo,coord);
    if(mode!=1){ empiri(); if(moperr) return; }
    rxyz.assign((numat*(numat+1))/2+1,0.0);
    if(maxtxt==26 && keywrd.find(" GEO-OK")==std::string::npos && keywrd.find(" 0SCF")==std::string::npos){
        l=0;
        for(i=1;i<=numat;++i){
            j=nat[i];
            if(j!=1 && (j<6||j>8)) continue;
            l++;
            if(l==20) break;
        }
        if(l>0){ mopend("Atom name does not match atom label. To suppress this error, add GEO-OK"); return; }
    }
    rmin=100.0; l=0;
    for(i=1;i<=numat;++i) for(j=1;j<=i;++j){
        ++l;
        rxyz[l]=std::sqrt(std::pow(coord[0][i]-coord[0][j],2)+std::pow(coord[1][i]-coord[1][j],2)+std::pow(coord[2][i]-coord[2][j],2));
        if(!(rmin>rxyz[l] && i!=j && (nat[i]<103||nat[j]<103))) continue;
        iminr=i; jminr=j; rmin=rxyz[l];
    }
    setcup(); if(moperr) return;
    if(keywrd.find(" ADD-H")==std::string::npos || keywrd.find(" PDBOUT")!=std::string::npos){
set_up_dentate(); check_CVS(false); check_h(i);
    }
    if(use_ref_geo) big_swap(0,1);
    if(id==0){
        if((keywrd.find(" ADD-H")==std::string::npos && !mozyme) || numat<200){
symtrz(c.data(), pdiag.data(), 1, false); if(moperr) return;
        }
    }
    mol_weight=0; for(i=1;i<=numat;++i) mol_weight+=atmass[i];
    n9=0;n4=0;n1=0;
    for(i=1;i<=natoms;++i){
        j=labels[i];
        if(j>0 && j!=99 && j!=107){
            k=natorb[j];
            if(k==1) ++n1; else if(k==4) ++n4; else if(k==9) ++n9;
        }
    }
    ispd=n9;
    if(id==0) n2elec=2025*n9+100*n4+n1+2025*(n9*(n9-1))/2+450*n9*n4+45*n9*n1+100*(n4*(n4-1))/2+10*n4*n1+(n1*(n1-1))/2+10;
    else      n2elec=2025*n9+100*n4+n1+2025*(n9*(n9+1))/2+450*n9*n4+45*n9*n1+100*(n4*(n4+1))/2+10*n4*n1+(n1*(n1+1))/2+10;
    lm61=45*n9+10*n4+n1;
    norbs=nlast[numat];
    mpack=(norbs*(norbs+1))/2;
    if(keywrd.find(" 0SCF")!=std::string::npos) mpack=1;
    sing=keywrd.find(" SING")!=std::string::npos; doub=keywrd.find(" DOUB")!=std::string::npos;
    trip=keywrd.find(" TRIP")!=std::string::npos; quar=keywrd.find(" QUAR")!=std::string::npos;
    quin=keywrd.find(" QUIN")!=std::string::npos; sext=keywrd.find(" SEXT")!=std::string::npos;
    sept=keywrd.find(" SEPT")!=std::string::npos; octe=keywrd.find(" OCTE")!=std::string::npos;
    none=keywrd.find(" NONE")!=std::string::npos; exci=keywrd.find(" EXCI")!=std::string::npos;
    birad = exci || keywrd.find("BIRAD")!=std::string::npos;
    nalpha=0; nbeta=0;
    nelecs=(int)std::max(elecs,0.0);
    nelecs=std::min(2*norbs,nelecs);
    if(!rhf && !uhf){ if(nelecs%2==0) rhf=true; else uhf=true; }
    odd=(nelecs%2!=0);
    msdel=100000;
    if(sing){ if(odd){mopend("SINGLET SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT");return;} msdel=0; }
    else if(doub){ if(!odd){mopend("DOUBLET SPECIFIED WITH EVEN NUMBER OF ELECTRONS, CORRECT FAULT");return;} msdel=1; }
    else if(trip){ if(odd){mopend("TRIPLET...");return;} msdel=2; }
    else if(quar){ if(!odd){mopend("QUARTET...");return;} msdel=3; }
    else if(quin){ if(odd){mopend("QUINTET...");return;} msdel=4; }
    else if(sext){ if(!odd){mopend("SEXTET...");return;} msdel=5; }
    else if(sept){ if(odd){mopend("SEPTET...");return;} msdel=6; }
    else if(octe){ if(!odd){mopend("OCTET...");return;} msdel=7; }
    else if(none){ if(odd){mopend("NONET...");return;} msdel=8; }
    i=keywrd.find(" MS");
    if(i!=-1){
        if(msdel<100000){ mopend("MS cannot be used with other keywords that define spin"); return; }
        int ims=(int)(2*reada(keywrd,keywrd.find(" MS")));
        if((odd == (ims%2==0)) && keywrd.find(" 0SCF")==std::string::npos){ mopend("Value of MS not consistent with number of electrons"); return; }
        msdel=(int)(2.0*reada(keywrd,keywrd.find(" MS")));
        if((nelecs+msdel)%2==1){ mopend("Impossible value of MS"); return; }
    }
    if(msdel>99999) msdel=0;
    if(uhf){
        nbeta=(nelecs-msdel)/2; nalpha=nelecs-nbeta;
        nopen=0;nclose=0;fract=0;nalpha_open=nalpha;nbeta_open=nbeta;
        i=keywrd.find("OPEN(");
        if(i!=-1){
            std::string sub=keywrd.substr(i,11);
            j=sub.find(',')+i-1;
            nmos=(int)reada(keywrd,j);
            ielec=(int)reada(keywrd,keywrd.find("OPEN(")+5);
            fract=(double)ielec/nmos;
        } else { ielec=0; nmos=0; }
        nalpha-=ielec; nalpha_open=nalpha+nmos; nbeta_open=nbeta;
    } else {
        ielec=0;nmos=0;fract=0;
        if(exci||birad){
            if(odd){mopend("SYSTEM SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT");return;}
            ielec=2;nmos=2;
        } else if((nelecs/2)*2!=nelecs){ ielec=1;nmos=1; }
        i=keywrd.find("OPEN(");
        if(i!=-1){
            std::string sub=keywrd.substr(i,11);
            j=sub.find(',')+i-1;
            nmos=(int)reada(keywrd,j);
            ielec=(int)reada(keywrd,keywrd.find("OPEN(")+5);
        }
        nclose=nelecs/2; nopen=nelecs-nclose*2;
        if(ielec!=0 && norbs>0){
            nclose-=ielec/2; nopen=nmos;
            if(nclose+nopen>norbs){ mopend("NUMBER OF DOUBLY FILLED PLUS PARTLY FILLED LEVELS GREATER THAN TOTAL NUMBER OF ORBITALS"); return; }
            fract=ielec*1.0/nmos;
        }
        nopen=nopen+nclose;
    }
    if(msdel!=0 && !uhf){
        ndoubl=99; nmos=0;
        i=keywrd.find("C.I.=(");
        if(i!=-1){ std::string sub=keywrd.substr(i,11); j=sub.find(',')+i-1; ndoubl=(int)reada(keywrd,j); nmos=(int)reada(keywrd,keywrd.find("C.I.=(")+5); }
        else if(keywrd.find("C.I.=")!=-1){ nmos=(int)reada(keywrd,keywrd.find("C.I.=")+5); }
        else { nmos=nopen-nclose; nmos=std::min(norbs,nmos); }
        int jj;
        if(ndoubl==99) jj=std::max(std::min((nclose+nopen+1)/2-(nmos-1)/2,norbs-nmos+1),1);
        else jj=nclose-ndoubl+1;
        ne=(int)(std::max(0.0,(double)(nclose-jj+1))*2.0+std::max(0.0,(double)(nopen-nclose)*fract)+0.5);
        nupp=(ne+msdel)/2; ndown=ne-nupp;
        if(nmos==0){ nupp=0; ndown=0; }
        if(nupp*ndown<0||nupp>nmos||ndown>nmos||nmos==0){ mopend("SPECIFIED SPIN COMPONENT NOT SPANNED BY ACTIVE SPACE"); return; }
    }
    halfe=(nopen>nclose && std::abs(fract-2.0)>1e-20 && std::abs(fract)>1e-20) || keywrd.find("C.I.")!=std::string::npos;
    yy=(double)kharge/(norbs+1e-10);
    for(i=1;i<=norbs;++i) pdiag[i]=0;
    for(i=1;i<=numat;++i){
        ni=nat[i];
        if(nlast[i]-nfirst[i]==-1) continue;
        l=nfirst[i]-1;
        if(nlast[i]-nfirst[i]==0){ l++; pdiag[l]=tore[ni]-yy; }
        else if(nlast[i]-nfirst[i]==3){ w=tore[ni]*0.25-yy; for(int kk=1;kk<=4;++kk) pdiag[l+kk]=w; }
        else {
            if(ni<21||(ni>30&&ni<39)||(ni>48&&ni<57)){
                w=tore[ni]*0.25-yy; for(int kk=1;kk<=4;++kk) pdiag[l+kk]=w;
                l+=4; for(int kk=1;kk<=5;++kk) pdiag[l+kk]=-yy;
            } else if(ni<99){
                sum=tore[ni]-9*yy;
                l++; pdiag[l]=std::max(0.0,std::min(sum,2.0)); sum-=2.0;
                if(sum>0){
                    l+=3;
                    for(j=1;j<=5;++j){ l++; pdiag[l]=std::max(0.0,std::min(sum*0.2,2.0)); }
                    sum-=10.0;
                    if(sum>0){ l-=8; for(int kk=1;kk<=3;++kk) pdiag[l+kk]=sum/3.0; }
                }
            }
        }
    }
    setup_nhco(ii);
    for(i=1;i<=numat;++i){
        if(nat[i]==7 && nbonds[i]==3){
            j=0;
            if(nat[ibonds[1][i]]==1) j=1;
            if(nat[ibonds[2][i]]==1) j++;
            if(nat[ibonds[3][i]]==1) j++;
            if(j<2) N_3_present=true;
        }
        if(N_3_present) break;
    }
    for(i=1;i<=numat;++i){
        if(nat[i]==14){
            for(j=1;j<=nbonds[i];++j){
                if(nat[ibonds[j][i]]==8){
                    l=ibonds[j][i];
                    for(k=1;k<=nbonds[l];++k) if(nat[ibonds[k][l]]==1) Si_O_H_present=true;
                }
            }
        }
        if(Si_O_H_present) break;
    }
    if(mode!=1 && !debug) return;
}