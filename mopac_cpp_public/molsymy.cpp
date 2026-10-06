// molsymy.cpp — C++ translation of MOPAC point-group symmetry driver.
#include "molsymy.h"
#include "symmetry_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "bldsym.h"
#include "rsp.h"
#include "symopr.h"
#include "rotmol.h"
#include "geout.h"
#include "bangle.h"
#include "mopend.h"
#include "mult33.h"
#include <cmath>
#include <vector>
#include <string>
using namespace common_arrays_C;
using namespace molkst_C;
using namespace symmetry_C;

void chi(double toler, double* coord, int ioper, int& iqual);
void orient(int numat, double* coord, double* r);
void plato(double* coord, double* r, bool& cubic);
void cartab();

void molsym(double* coord, int& ierror, double* r) {
    auto C = [&](int i, int a) -> double& { return coord[(i - 1) + (a - 1) * 3]; }; // coord(3,numat) col-major
    auto R = [&](int i, int j) -> double& { return r[(i - 1) + (j - 1) * 3]; };

    int icyc[7] = {0,0,0,0,0,0,0};
    int i, j, k, l, ij, iturn, iqual = 0, kndex, icheck, naxes, iz, ix, iy;
    double f[7] = {0,0,0,0,0,0,0};
    double ew[4] = {0,0,0,0}, help[4] = {0,0,0,0};
    double rhelp[4][4];
    double shift[4] = {0,0,0,0};
    double toler, wmol, sum, rxy, tole, distxy, rmin, sina, cosa, theta, total;
    bool linear, cubic, axis, sphere, debug, reorie, iscube;

    debug = keywrd.find(" MOLSYM") != std::string::npos;
    reorie = keywrd.find(" NOREOR") == std::string::npos;
    toler = 0.1; ierror = 0; name = "????";
    for (i = 1; i <= 3; ++i) { for (j = 1; j <= 3; ++j) R(i,j) = 0; R(i,i) = 1; }
    // cub reset
    for (i = 1; i <= 3; ++i) { for (j = 1; j <= 3; ++j) cub[i][j] = 0; cub[i][i] = 1; }

    for (i = 1; i <= 18; ++i) { bldsym(i, i); ielem[i] = 0; }
    ielem[19] = 0; ielem[20] = 0;

    wmol = 0; for (i=1;i<=3;++i) shift[i]=0;
    for (i = 1; i <= numat; ++i) {
        wmol += nat[i];
        for (j = 1; j <= 3; ++j) shift[j] += nat[i] * C(j, i);
    }
    ij = 0;
    for (i = 1; i <= 3; ++i) {
        shift[i] /= wmol;
        for (j = 1; j <= numat; ++j) C(i, j) -= shift[i];
        for (j = 1; j <= i; ++j) {
            ++ij; f[ij] = ij * 1e-8;
            for (k = 1; k <= numat; ++k) f[ij] += nat[k] * C(i,k) * C(j,k);
        }
    }
    iscube = true;

    // --- outer do loop (lines 88..222) ---
    do {
        if (!reorie) {
            sum = 0; for (i=1;i<=5;++i) sum += std::fabs(f[i]);
            linear = sum < 0.01;
            sphere = linear && std::fabs(f[6]) < 0.01;
            cubic = std::fabs(f[1]-f[3])<0.01 && std::fabs(f[1]-f[6])<0.01 &&
                    std::fabs(f[2])+std::fabs(f[4])+std::fabs(f[5])<0.01;
            for (i=1;i<=3;++i){for(j=1;j<=3;++j)R(i,j)=0;R(i,i)=1;}
        } else {
            rsp(f, 3, ew, r);
            sum = R(1,1)*(R(2,2)*R(3,3)-R(3,2)*R(2,3)) +
                  R(1,2)*(R(2,3)*R(3,1)-R(2,1)*R(3,3)) +
                  R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1));
            if (sum > 1.0) {
                sum = 1.0;
                for (j=1;j<=3;++j){ if(R(j,j)>=sum) continue; sum=R(j,j); i=j; }
                for (k=1;k<=3;++k) R(k,i) = -R(k,i);
            }
            R(1,3)=R(2,1)*R(3,2)-R(3,1)*R(2,2);
            R(2,3)=R(3,1)*R(1,2)-R(1,1)*R(3,2);
            R(3,3)=R(1,1)*R(2,2)-R(2,1)*R(1,2);
            linear = ew[1] < 1e-2;               // F90 ew(2)
            sphere = ew[2] < 1e-2;               // F90 ew(3)
            cubic = ew[2]-ew[0] < 5e-3*std::max(ew[2],40.0);   // F90 ew(3)-ew(1)
        }
        if (sphere) {
            ielem[7]=1; ielem[8]=1; ielem[10]=1; ielem[20]=1;
            for(i=1;i<=3;++i)for(j=1;j<=3;++j)cub[i][j]=0;
            cub[1][2]=1; cub[2][3]=1; cub[3][1]=1;
            goto L270;
        } else if (linear) {
            symopr(numat, coord, 1, r);
            ielem[20] = 1;
            goto L250;
        }
        if (reorie) {
            if (!cubic && ew[3]-ew[2] < 1e-2*ew[3]) {
                for (i=1;i<=3;++i){ rxy=-R(i,1); R(i,1)=R(i,3); R(i,3)=rxy; }
                rxy=ew[1]; ew[1]=ew[3]; ew[3]=rxy;
            }
            axis = std::fabs(ew[0]-ew[1]) < 0.01*ew[1];
        } else {
            axis = std::fabs(f[1]-f[3]) < 0.01*f[6];
        }
        symopr(numat, coord, 1, r);
        if (!cubic) goto L1100;
        ielem[19]=1; plato(coord, r, cubic);
        if (moperr) return;
        if (!iscube || cubic) goto L1000;
        for (i=1;i<=6;++i) f[i]=0;
        ij=0;
        for (i=1;i<=3;++i) for (j=1;j<=i;++j) {
            ++ij; f[ij]=ij*1e-8;
            for (k=1;k<=numat;++k){
                sum=0;
                for (l=1;l<=numat;++l)
                    sum += std::sqrt(std::pow(C(1,k)-C(1,l),2)+std::pow(C(2,k)-C(2,l),2)+std::pow(C(3,k)-C(3,l),2));
                f[ij]+=sum*C(i,k)*C(j,k);
            }
        }
        f[1]=f[1]-0.8*f[1]; f[3]=f[3]-0.8*f[1]; f[6]=f[6]-0.8*f[1];
        ielem[19]=0; iscube=false;
    } while (false);
    // fall-through after the outer loop (non-sphere/linear pseudo-cubic exit)
    ielem[7]=1; ielem[8]=1; ielem[10]=1; ielem[20]=1;
    for(i=1;i<=3;++i)for(j=1;j<=3;++j)cub[i][j]=0;
    cub[1][2]=1; cub[2][3]=1; cub[3][1]=1;
    goto L270;

L1000:
    if (moperr) return;
L1100:
    if (axis) {
L130:
        iturn = 7; j = 0;
        for (i = 8; i <= 18; ++i) {
            if (i < 14) tole = toler*9.0/std::pow(i-5,2);
            else tole = toler*16.0/std::pow(2*i-24,2);
            chi(tole, coord, i, iqual);
            if (ielem[i]!=1 || i>=14) continue;
            if (iturn > 9) { toler*=0.5; goto L130; }
            iturn = i;
        }
        if (ielem[14]+ielem[15]+ielem[17]>1 || ielem[15]+ielem[16]+ielem[17]>1 ||
            ielem[16]+ielem[17]+ielem[18]>1) { toler*=0.5; goto L130; }
        iturn -= 5;
        // find adjacent equivalent pair off the Z axis to define X axis
        kndex = 0;
        for (i = 1; i <= numat; ++i) {
            distxy = C(1,i)*C(1,i)+C(2,i)*C(2,i);
            if (distxy < toler) continue;
            rmin = 1000.0;
            for (j = i+1; j <= numat; ++j) {
                if (std::fabs(std::fabs(C(3,i))-std::fabs(C(3,j))) > 0.2) continue;
                rxy = C(1,j)*C(1,j)+C(2,j)*C(2,j);
                if (std::fabs(rxy-distxy)>toler || nat[i]!=nat[j]) continue;
                if (nbonds[i]!=nbonds[j]) continue;
                std::vector<bool> used(nbonds[i]+1, false);
                for (k=1;k<=nbonds[i];++k){
                    bool matched=false;
                    for (l=1;l<=nbonds[j];++l)
                        if (nat[ibonds[l][j]]==nat[ibonds[k][i]] && !used[l]) { used[l]=true; matched=true; break; }
                    (void)matched;
                }
                l=0; for(k=1;k<=nbonds[i];++k) if(used[k]) ++l;
                if (l < nbonds[i]) continue;
                rxy = std::pow(C(1,i)-C(1,j),2)+std::pow(C(2,i)-C(2,j),2);
                if (rxy > rmin) continue;
                kndex = j; rmin = rxy;
            }
            if (kndex >= 1) break;
        }
        if (kndex < 1) { axis=false; goto L190; }
        help[1]=C(1,i)+C(1,kndex); help[2]=C(2,i)+C(2,kndex);
        distxy = std::sqrt(help[1]*help[1]+help[2]*help[2]);
        sina=help[2]/distxy; cosa=help[1]/distxy;
        rotmol(numat, coord, sina, cosa, 1, 2, r);
        chi(toler, coord, 5, iqual);
        if (ielem[5] != 1) {
            chi(toler, coord, 1, iqual);
            if (ielem[1] != 0) {
                theta = 1.5707963268/(double)iturn;
                sina=std::sin(theta); cosa=std::cos(theta);
                icheck=0;
L180:
                rotmol(numat, coord, sina, cosa, 1, 2, r);
                if (icheck > 0) { axis = (iturn != 2); goto L190; }
                chi(toler, coord, 5, iqual);
                if (ielem[5] > 0) goto L190;
                icheck=1; sina=-sina; goto L180;
            }
        }
    }
L190:
    if (cubic) orient(numat, coord, r);
    if (!axis) {
        toler = 0.2;
        for (i=1;i<=6;++i){ chi(toler, coord, i, iqual); icyc[i]=(1+iqual)*ielem[i]; }
        if (reorie) {
            naxes = ielem[1]+ielem[2]+ielem[3];
            if (naxes <= 1) {
                iz=1; if(ielem[1]==1) goto L220;
                iz=2; if(ielem[2]==1) goto L220;
                iz=3; if(ielem[3]==1) goto L220;
                if(icyc[5]>icyc[4]) iz=2;
                if(icyc[6]>icyc[7-iz]) iz=1;
            } else {
                iz=1;
                if(icyc[2]>icyc[1]) iz=2;
                if(icyc[3]>icyc[iz]) iz=3;
            }
L220:
            icyc[7-iz] = -1;
            ix=1;
            if(icyc[5]>icyc[6]) ix=2;
            if(icyc[4]>icyc[7-ix]) ix=3;
            iy = 6 - ix - iz;
            for(k=1;k<=3;++k){ rhelp[k][1]=R(k,ix); rhelp[k][2]=R(k,iy); }
            rhelp[1][3]=R(2,ix)*R(3,iy)-R(3,ix)*R(2,iy);
            rhelp[2][3]=R(3,ix)*R(1,iy)-R(1,ix)*R(3,iy);
            rhelp[3][3]=R(1,ix)*R(2,iy)-R(2,ix)*R(1,iy);
            symopr(numat, coord, -1, r);
            for(k=1;k<=3;++k)for(j=1;j<=3;++j)R(k,j)=rhelp[k][j];
            symopr(numat, coord, 1, r);
        }
    }
L250:
    for (i=1;i<=7;++i) chi(toler, coord, i, iqual);
L270:
    symopr(numat, coord, -1, r);
    total = ew[1]+ew[2]+ew[3];
    for (i=1;i<=3;++i) ew[i] = total - ew[i];
    cartab();
}

// chi: test whether symmetry operation ioper leaves the geometry invariant.
void chi(double toler, double* coord, int ioper, int& iqual) {
    auto C = [&](int i, int a) -> double& { return coord[(i-1) + (a-1)*3]; };
    int iresul = 1;
    iqual = 0;
    for (int i = 1; i <= numat; ++i) {
        double help1 = C(1,i)*elem[1][1][ioper] + C(2,i)*elem[1][2][ioper] + C(3,i)*elem[1][3][ioper];
        double help2 = C(1,i)*elem[2][1][ioper] + C(2,i)*elem[2][2][ioper] + C(3,i)*elem[2][3][ioper];
        double help3 = C(1,i)*elem[3][1][ioper] + C(2,i)*elem[3][2][ioper] + C(3,i)*elem[3][3][ioper];
        bool found = false;
        for (int j = 1; j <= numat; ++j) {
            if (nat[i] != nat[j]) continue;
            if (std::fabs(C(1,j)-help1) > toler) continue;
            if (std::fabs(C(2,j)-help2) > toler) continue;
            if (std::fabs(C(3,j)-help3) > toler) continue;
            jelem[ioper][i] = j;
            if (i == j) ++iqual;
            found = true; break;
        }
        if (!found) iresul = 0;
    }
    ielem[ioper] = iresul;
}

// makopr: build the full set of symmetry operations for the identified group.
void makopr(int numat, double* coord, int& ierror, double* r) {
    symopr(numat, coord, 1, r);
    if (nclass < 2) return;
    for (int i = 2; i <= nclass; ++i) bldsym(jy[i], i);
    double toler = 0.2;
    int iqual = 0;
    for (int i = 2; i <= nclass; ++i) {
        chi(toler, coord, i, iqual);
        if (ielem[i] >= 1) continue;
        ierror = 5;
    }
    symopr(numat, coord, -1, r);
}

// orient: for cubic systems, rotate to align principal axes; probe S4/C5 axes.
void orient(int numat, double* coord, double* r) {
    double toler = 0.1;
    double wink[3] = {0, 0.955316618125, 0.65235813978437};
    double wink2 = 0.0, sina, cosa, sinb, cosb;
    if (ielem[8] >= 1) {
        int iqual = 0;
        for (int i = 1; i <= 2; ++i) {
            int jota = 18 - 4*i;
            wink2 = wink[i];
            sina = std::sin(wink2); cosa = std::cos(wink2);
            rotmol(numat, coord, sina, cosa, 1, 3, r);
            chi(toler, coord, jota, iqual);
            if (ielem[jota] > 0) break;
            if (i == 1) { chi(toler, coord, 3, iqual); if (ielem[3] == 1) break; }
            wink2 = -wink2;
            sinb = std::sin(2.0*wink2); cosb = std::cos(2.0*wink2);
            rotmol(numat, coord, sinb, cosb, 1, 3, r);
            chi(toler, coord, jota, iqual);
            if (ielem[jota] > 0) break;
            if (i == 1) { chi(toler, coord, 3, iqual); if (ielem[3] == 1) break; }
            rotmol(numat, coord, sina, cosa, 1, 3, r);
        }
        chi(toler, coord, 9, iqual);
        if (ielem[10] > 0) chi(toler, coord, 17, iqual);
    } else {
        wink2 = -wink[1];
        if (ielem[10] > 0) wink2 = -wink[2];
        sina = -std::sin(wink2); cosa = std::cos(wink2);
        rotmol(numat, coord, sina, cosa, 1, 3, r);
        int iqual = 0;
        chi(toler, coord, 8, iqual);
        rotmol(numat, coord, -sina, cosa, 1, 3, r);
        if (ielem[8] <= 0) {
            if (ielem[9] <= 0) wink2 = -wink2;
            else rotmol(numat, coord, 0.707106781186, 0.707106781186, 1, 2, r);
        }
    }
    int j = 0; for (int i = 1; i <= 17; ++i) j += ielem[i];
    if (j == 2 && ielem[1]+ielem[8] == 2) return;
    cub[1][1] = std::cos(wink2);
    cub[3][3] = cub[1][1];
    cub[1][3] = std::sin(wink2);
    cub[3][1] = -cub[1][3];
    double cub3[9];
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j) cub3[(j - 1) * 3 + (i - 1)] = cub[i][j];
    mult33(cub3, 8);
    mult33(cub3, 15);
    int iqual = 0;
    chi(toler, coord, 8, iqual);
    j = iqual;
    chi(toler, coord, 15, iqual);
    if (j + iqual == 0) return;
}

// plato: identify Platonic/polygon axis and build the orienting unitary matrix r.
void plato(double* coord, double* r, bool& cubic) {
    auto C = [&](int i, int a) -> double& { return coord[(i-1)+(a-1)*3]; };
    auto R = [&](int i, int j) -> double& { return r[(i-1)+(j-1)*3]; };
    std::vector<int> near(numat+1, 0);
    int ipoly[4] = {0,0,0,0};
    std::vector<std::vector<double>> xyz(4, std::vector<double>(4, 0.0));
    double toler = 0.1, xmin, dist, r2j, angle=0, sum, buff, buff1, rmin, ymin, dist1, dist2;
    double allr[4] = {0,0,0,0};
    int i, l, j, ii, k, i1, j1, m, i3, i2, ligand_type=0, store_ligand=0, j2, k1, k2;

    symopr(numat, coord, -1, r);
    rmin = 0.1; ymin = 0; k = 100;
    // shell loop
    bool found_shell = false;
    while (!found_shell) {
        xmin = 0; ligand_type = 0;
        for (i = 1; i <= numat; ++i) {
            xmin = C(1,i)*C(1,i)+C(2,i)*C(2,i)+C(3,i)*C(3,i);
            if (xmin >= rmin) { ligand_type = nat[i]; break; }
        }
        if (xmin < rmin) { xmin = ymin; found_shell = true; continue; }
        l = 0;
        for (i = 1; i <= numat; ++i) {
            dist = C(1,i)*C(1,i)+C(2,i)*C(2,i)+C(3,i)*C(3,i);
            if (std::fabs(dist-xmin) < toler && nat[i]==ligand_type) ++l;
        }
        if (l < k) { k = l; ymin = xmin; store_ligand = ligand_type; }
        if (l < 5) { found_shell = true; continue; }
        rmin = xmin + 2.0*toler;
    }
    ligand_type = store_ligand; l = 0;
    for (i = 1; i <= numat; ++i) {
        dist = C(1,i)*C(1,i)+C(2,i)*C(2,i)+C(3,i)*C(3,i);
        if (dist >= toler && dist <= xmin+toler)
            if (std::fabs(dist-xmin) < toler && nat[i]==ligand_type) { ++l; near[l]=i; }
    }
    if (l == 0) return;
    if (l == 1) { cubic = false; return; }
    if (l==4||l==6||l==8||l==12||l==20) {
        xmin = 10000.0; j = near[1];
        for (ii = 2; ii <= l; ++ii) {
            i = near[ii];
            dist = std::pow(C(1,j)-C(1,i),2)+std::pow(C(2,j)-C(2,i),2)+std::pow(C(3,j)-C(3,i),2);
            if (dist < xmin) xmin = dist;
        }
        j1 = near[1]; j2 = near[2]; k1 = 0; k2 = 0;
        for (ii = 1; ii <= l; ++ii) {
            i = near[ii];
            dist1 = std::pow(C(1,j1)-C(1,i),2)+std::pow(C(2,j1)-C(2,i),2)+std::pow(C(3,j1)-C(3,i),2);
            dist2 = std::pow(C(1,j2)-C(1,i),2)+std::pow(C(2,j2)-C(2,i),2)+std::pow(C(3,j2)-C(3,i),2);
            if (std::fabs(dist2-xmin) < toler) ++k2;
            if (std::fabs(dist1-xmin) < toler) ++k1;
        }
        if (k1==k2) k = k1;
        else if (k1>2||k2>2) { cubic=false; return; }
        else k = 0;
        if ((l==4&&k==3)||(l==6&&k==4)||(l==8&&k==3)||(l==12&&k==5)||(l==20&&k==3)) {
            i = near[l];
            double dd = std::sqrt(std::pow(C(1,i),2)+std::pow(C(2,i),2)+std::pow(C(3,i),2));
            R(1,3)=C(1,i)/dd; R(2,3)=C(2,i)/dd; R(3,3)=C(3,i)/dd;
            goto L1100;
        }
    }
    i1 = near[1];
    for (i = 1; i <= 3; ++i) {
        xmin = 10000.0; ipoly[i] = 0;
        for (j1 = 2; j1 <= l; ++j1) {
            j = near[j1];
            bool skip = false;
            for (m = 1; m <= i-1; ++m) if (ipoly[m]==j) { skip=true; break; }
            if (skip) continue;
            r2j = std::pow(C(1,i1)-C(1,j),2)+std::pow(C(2,i1)-C(2,j),2)+std::pow(C(3,i1)-C(3,j),2);
            if (xmin > r2j) {
                xmin = r2j; ipoly[i]=j; allr[i]=r2j;
                if (xmin < 0.01) {
                geout(6);
                mopend("Geometry is severely faulty");
                return;
            }
            }
        }
    }
    for (i = 1; i <= 3; ++i) allr[i] = std::sqrt(allr[i]);
    if (ipoly[1]*ipoly[2]*ipoly[3]==0) { cubic=false; return; }
    {
        int ii=0, jj=0;
        bool found = false;
        for (i = 1; i <= 2 && !found; ++i)
            for (j = i+1; j <= 3; ++j)
                if (std::fabs(allr[i]-allr[j]) < 0.01) { ii=i; jj=j; found=true; break; }
        if (!found) { cubic=false; return; }
        i = ii; j = jj;
    }
    for (k = 1; k <= 3; ++k) {
        xyz[k][3] = C(k, ipoly[i]);
        xyz[k][2] = C(k, ipoly[j]);
        xyz[k][1] = C(k, i1);
    }
    i3 = ipoly[i]; i2 = ipoly[j];
    ipoly[1]=i3; ipoly[2]=i2;
    bangle(xyz, 3, 1, 2, angle);
    if (std::fabs(angle-1.0472) < 0.1) {
        for (i=1;i<=3;++i) R(i,3)=C(i,i1)+C(i,i2)+C(i,i3);
    } else if (std::fabs(angle-1.5707963) < 0.1) {
        for (i=1;i<=3;++i) R(i,3)=C(i,i2)+C(i,i3);
    } else if (std::fabs(angle-1.885) < 0.1) {
        for (i=1;i<=3;++i) R(i,3)=1.6180341*(C(i,i2)+C(i,i3))-C(i,i1);
    } else { cubic=false; return; }
    L1100:
    sum = std::sqrt(R(1,3)*R(1,3)+R(2,3)*R(2,3)+R(3,3)*R(3,3));
    for (i=1;i<=3;++i) R(i,3) /= sum;
    buff = std::sqrt(R(1,3)*R(1,3)+R(2,3)*R(2,3));
    buff1 = std::sqrt(R(1,3)*R(1,3)+R(3,3)*R(3,3));
    if (buff <= buff1) { R(1,1)=R(3,3)/buff1; R(2,1)=0; R(3,1)=-R(1,3)/buff1; }
    else { R(1,1)=R(2,3)/buff; R(2,1)=-R(1,3)/buff; R(3,1)=0; }
    R(1,2)=R(2,3)*R(3,1)-R(2,1)*R(3,3);
    R(2,2)=R(3,3)*R(1,1)-R(3,1)*R(1,3);
    R(3,2)=R(1,3)*R(2,1)-R(1,1)*R(2,3);
    symopr(numat, coord, 1, r);
}

// cartab data arrays come from symmetry_C.h / symmetry_tables.cpp.
namespace cartab_local {
    int nope[58], nrep[58];
    std::string classn[21] = {"","C2(x)","C2(y)","C2(z)","Sigma(XY)","Sigma(v)","Sigma(d)","Inversion"," C3"," C4"," C5"," C6"," C7"," C8"," S4"," S6"," S8"," S10"," S12"," ?? ","C(Inf)"};
    bool debug=false, first=true, large=false;
    int icalcn=0, prnt=1;
}

bool symdec(int n1, const int* ielem) {
    int in1 = n1;
    for (int i = 1; i <= 20; ++i) {
        int ibin = in1 % 2;
        if (ielem[i] != 1 && ibin == 1) return false;
        in1 = in1 / 2;
    }
    return true;
}

void cartab() {
    using namespace cartab_local;
    if (numcal != icalcn) {
        icalcn = numcal;
        debug = keywrd.find("CARTAB") != std::string::npos;
        large = keywrd.find("LARGE") != std::string::npos && debug;
    }
    large = large && (prnt == 1);
    if (large) prnt = 0;
    if (first) {
        first = false;
        nope[1] = 1; nrep[1] = 1;
        for (int i = 2; i <= 57; ++i) {
            nope[i] = nope[i-1] + nallop[nope[i-1]] + 4;
            nrep[i] = nrep[i-1] + nallop[nope[i-1]+1] + 1;
        }
        int j = 1;
        for (int i = 1; i <= ntbs; ++i) { int k = ntab[i]; ntab[i] = j; j += k; }
    }
    int igroup = 1;
    for (int g = 57; g >= 1; --g) {
        if (!symdec(nallop[nope[g]+3], ielem)) continue;
        igroup = g; break;
    }
    if (keywrd.find(" NOSYM") != std::string::npos) igroup = 1;
    int istart = nope[igroup];
    int igp = nrep[igroup];
    name = allrep[igp];
    nclass = nallop[istart] + 1;
    nirred = nallop[istart+1];
    for (int k = 1; k <= nclass; ++k) group[(k-1)*20 + 0] = 1.0;
    for (int k = 2; k <= nclass; ++k) jy[k] = nallop[istart+4+k-2];
    for (int i = 1; i <= nirred; ++i) jx[i] = allrep[igp + i];
    istart = ntab[nallop[istart+2]] - 1;
    for (int i = 2; i <= nirred; ++i) {
        for (int k = 1; k <= nclass; ++k) {
            ++istart;
            double buff = nallg[istart];
            if (buff >= 10.0 && (name != "R3" || (k != 1 && k != nclass))) {
                int nzz = nallg[istart], nz = nzz/10;
                double fz = nz, fn = nzz - 10*nz;
                buff = 2.0*std::cos(6.283185307179*fn/fz);
            }
            group[(k-1)*20 + (i-1)] = buff;
        }
    }
    jy[1] = 0;
}