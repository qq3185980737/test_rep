// H_bonds4.cpp
#include "H_bonds4.h"
#include <algorithm>
#include <cmath>
#include <string>
#include "H_bond_correction_bits.h"
#include "bangle.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using common_arrays_C::coord;
using common_arrays_C::nat;

// energy_corr_h4: PM6-DH4/PM7 hydrogen-bond energy and gradient.
// grad_h4 is a flat (3, numat) workspace, column-major.
double energy_corr_h4(bool l_grad, double* grad_h4) {
    using namespace molkst_C;
    const double HB_R_0 = 1.5, HR_R_CUTOFF = 5.5;
    const double M_PI = 3.14159265358979323846;

    // Method-dependent parameters (all identical across branches here).
    double para_oh_o = 2.32, para_oh_n = 3.10, para_nh_o = 1.07, para_nh_n = 2.01;
    double multiplier_wh_o = 0.42, multiplier_nh4 = 3.61, multiplier_coo = 1.41;

    bool prt = (keywrd.find(" 0SCF ") != std::string::npos ||
                keywrd.find(" PRT ") != std::string::npos) &&
               (keywrd.find(" DISP(") != std::string::npos);

    int n = numat;
    for (int a = 1; a <= n; ++a)
        for (int c = 0; c < 3; ++c) grad_h4[(a - 1) * 3 + c] = 0.0;

    auto gset = [&](int at, double gx, double gy, double gz) {
        grad_h4[(at - 1) * 3 + 0] = gx;
        grad_h4[(at - 1) * 3 + 1] = gy;
        grad_h4[(at - 1) * 3 + 2] = gz;
    };
    auto gadd = [&](int at, double gx, double gy, double gz) {
        grad_h4[(at - 1) * 3 + 0] += gx;
        grad_h4[(at - 1) * 3 + 1] += gy;
        grad_h4[(at - 1) * 3 + 2] += gz;
    };

    double e_corr_sum = 0.0;
    N_Hbonds = 0;
    for (int i = 1; i <= n; ++i) {
        if (nat[i] != 7 && nat[i] != 8) continue;
        for (int j = 1; j <= i - 1; ++j) {
            if (nat[j] != 7 && nat[j] != 8) continue;
            double rda = distance_hb(i, j);
            if (rda < HB_R_0 || rda > HR_R_CUTOFF) continue;

            for (int hi = 1; hi <= n; ++hi) {
                if (nat[hi] != 1) continue;
                double rih = distance_hb(i, hi);
                double rjh = distance_hb(j, hi);
                double angle = 0.0;
                bangle(coord, i, hi, j, angle);
                angle = M_PI - angle;
                if (angle >= M_PI * 0.5) continue;

                int di, ai; double rdh, rah;
                if (rih < rjh) { di = i; ai = j; rdh = rih; rah = rjh; }
                else           { di = j; ai = i; rdh = rjh; rah = rih; }

                // Radial term.
                double e_radial =
                    -0.00303407407407313510*std::pow(rda,7)
                    +0.07357629629627092382*std::pow(rda,6)
                    -0.70087111111082800452*std::pow(rda,5)
                    +3.25309629629461749545*std::pow(rda,4)
                    -7.20687407406838786983*std::pow(rda,3)
                    +5.31754666665572184314*rda*rda
                    +3.40736000001102778967*rda
                    -4.68512000000450434811;
                double drdx=0,drdy=0,drdz=0,radx=0,rady=0,radz=0;
                if (l_grad) {
                    double dradial =
                        -0.02123851851851194655*std::pow(rda,6)
                        +0.44145777777762551519*std::pow(rda,5)
                        -3.50435555555413991158*std::pow(rda,4)
                        +13.01238518517846998179*std::pow(rda,3)
                        -21.62062222220516360949*rda*rda
                        +10.63509333331144368628*rda
                        +3.40736000001102778967;
                    drdx=(coord[0][di]-coord[0][ai])/rda*dradial;
                    drdy=(coord[1][di]-coord[1][ai])/rda*dradial;
                    drdz=(coord[2][di]-coord[2][ai])/rda*dradial;
                    radx=-drdx; rady=-drdy; radz=-drdz;
                }

                // Angular term.
                double a = angle/(M_PI*0.5);
                double ax = -20.0*std::pow(a,7)+70.0*std::pow(a,6)-84.0*std::pow(a,5)+35.0*std::pow(a,4);
                double e_angular = 1.0 - ax*ax;
                double dadx=0,dady=0,dadz=0,daax=0,daay=0,daaz=0,dahx=0,dahy=0,dahz=0;
                if (l_grad) {
                    double xd = (-140.0*std::pow(a,6)+420.0*std::pow(a,5)-420.0*std::pow(a,4)+140.0*std::pow(a,3))/(M_PI*0.5);
                    double d_angular = -2.0*xd*ax;
                    double dd = (coord[0][di]-coord[0][hi])*(coord[0][ai]-coord[0][hi])
                              + (coord[1][di]-coord[1][hi])*(coord[1][ai]-coord[1][hi])
                              + (coord[2][di]-coord[2][hi])*(coord[2][ai]-coord[2][hi]);
                    double xx = -d_angular/std::sqrt(1.0 - (dd*dd)/(rdh*rdh)/(rah*rah));
                    dadx = -xx*((coord[0][ai]-coord[0][hi])/(rdh*rah) - (coord[0][di]-coord[0][hi])*dd/(rdh*rdh*rdh*rah));
                    dady = -xx*((coord[1][ai]-coord[1][hi])/(rdh*rah) - (coord[1][di]-coord[1][hi])*dd/(rdh*rdh*rdh*rah));
                    dadz = -xx*((coord[2][ai]-coord[2][hi])/(rdh*rah) - (coord[2][di]-coord[2][hi])*dd/(rdh*rdh*rdh*rah));
                    daax = -xx*((coord[0][di]-coord[0][hi])/(rdh*rah) - (coord[0][ai]-coord[0][hi])*dd/(rah*rah*rah*rdh));
                    daay = -xx*((coord[1][di]-coord[1][hi])/(rdh*rah) - (coord[1][ai]-coord[1][hi])*dd/(rah*rah*rah*rdh));
                    daaz = -xx*((coord[2][di]-coord[2][hi])/(rdh*rah) - (coord[2][ai]-coord[2][hi])*dd/(rah*rah*rah*rdh));
                    dahx = -dadx-daax; dahy = -dady-daay; dahz = -dadz-daaz;
                }

                // Energy coefficient.
                double e_para;
                if (nat[di]==8 && nat[ai]==8) e_para=para_oh_o;
                else if (nat[di]==8 && nat[ai]==7) e_para=para_oh_n;
                else if (nat[di]==7 && nat[ai]==8) e_para=para_nh_o;
                else e_para=para_nh_n;

                // Bond switching.
                double e_bond_switch=1.0, bsdx=0,bsdy=0,bsdz=0,bsax=0,bsay=0,bsaz=0,bshx=0,bshy=0,bshz=0;
                if (rdh > 1.15) {
                    double rdhs = rdh - 1.15;
                    double ravgs = 0.5*rdh + 0.5*rah - 1.15;
                    double sx = rdhs/ravgs;
                    e_bond_switch = 1.0 - (-20.0*std::pow(sx,7)+70.0*std::pow(sx,6)-84.0*std::pow(sx,5)+35.0*std::pow(sx,4));
                    if (l_grad) {
                        double dbs = -(-140.0*std::pow(sx,6)+420.0*std::pow(sx,5)-420.0*std::pow(sx,4)+140.0*std::pow(sx,3));
                        double xx = dbs/ravgs;
                        double x2 = -0.5*dbs*sx/ravgs;
                        bsdx=(coord[0][di]-coord[0][hi])*(xx+x2)/rdh;
                        bsdy=(coord[1][di]-coord[1][hi])*(xx+x2)/rdh;
                        bsdz=(coord[2][di]-coord[2][hi])*(xx+x2)/rdh;
                        bsax=(coord[0][ai]-coord[0][hi])*x2/rah;
                        bsay=(coord[1][ai]-coord[1][hi])*x2/rah;
                        bsaz=(coord[2][ai]-coord[2][hi])*x2/rah;
                        bshx=-bsdx-bsax; bshy=-bsdy-bsay; bshz=-bsdz-bsaz;
                    }
                }

                // Water scaling.
                double e_scale_w=1.0, sign_wat=1.0;
                if (nat[di]==8 && nat[ai]==8) {
                    double hydrogens=0, others=0;
                    for (int k=1;k<=n;++k) {
                        if (nat[k]==1) hydrogens += cvalence_contribution(di,k);
                        else           others  += cvalence_contribution(di,k);
                    }
                    if (hydrogens >= 1.0) {
                        double slope = multiplier_wh_o - 1.0;
                        double v = hydrogens, fv=0;
                        if (v>1.0 && v<=2.0) { fv=v-1.0; sign_wat=1.0; }
                        if (v>2.0 && v<3.0)  { fv=3.0-v; sign_wat=-1.0; }
                        double fv2 = std::max(1.0-others, 0.0);
                        e_scale_w = 1.0 + slope*fv*fv2;
                    }
                }

                // Charged groups: NR4+ and COO-.
                double e_scale_chd=1.0, e_scale_cha=1.0;
                int cc=-1, o1=ai, o2=-1;
                double f_o1=0,f_o2=0,f_cc=0, cv_o1=0, cv_o2=0, cv_cc=0;
                if (nat[di]==7) {
                    double slope = multiplier_nh4 - 1.0;
                    double v=0;
                    for (int k=1;k<=n;++k) v += cvalence_contribution(di,k);
                    if (v>3.0) v-=3.0; else v=0.0;
                    e_scale_chd = 1.0 + slope*v;
                }
                if (nat[ai]==8) {
                    double slope = multiplier_coo - 1.0;
                    double cdist=9.9e9; cv_o1=0;
                    for (int k=1;k<=n;++k) {
                        double v = cvalence_contribution(o1,k);
                        cv_o1 += v;
                        double s = distance_hb(o1,k);
                        if (v>0.0 && nat[k]==6 && s<cdist) { cdist=s; cc=k; }
                    }
                    if (cc != -1) {
                        double odist=9.9e9; cv_cc=0;
                        for (int k=1;k<=n;++k) {
                            double v = cvalence_contribution(cc,k);
                            cv_cc += v;
                            double s = distance_hb(cc,k);
                            if (v>0.0 && k!=o1 && nat[k]==8 && s<odist) { odist=s; o2=k; }
                        }
                    }
                    if (o2 != -1) {
                        cv_o2=0;
                        for (int k=1;k<=n;++k) cv_o2 += cvalence_contribution(o2,k);
                        f_o1 = std::max(1.0-std::abs(1.0-cv_o1),0.0);
                        f_o2 = std::max(1.0-std::abs(1.0-cv_o2),0.0);
                        f_cc = std::max(1.0-std::abs(3.0-cv_cc),0.0);
                        e_scale_cha = 1.0 + slope*f_o1*f_o2*f_cc;
                    }
                }

                double e_corr = e_para*e_radial*e_angular*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha;
                if (e_corr < -1.0) ++N_Hbonds;
                if (prt) prt_hbonds(di, hi, ai, e_corr);
                e_corr_sum += e_corr;

                double sA=e_para*e_angular*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha;
                double sR=e_para*e_radial*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha;
                double sB=e_para*e_radial*e_angular*e_scale_w*e_scale_chd*e_scale_cha;
                double sC=e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_cha;
                gadd(di, drdx*sA, drdy*sA, drdz*sA);
                gadd(ai, radx*sA, rady*sA, radz*sA);
                gadd(di, dadx*sR, dady*sR, dadz*sR);
                gadd(ai, daax*sR, daay*sR, daaz*sR);
                gadd(hi, dahx*sR, dahy*sR, dahz*sR);
                gadd(di, bsdx*sB, bsdy*sB, bsdz*sB);
                gadd(ai, bsax*sB, bsay*sB, bsaz*sB);
                gadd(hi, bshx*sB, bshy*sB, bshz*sB);

                // Water scaling gradient.
                if (l_grad && std::abs(e_scale_w-1.0)>1e-14) {
                    double slope = multiplier_wh_o - 1.0;
                    for (int k=1;k<=n;++k) {
                        if (k==di) continue;
                        double x = distance_hb(di,k);
                        double xd = (nat[k]==1)
                            ? cvalence_contribution_d(di,k)*sign_wat
                            : cvalence_contribution_d(di,k);
                        double gx=-(coord[0][di]-coord[0][k])*xd*slope/x;
                        double gy=-(coord[1][di]-coord[1][k])*xd*slope/x;
                        double gz=-(coord[2][di]-coord[2][k])*xd*slope/x;
                        gadd(di, -gx*sC, -gy*sC, -gz*sC);
                        gadd(k,   gx*sC,  gy*sC,  gz*sC);
                    }
                }
                // NR4+ gradient.
                if (l_grad && std::abs(e_scale_chd-1.0)>1e-14) {
                    double slope = multiplier_nh4 - 1.0;
                    for (int k=1;k<=n;++k) {
                        if (k==di) continue;
                        double x = distance_hb(di,k);
                        double xd = cvalence_contribution_d(di,k);
                        double gx=-(coord[0][di]-coord[0][k])*xd*slope/x;
                        double gy=-(coord[1][di]-coord[1][k])*xd*slope/x;
                        double gz=-(coord[2][di]-coord[2][k])*xd*slope/x;
                        double s_ = e_para*e_radial*e_angular*e_bond_switch*e_scale_cha*e_scale_w;
                        gadd(di, -gx*s_, -gy*s_, -gz*s_);
                        gadd(k,   gx*s_,  gy*s_,  gz*s_);
                    }
                }
                // COO- gradient.
                if (!l_grad || f_o1*f_o2*f_cc==0.0) continue;
                double slope = multiplier_coo - 1.0;
                double sO = e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w;
                for (int k=1;k<=n;++k) {
                    if (k==o1) continue;
                    double xd = cvalence_contribution_d(o1,k);
                    if (xd==0.0) continue;
                    double x = distance_hb(o1,k);
                    if (cv_o1>1.0) xd=-xd;
                    xd *= f_o2*f_cc;
                    double gx=-(coord[0][o1]-coord[0][k])*xd*slope/x;
                    double gy=-(coord[1][o1]-coord[1][k])*xd*slope/x;
                    double gz=-(coord[2][o1]-coord[2][k])*xd*slope/x;
                    gadd(o1, -gx*sO, -gy*sO, -gz*sO);
                    gadd(k,   gx*sO,  gy*sO,  gz*sO);
                }
                if (o2 != -1)
                for (int k=1;k<=n;++k) {
                    if (k==o2) continue;
                    double xd = cvalence_contribution_d(o2,k);
                    if (xd==0.0) continue;
                    double x = distance_hb(o2,k);
                    if (cv_o2>1.0) xd=-xd;
                    xd *= f_o1*f_cc;
                    double gx=-(coord[0][o2]-coord[0][k])*xd*slope/x;
                    double gy=-(coord[1][o2]-coord[1][k])*xd*slope/x;
                    double gz=-(coord[2][o2]-coord[2][k])*xd*slope/x;
                    gadd(o2, -gx*sO, -gy*sO, -gz*sO);
                    gadd(k,   gx*sO,  gy*sO,  gz*sO);
                }
                for (int k=1;k<=n;++k) {
                    if (k==cc) continue;
                    double xd = cvalence_contribution_d(cc,k);
                    if (xd==0.0) continue;
                    double x = distance_hb(cc,k);
                    if (cv_cc>3.0) xd=-xd;
                    xd *= f_o1*f_o2;
                    double gx=-(coord[0][cc]-coord[0][k])*xd*slope/x;
                    double gy=-(coord[1][cc]-coord[1][k])*xd*slope/x;
                    double gz=-(coord[2][cc]-coord[2][k])*xd*slope/x;
                    gadd(cc, -gx*sO, -gy*sO, -gz*sO);
                    gadd(k,   gx*sO,  gy*sO,  gz*sO);
                }
            }
        }
    }
    return e_corr_sum;
}

static const double covalent_radii[119] = {
    0.0,
    0.37,0.32,1.34,0.9,0.82,0.77,0.75,0.73,0.71,0.69,
    1.54,1.3,1.18,1.11,1.06,1.02,0.99,0.97,
    1.96,1.74,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,
    1.26,1.22,1.19,1.16,1.14,1.1,2.11,1.92,1.62,1.48,1.37,1.45,
    1.56,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33,1.3,
    2.25,1.98,1.69
};

// Smooth covalent-bonding contribution: 1 inside r0, 0 beyond 1.6*r0.
double cvalence_contribution(int atom_a, int atom_b) {
    using namespace common_arrays_C;
    double ri = covalent_radii[nat[atom_a]];
    double rj = covalent_radii[nat[atom_b]];
    double r0 = ri + rj;
    double r1 = r0 * 1.6;
    double r = distance_hb(atom_a, atom_b);
    if (r == 0.0 || r >= r1) return 0.0;
    if (r <= r0) return 1.0;
    double x = (r - r0) / (r1 - r0);
    return 1.0 - (-20.0*x*x*x*x*x*x*x + 70.0*x*x*x*x*x*x - 84.0*x*x*x*x*x + 35.0*x*x*x*x);
}

double cvalence_contribution_d(int atom_a, int atom_b) {
    using namespace common_arrays_C;
    double ri = covalent_radii[nat[atom_a]];
    double rj = covalent_radii[nat[atom_b]];
    double r0 = ri + rj;
    double r1 = r0 * 1.6;
    double r = distance_hb(atom_a, atom_b);
    if (r == 0.0 || r >= r1 || r <= r0) return 0.0;
    double x = (r - r0) / (r1 - r0);
    return -(-140.0*x*x*x*x*x*x + 420.0*x*x*x*x*x - 420.0*x*x*x*x + 140.0*x*x*x) / (r1 - r0);
}

double poly(double r, bool l_grad, double& dpoly) {
    if (r <= 1.0) {
        dpoly = 0.0;
        return 25.46293603147693;
    } else if (r < 1.5) {
        double p = -2714.952351603469651*std::pow(r,5)
                   +17103.650110591705015*std::pow(r,4)
                   -42511.857982217959943*std::pow(r,3)
                   +52063.196799138342612*std::pow(r,2)
                   -31430.658335972289933*r
                   +7516.084696095140316;
        if (l_grad)
            dpoly = -2714.952351603469651*5*std::pow(r,4)
                    +17103.650110591705015*4*std::pow(r,3)
                    -42511.857982217959943*3*std::pow(r,2)
                    +52063.196799138342612*2*r
                    -31430.658335972289933;
        return p;
    }
    double er = std::exp(-1.53965*std::pow(r,1.72905));
    if (l_grad)
        dpoly = -1.53965*1.72905*std::pow(r,0.72905)*118.7326*er;
    return 118.7326*er;
}

double H_bonds4(bool l_grad, double* dxyz) {
    int n = (int)common_arrays_C::nat.size() - 1;
    std::vector<double> grad_h4((n + 1) * 3, 0.0);
    double e = energy_corr_h4(l_grad, grad_h4.data());
    if (l_grad) {
        for (int a = 1; a <= n; ++a)
            for (int c = 0; c < 3; ++c)
                dxyz[c * (n + 1) + a] += grad_h4[(a - 1) * 3 + c];
    }
    return e;
}


// energy_corr_hh_rep: extra H-H repulsion (Rezac & Hobza 2012).
// Iterates over all H-H pairs; poly(r) gives the correction, and the
// gradient accumulates into a column-major (3, numat) workspace.
double energy_corr_hh_rep(bool l_grad, double* dxyz) {
    using namespace molkst_C;
    int n = (int)common_arrays_C::nat.size() - 1;
    std::vector<double> grad_hh((n + 1) * 3, 0.0);
    double e_corr_sum = 0.0;
    for (int i = 1; i <= n; ++i) {
        if (nat[i] != 1) continue;
        for (int j = 1; j <= i - 1; ++j) {
            if (nat[j] != 1) continue;
            double r = distance_hb(i, j);
            double d_rad = 0.0;
            e_corr_sum += poly(r, l_grad, d_rad);
            if (l_grad) {
                double gx = -(coord[0][i] - coord[0][j]) / r * d_rad;
                double gy = -(coord[1][i] - coord[1][j]) / r * d_rad;
                double gz = -(coord[2][i] - coord[2][j]) / r * d_rad;
                grad_hh[(i - 1) * 3 + 0] -= gx; grad_hh[(i - 1) * 3 + 1] -= gy; grad_hh[(i - 1) * 3 + 2] -= gz;
                grad_hh[(j - 1) * 3 + 0] += gx; grad_hh[(j - 1) * 3 + 1] += gy; grad_hh[(j - 1) * 3 + 2] += gz;
            }
        }
    }
    if (l_grad) {
        for (int a = 1; a <= n; ++a)
            for (int c = 0; c < 3; ++c)
                dxyz[c * (n + 1) + a] += grad_hh[(a - 1) * 3 + c];
        // F90 prints "HH REPULSION" + per-atom gradient table to iw under
        // " DERIV"; print-only branch omitted (no behavior change).
    }
    return e_corr_sum;
}
