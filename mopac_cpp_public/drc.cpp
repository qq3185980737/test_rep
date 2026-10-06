// drc.cpp — C++ translation of MOPAC 2016 "drc.F90".
// Dynamic Reaction Coordinate (DRC) / IRC driver: init velocities from
// startv/startk, zero net momentum, integrate Newton's equations over
// time steps (velo1/2/3 series), damp kinetic energy, track gradients and
// errors, call compfg for energies, write restart records on timeout.
// Fortran 1-based indexing preserved.
#include "drc.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <array>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace elemts_C;

extern "C" void to_screen_(const char*);
extern double second(int);
extern double reada(const std::string& string, int istart);
extern void compfg(const std::vector<double>& xparam, bool int_flag,
                   double& escf, bool fulscf, std::vector<double>& grad,
                   bool lgrad);
extern void prtdrc(double deltat, std::vector<double>& xparam,
                   std::vector<double>& ref, double ekin,
                   double& gtot, double& etot, std::vector<double>& velo0,
                   const std::vector<std::array<int, 2>>& mcoprt, int ncoprt,
                   bool parmax);
extern void gmetry(std::vector<std::vector<double>>& geo,
                   std::vector<std::vector<double>>& coord);
extern void l_control(const std::string& keywrd, int len, int n);
extern void mopend(const char* message);
extern void den_in_out(int mode);
extern double ddot(int n, const double* x, int incx, const double* y,
                   int incy);

void drc(std::vector<double>& startv, const std::vector<double>& startk) {
    bool addk = true;
    double ekin = 0.0, escf_old = 10.0, escf_diff = 20.0, elost1 = 0.0,
           etold = 0.0, dlold2 = 0.0;
    std::vector<double> past10(11, 0.0);
    double average_old_hof = 0.0, average_new_hof = 0.0;
    double tnow = second(1), oldtim = second(1), delold = 10.0;
    int percent = 0;
    double gtot = 0.0, damp = 0.99;
    gnorm = 0.0;
    const int iw00 = chanel_C::iw0;
    chanel_C::iw0 = -1;
    int iloop = 1;
    double escf_min = 1e20, start_hof_local = 1e20;
    int n_escf = 0;
    l_control("LDRC_FIRST", 10, 1);
    double minstep;
    if (nopen != nclose && keywrd.find(" IRC") != std::string::npos)
        minstep = 5e-16;
    else
        minstep = 1e-16;

    grad.assign(3 * numat + 1, 0.0);
    errfn.assign(3 * numat + 1, 0.0);

    double accu;
    if (keywrd.find(" PREC") != std::string::npos)
        accu = 0.25;
    else
        accu = 1.0;
    double stepx, stepxx;
    if (keywrd.find(" DRC") == std::string::npos) {
        if (keywrd.find(" X-PRIORITY=") != std::string::npos)
            stepx = reada(keywrd, static_cast<int>(keywrd.find("X-PRIO")) + 5);
        else
            stepx = 0.05;
        stepx *= 0.2;
        stepxx = 0.01;
    } else {
        stepx = 0.0;
        stepxx = 0.0;
    }
    int ilim;
    if (keywrd.find(" EPS") != std::string::npos) {
        ilim = 300;
        accu *= 10;
    } else {
        ilim = 30;
    }
    double gnlim = 1.0;
    past10[5] = 100.0;
    const bool debug = keywrd.find(" DEBUG") != std::string::npos;
    int i = static_cast<int>(keywrd.find("GNORM"));
    if (i != -1) gnlim = reada(keywrd, i);
    int n_min;
    if (gnlim < 0.9)
        n_min = 80;
    else
        n_min = 40;

    int iskin = 0;
    double addonk = 0.0;
    if (keywrd.find("KINE") != std::string::npos) {
        iskin = 1;
        addonk = reada(keywrd, static_cast<int>(keywrd.find("KINE")));
        if (addonk < 0.0) {
            for (auto& v : startv) v = -v;
        }
        std::printf("\n\n          EXCESS KINETIC ENERGY ENTERED INTO SYSTEM =%12.6f\n",
                    std::abs(addonk));
        if (addonk < 0.0) {
            addonk = -addonk;
            std::printf("          KINETIC ENERGY SUPPLIED WAS NEGATIVE, SO INITIAL VELOCITY IS REVERSED\n");
        }
    }
    const bool velred = keywrd.find("VELO") != std::string::npos;
    if (ddot(3 * numat, &startv[0], 1, &startv[0], 1) > 0.001) {
        std::printf("          INITIAL VELOCITY IN DRC (Angstroms/Femtosecond)\n");
        for (int iv = 1; iv <= numat * 3; iv += 3)
            std::printf("%16.5f%16.5f%16.5f\n", startv[iv - 1], startv[iv],
                        startv[iv + 1]);
        for (int iv = 0; iv < numat * 3; ++iv) startv[iv] = -startv[iv];
    }
    bool let = velred && keywrd.find(" IRC") == std::string::npos;
    if (keywrd.find(" SYMM") != std::string::npos) ndep = 0;

    if ((keywrd.find(" XYZ") == std::string::npos &&
         keywrd.find("VELO") == std::string::npos) ||
        keywrd.find("IRC=") != std::string::npos) {
        gmetry(geo, coord);
        int l0 = 0;
        for (int j = 1; j <= 3; ++j)
            for (int iat = 1; iat <= numat; ++iat) {
                geo[j][iat] = coord[j - 1][iat];
                coord[j - 1][iat] = 0.0;
            }
        na.assign(numat + 1, 0);
    }
    for (int iat = 1; iat <= numat; ++iat) {
        if (atmass[iat] >= 0.1) continue;
        std::printf(" ATOMIC MASS OF ATOM%3d TOO SMALL\n", iat);
        return;
    }

    double half;
    bool parmax = false;
    std::vector<std::vector<int>> mcoprt(3, std::vector<int>(3 * numat + 1));
    int ncoprt = 0;
    if (keywrd.find(" DRC") != std::string::npos) {
        parmax = (loc[1][1] != 0);
        if (parmax) {
            for (int iv = 1; iv <= nvar; ++iv) {
                mcoprt[1][iv] = loc[1][iv];
                mcoprt[2][iv] = loc[2][iv];
            }
            ncoprt = nvar;
        }
    }
    ncoprt = 0;

    int l = 0;
    std::vector<std::vector<double>> georef(4, std::vector<double>(numat + 1));
    loc[1].assign(3 * numat + 1, 0);
    loc[2].assign(3 * numat + 1, 0);
    xparam.assign(3 * numat + 1, 0.0);
    for (int iat = 1; iat <= numat; ++iat) {
        for (int j = 1; j <= 3; ++j) {
            loc[1][l + j] = iat;
            loc[2][l + j] = j;
            georef[j][iat] = geo[j][iat];
            xparam[l + j] = geo[j][iat];
        }
        l += 3;
    }
    nvar = numat * 3;

    if (keywrd.find("DRC=") != std::string::npos) {
        half = reada(keywrd, static_cast<int>(keywrd.find("DRC=")));
        std::printf("\n\n          DAMPING FACTOR FOR KINETIC ENERGY =%12.6f\n", half);
    } else if (keywrd.find(" DRC") == std::string::npos) {
        half = 0.0;
        if (addonk > 1e-4) {
            for (int iv = 0; iv < 3 * numat; ++iv) startv[iv] *= addonk;
        }
    } else {
        half = 1e6;
    }
    const bool letot_0 = keywrd.find("IRC=") == std::string::npos && !let;
    bool letot = letot_0;
    half = std::copysign(std::max(1e-6, std::abs(half)), half);

    double deltat = 1e-16;
    if (keywrd.find("IRC=") != std::string::npos)
        deltat = std::min(1e-13, deltat * std::pow(numat * 0.25, 4));
    double quadr = 1.0, etot = 0.0;
    escf = 0.0;
    double const_ = 1.0;

    std::vector<double> velo0(3 * numat + 1, 0.0), velo1(3 * numat + 1, 0.0),
        velo2(3 * numat + 1, 0.0), velo3(3 * numat + 1, 0.0),
        gerror(3 * numat + 1, 0.0), grold(3 * numat + 1, 0.0),
        grold2(3 * numat + 1, 0.0);
    const bool restart = keywrd.find("RESTART") != std::string::npos &&
                         keywrd.find("IRC=") == std::string::npos;

    if (restart) {
        std::ifstream rf(chanel_C::restart_fn,
                         std::ios::binary | std::ios::in | std::ios::ate);
        if (!rf.is_open()) {
            std::printf(" Restart file either does not exist or is not available for reading\n");
            mopend("Restart file either does not exist or is not available for reading");
            return;
        }
        rf.seekg(0);
        int io_stat = 0;
        for (int iv = 1; iv <= nvar; ++iv)
            if (!rf.read(reinterpret_cast<char*>(&xparam[iv]), sizeof(double)))
                io_stat = -1;
        for (int iv = 1; iv <= nvar; ++iv)
            if (!rf.read(reinterpret_cast<char*>(&velo0[iv]), sizeof(double)))
                io_stat = -1;
        for (int iv = 1; iv <= nvar; ++iv)
            if (!rf.read(reinterpret_cast<char*>(&grad[iv]), sizeof(double)))
                io_stat = -1;
        for (int iv = 1; iv <= nvar; ++iv)
            if (!rf.read(reinterpret_cast<char*>(&grold[iv]), sizeof(double)))
                io_stat = -1;
        for (int iv = 1; iv <= nvar; ++iv)
            if (!rf.read(reinterpret_cast<char*>(&grold2[iv]), sizeof(double)))
                io_stat = -1;
        double delold_r, deltat_r, dlold2_r;
        int i_r;
        char letot_c;
        double elost1_r, gtot_r;
        if (!rf.read(reinterpret_cast<char*>(&etot), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&escf), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&ekin), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&delold_r), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&deltat_r), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&dlold2_r), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&i_r), sizeof(int)) ||
            !rf.read(reinterpret_cast<char*>(&gnorm), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&letot_c), sizeof(char)) ||
            !rf.read(reinterpret_cast<char*>(&elost1_r), sizeof(double)) ||
            !rf.read(reinterpret_cast<char*>(&gtot_r), sizeof(double)))
            io_stat = -1;
        rf.close();
        if (io_stat != 0) {
            std::printf(" Restart file is corrupt\n");
            mopend("Restart file is corrupt");
            return;
        }
        delold = delold_r;
        deltat = deltat_r;
        dlold2 = dlold2_r;
        iloop = i_r - 1;
        letot = letot_c != 0;
        elost1 = elost1_r;
        gtot = gtot_r;
        std::printf("\n\n          CALCULATION RESTARTED, CURRENT KINETIC ENERGY=%10.5f\n\n",
                    ekin);
        goto restart_point;
    } else {
        velo0.assign(3 * numat + 1, 0.0);
        grold2.assign(3 * numat + 1, 0.0);
        grold.assign(3 * numat + 1, 0.0);
        grad.assign(3 * numat + 1, 0.0);
        if (keywrd.find("IRC=") != std::string::npos || velred) {
            int k;
            if (keywrd.find("IRC=") != std::string::npos)
                k = static_cast<int>(std::lround(
                    reada(keywrd, static_cast<int>(keywrd.find("IRC=")))));
            else
                k = 1;
            double one;
            if (k < 0) {
                k = -k;
                one = -1.0;
            } else {
                one = 1.0;
            }
            const int kl = (k - 1) * nvar;
            double summ = 0.0;
            std::vector<double> velo1c(4, 0.0);
            double summas = 0.0;
            i = 0;
            for (int iat = 1; iat <= numat; ++iat) {
                const double ams = atmass[iat];
                summas += ams;
                velo0[i + 1] = startv[kl + i] * one;
                velo1c[1] += velo0[i + 1] * ams;
                velo0[i + 2] = startv[kl + i + 1] * one;
                velo1c[2] += velo0[i + 2] * ams;
                velo0[i + 3] = startv[kl + i + 2] * one;
                velo1c[3] += velo0[i + 3] * ams;
                i += 3;
            }
            for (int j = 1; j <= 3; ++j) velo1c[j] = -velo1c[j] / summas;
            i = 0;
            if (addonk > 1e-5 || !velred) {
                for (int iat = 1; iat <= numat; ++iat) {
                    const double ams = atmass[iat];
                    for (int i1 = 1; i1 <= 3; ++i1) {
                        ++i;
                        velo0[i] += velo1c[i1];
                        summ += velo0[i] * velo0[i] * ams;
                    }
                }
            } else {
                for (int iat = 1; iat <= numat; ++iat) {
                    const double ams = atmass[iat];
                    for (int i1 = 1; i1 <= 3; ++i1) {
                        ++i;
                        summ += velo0[i] * velo0[i] * ams;
                    }
                }
            }
            if (addonk < 1e-5 && velred)
                addonk = 0.5 * summ / 4.184e10;
            if (addonk < 1e-5 && !velred) {
                if (std::abs(half) > 1e-3 && startk[k] > 105.0) {
                    std::printf(" BY DEFAULT, ONE QUANTUM OF ENERGY, EQUIVALENT TO%10.3f CM(-1)\n WILL BE USED TO START THE DRC\n",
                                startk[k]);
                    addonk = startk[k] * 2.8585086e-3;
                    std::printf(" THIS REPRESENTS AN ENERGY OF%7.2f KCALS/MOLE\n",
                                addonk);
                } else if (std::abs(half) > 1e-3) {
                    std::printf(" THE VIBRATIONAL FREQUENCY (%9.2fCM(-1)) IS TOO SMALL\n FOR ONE QUANTUM TO BE USED\n",
                                startk[k]);
                    std::printf(" INSTEAD 0.3KCAL/MOLE WILL BE USED TO START THE IRC\n");
                    addonk = 0.3;
                } else {
                    addonk = 0.3;
                }
            }
            if (summ < 1e-4) {
                std::printf(" SYSTEM IS APPARENTLY NOT MOVING!\n");
                return;
            }
            if (half > 0.1 && half < 10000.0)
                addonk *= (1.0 + 0.06972 / half);
            if (half < -0.1 && half > -10000.0)
                addonk *= (1.0 + 0.06886 / half);
            const double summ_s = std::sqrt(addonk / (0.5 * summ / 4.184e10));
            addk = false;
            if (summ_s > 1e-10) {
                for (int iv = 1; iv <= nvar; ++iv) velo0[iv] *= summ_s;
                if (half > 1e-3) addonk = 0.0;
            }
        }
    }
restart_point:;
    int bigcycles;
    if (keywrd.find(" BIGCYCLES") != std::string::npos)
        bigcycles = static_cast<int>(std::lround(
            reada(keywrd, static_cast<int>(keywrd.find(" BIGCYCLES"))))) * 4;
    else
        bigcycles = -1;
    int maxcyc, jloop_lim;
    if (keywrd.find(" CYCLES") != std::string::npos) {
        maxcyc = 1000000;
        jloop_lim = static_cast<int>(std::lround(
            reada(keywrd, static_cast<int>(keywrd.find(" CYCLES")))));
    } else {
        maxcyc = 4999;
        jloop_lim = 100000;
    }
    int iupper = iloop + maxcyc;
    const int ilp = iloop;
    double one = 0.0;
    if (restart) one = 1.0;
    gerror.assign(3 * numat + 1, 0.0);

    for (iloop = ilp; iloop <= iupper; ++iloop) {
        const double error = etot - (ekin + escf);
        if (iloop > 2) {
            quadr = 1.0 + error / (ekin * const_ + 0.001) * 0.5;
            quadr = std::min(1.3, std::max(0.8, quadr));
        } else {
            quadr = 1.0;
        }
        if ((let || ekin > 0.2) && addk) {
            etot += addonk;
            addk = false;
            addonk = 0.0;
        }
        const_ = std::max(1e-36, std::pow(0.5, deltat * 1e15 / half));
        const_ = std::sqrt(const_);
        double velvec = 0.0;
        ekin = 0.0;
        const double delta1 = delold + dlold2;
        double elost = 0.0;
        for (int iv = 1; iv <= nvar; ++iv) startv[iv - 1] = xparam[iv];
        if (iloop > 3) {
            for (int iv = 1; iv <= nvar; ++iv) {
                velo1[iv] = 1.0 / atmass[loc[1][iv]] * grad[iv];
                velo3[iv] = 2.0 / atmass[loc[1][iv]] *
                    (delta1 * (grold[iv] - grad[iv]) -
                     delold * (grold2[iv] - grad[iv])) /
                    (delta1 * (delold * delold * 1e30) -
                     delold * (delta1 * delta1 * 1e30));
                velo2[iv] = 1.0 / atmass[loc[1][iv]] *
                    (grad[iv] - grold[iv] -
                     0.5 * velo3[iv] * (1e30 * delold * delold)) /
                    (delold * 1e15);
            }
            for (;;) {
                for (int iv = 1; iv <= nvar; ++iv) {
                    startv[iv - 1] =
                        xparam[iv] -
                        1e8 * (deltat * velo0[iv] * one +
                               0.5 * deltat * deltat * velo1[iv] +
                               0.16666 * (deltat * deltat * 1e15) * deltat *
                                   velo2[iv] +
                               0.0416666 * deltat * deltat *
                                   (1e30 * deltat * deltat) * velo3[iv]);
                }
                double sum = 0.0;
                for (int iv = 1; iv <= nvar; ++iv)
                    sum += (xparam[iv] - startv[iv - 1]) *
                           (xparam[iv] - startv[iv - 1]);
                if (stepxx < 1e-5) break;   // drc.f90 L485
                sum = stepxx / std::sqrt(sum);
                if (sum > 0.8 && sum < 1.25) break;
                deltat *= std::max(0.99, std::min(1.01, sum));
            }
            for (int iv = 1; iv <= nvar; ++iv) xparam[iv] = startv[iv - 1];
            double sum = 0.0;
            for (int iat = 1; iat <= numat; ++iat)
                for (int j = 1; j <= 3; ++j)
                    sum += (geo[j][iat] - georef[j][iat]) *
                           (geo[j][iat] - georef[j][iat]);
            if (sum > 0.01) gnorm = 4.0;
            for (int iv = 1; iv <= nvar; ++iv) {
                velvec += velo0[iv] * velo0[iv];
                velo0[iv] = velo0[iv] + deltat * velo1[iv] +
                            0.5 * deltat * deltat * velo2[iv] * 1e15 +
                            0.166666 * deltat * (1e30 * deltat * deltat) *
                                velo3[iv];
                if (let || gnorm > 3.0) {
                    let = true;
                    elost += velo0[iv] * velo0[iv] * atmass[loc[1][iv]] *
                             (1 - const_ * const_);
                    velo0[iv] *= const_ * quadr;
                }
                ekin += velo0[iv] * velo0[iv] * atmass[loc[1][iv]];
            }
        } else {
            for (int iv = 1; iv <= nvar; ++iv) {
                velo1[iv] = 1.0 / atmass[loc[1][iv]] * grad[iv];
                velo2[iv] = 1.0 / atmass[loc[1][iv]] *
                            (grad[iv] - grold[iv]) / (1e15 * delold);
                velo3[iv] = 0.0;
                xparam[iv] =
                    xparam[iv] -
                    1e8 * (deltat * velo0[iv] * one +
                           0.5 * deltat * deltat * velo1[iv] +
                           0.16666 * (deltat * deltat * 1e15) * deltat *
                               velo2[iv] +
                           0.0416666 * deltat * deltat *
                               (1e30 * deltat * deltat) * velo3[iv]);
                velvec += velo0[iv] * velo0[iv];
                velo0[iv] = velo0[iv] + deltat * velo1[iv] +
                            0.5 * deltat * deltat * velo2[iv] * 1e15 +
                            0.166666 * deltat * (1e30 * deltat * deltat) *
                                velo3[iv];
                if (let || gnorm > 3.0) {
                    let = true;
                    elost += velo0[iv] * velo0[iv] * atmass[loc[1][iv]] *
                             (1 - const_ * const_);
                    velo0[iv] *= const_ * quadr;
                }
                ekin += velo0[iv] * velo0[iv] * atmass[loc[1][iv]];
            }
        }
        one = 1.0;
        if (let || gnorm > 3.0) {
            if (!letot) {
                if (std::abs(half) < 1e-3) {
                    etot = escf + elost1;
                    addonk = 0.0;
                    elost1 = 0.0;
                    elost = 0.0;
                    deltat = 5e-16;
                } else if (iskin == 0) {
                    etot -= addonk;
                }
            }
            letot = true;
        }
        ekin = 0.5 * ekin / 4.184e10;
        if (letot && std::abs(half) > 0.00001)
            etot = etot - ekin / (const_ * const_) + ekin;
        elost1 += 0.5 * elost / 4.184e10;

        grold2.assign(grold.begin(), grold.end());
        grold.assign(grad.begin(), grad.end());
        grad.assign(3 * numat + 1, 0.0);

        if (std::abs(average_new_hof) > 1e-20)
            average_old_hof =
                damp * average_old_hof + (1.0 - damp) * escf;
        compfg(xparam, true, escf, true, grad, true);
        if (std::abs(average_new_hof) < 1e-20) {
            average_old_hof = escf + 5.0;
            average_new_hof = escf;
        }
        if (moperr) return;
        if (debug) {
            std::printf("\n\nPoint calculated:    %8d\n", iloop);
            std::printf("Heat of formation:    %13.5f\n", escf);
            if (start_hof_local > 1e19) start_hof_local = escf;
            std::printf("Calculated energy change:%10.5f\n", elost1);
            std::printf("Predicted energy change:%11.5f\n", etot - escf);
            std::printf("Cumulative error:       %11.5f\n", escf + elost1 - etot);
            std::printf("'Time' interval (fs):%13.4f\n", deltat * 1e15);
            std::printf("\nGeometry supplied to COMPFG in DRC\n");
            for (int iat = 1; iat <= numat; ++iat)
                std::printf(" %2s%15.4f%17.4f%17.4f\n",
                            elemnt[nat[iat]].c_str(),
                            xparam[(iat - 1) * 3 + 1],
                            xparam[(iat - 1) * 3 + 2],
                            xparam[(iat - 1) * 3 + 3]);
            std::printf("\nForces (gradients) acting on atoms\n");
            for (int iat = 1; iat <= numat; ++iat)
                std::printf(" %2s%15.4f%17.4f%17.4f\n",
                            elemnt[nat[iat]].c_str(),
                            grad[(iat - 1) * 3 + 1],
                            grad[(iat - 1) * 3 + 2],
                            grad[(iat - 1) * 3 + 3]);
            if (velvec > 1.0) {
                std::printf("\nVelocity of atoms, in cm/sec\n");
                for (int iat = 1; iat <= numat; ++iat)
                    std::printf(" %2s%15.4f%17.4f%17.4f\n",
                                elemnt[nat[iat]].c_str(),
                                -velo0[(iat - 1) * 3 + 1],
                                -velo0[(iat - 1) * 3 + 2],
                                -velo0[(iat - 1) * 3 + 3]);
            } else {
                std::printf("\nAcceleration of atoms, in 10^(-14)cm/(sec*sec)\n");
                for (int iat = 1; iat <= numat; ++iat)
                    std::printf(" %2s%15.4f%17.4f%17.4f\n",
                                elemnt[nat[iat]].c_str(),
                                -1e-14 * velo1[(iat - 1) * 3 + 1],
                                -1e-14 * velo1[(iat - 1) * 3 + 2],
                                -1e-14 * velo1[(iat - 1) * 3 + 3]);
            }
        }
        average_new_hof = damp * average_new_hof + (1.0 - damp) * escf;
        if (iloop > 2) {
            gnorm = 0.0;
            for (int iv = 1; iv <= nvar; iv += 3) {
                const double sum =
                    std::sqrt(ddot(3, &grad[iv], 1, &grad[iv], 1) /
                              (ddot(3, &velo0[iv], 1, &velo0[iv], 1) + 1e-20));
                for (int jv = 0; jv < 3; ++jv) {
                    gerror[iv + jv] += grad[iv + jv] + velo0[iv + jv] * sum;
                }
            }
            gnorm = std::sqrt(ddot(nvar, &gerror[1], 1, &gerror[1], 1));
            gtot = gnorm;
        }
        gnorm = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
        for (int iv = 1; iv <= nvar; ++iv) grad[iv] *= 4.184e18;
        if (iloop == 1) {
            for (int iv = 1; iv <= nvar; ++iv) grold[iv] = grad[iv];
        }
        dlold2 = delold;
        delold = deltat;
        double sum = 0.0;
        for (int iv = 1; iv <= nvar; ++iv)
            sum += std::pow((grad[iv] - grold[iv]) / 4.184e18, 2);
        if (std::abs(half) < 0.001) {
            deltat *= std::pow(
                std::min(2.0, 2e-4 * accu /
                                  (std::abs(escf + elost1 - etold) + 1e-20)),
                0.25);
            etold = escf + elost1;
            if (iloop > ilim && addonk < 1e-5 &&
                average_old_hof - average_new_hof < 0.0001) {
                chanel_C::iw0 = iw00;
                std::printf("\n\n IRC CALCULATION COMPLETE \n");
                return;
            }
            if (escf < escf_min || iloop < 10) {
                escf_min = escf;
                n_escf = 0;
            } else {
                ++n_escf;
                if (n_escf > n_min) {
                    std::printf("\n\n          POTENTIAL ENERGY HAS STOPPED DROPPING\n\n");
                    chanel_C::iw0 = iw00;
                    return;
                }
            }
        } else {
            deltat *= std::min(1.05, 10.0 * accu / (sum + 1e-4));
            deltat = std::min(deltat, 3e-15 * accu);
            past10[10] = gnorm;
            sum = 0.0;
            for (int ip = 1; ip <= 9; ++ip) {
                sum += std::abs(past10[ip] - past10[ip + 1]);
                past10[ip] = past10[ip + 1];
            }
            if (sum < gnlim) {
                std::printf("\n\n GRADIENT CONSTANT AND SMALL -- ASSUME ALL MOTION STOPPED\n");
                std::printf(" TO CONTINUE, USE KEYWORD 'GNORM=0'\n");
                chanel_C::iw0 = iw00;
                return;
            }
            deltat = std::min(deltat, 2e-15);
        }
        deltat = std::max(minstep, deltat);
        std::vector<double> ref1d(3 * numat + 1, 0.0);
        for (int i = 1; i <= numat; ++i)
            for (int d = 0; d < 3; ++d) ref1d[(i - 1) * 3 + d + 1] = georef[d + 1][i];
        std::vector<std::array<int, 2>> mcoprt2(ncoprt + 1);
        for (int k = 1; k <= ncoprt; ++k) {
            mcoprt2[k][0] = mcoprt[1][k]; mcoprt2[k][1] = mcoprt[2][k];
        }
        if (std::abs(half) < 0.00001) {
            prtdrc(deltat, xparam, ref1d, elost1, gtot, etot, velo0, mcoprt2,
                   ncoprt, parmax);
        } else {
            double dummy = 0.0;
            prtdrc(deltat, xparam, ref1d, ekin, dummy, etot, velo0, mcoprt2,
                   ncoprt, parmax);
            if (iloop > 10 && (escf - escf_old) / escf_diff < 0.0) {
                --bigcycles;
                if (bigcycles == 0) iupper = iloop;
            }
            escf_diff = escf - escf_old;
            escf_old = escf;
        }
        tnow = second(2);
        const double tcycle = tnow - oldtim;
        oldtim = tnow;
        tleft -= tcycle;
        if (iw00 > -1) {
            i = static_cast<int>(std::lround(100.0 * iloop / iupper));
            if (i != percent) {
                percent = i;
                char linebuf[32];
                std::snprintf(linebuf, sizeof(linebuf), "%4d%% of DRC/IRC done", i);
                line = linebuf;
                to_screen_(line.c_str());
            }
        }
        if (itemp_1 < jloop_lim && iloop != iupper && tleft >= 3 * tcycle)
            continue;
        if (tleft < 3 * tcycle) {
            std::ofstream rf(chanel_C::restart_fn,
                             std::ios::binary | std::ios::out |
                                 std::ios::trunc);
            if (rf.is_open()) {
                for (int iv = 1; iv <= nvar; ++iv)
                    rf.write(reinterpret_cast<const char*>(&xparam[iv]),
                             sizeof(double));
                for (int iv = 1; iv <= nvar; ++iv)
                    rf.write(reinterpret_cast<const char*>(&velo0[iv]),
                             sizeof(double));
                for (int iv = 1; iv <= nvar; ++iv)
                    rf.write(reinterpret_cast<const char*>(&grad[iv]),
                             sizeof(double));
                for (int iv = 1; iv <= nvar; ++iv)
                    rf.write(reinterpret_cast<const char*>(&grold[iv]),
                             sizeof(double));
                for (int iv = 1; iv <= nvar; ++iv)
                    rf.write(reinterpret_cast<const char*>(&grold2[iv]),
                             sizeof(double));
                const int i_next = iloop + 1;
                rf.write(reinterpret_cast<const char*>(&etot), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&escf), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&ekin), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&delold), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&deltat), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&dlold2), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&i_next), sizeof(int));
                rf.write(reinterpret_cast<const char*>(&gnorm), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&letot), sizeof(char));
                rf.write(reinterpret_cast<const char*>(&elost1), sizeof(double));
                rf.write(reinterpret_cast<const char*>(&gtot), sizeof(double));
                rf.close();
            }
            escf = -1e9;
            prtdrc(deltat, xparam, ref1d, ekin, elost, etot, velo0, mcoprt2,
                   ncoprt, parmax);
            den_in_out(1);
        }
        for (int j = 1; j <= 3; ++j)
            for (int iat = 1; iat <= numat; ++iat) {
                geo[j][iat] = coord[j - 1][iat];
                coord[j - 1][iat] = 0.0;
            }
        na.assign(numat + 1, 0);
        if ((itemp_1 < jloop_lim || iloop == iupper) && bigcycles < 0) {
            std::printf("\n          NUMBER OF CYCLES EXCEEDED, RESTART FILE WRITTEN\n");
            std::printf(" Number of cycles allowed in this run =%7d\n", maxcyc + 1);
            std::printf(" To increase the number of cycles, use keyword \"CYCLES=n\",\n");
            std::printf(" where \"n\" is the number of cycles you want to run\n");
        } else {
            std::printf("\n          RUNNING OUT OF TIME, RESTART FILE WRITTEN\n");
        }
        chanel_C::iw0 = iw00;
        return;
    }
    chanel_C::iw0 = iw00;
}
