// lbfgs_line_search.cpp — C++ translation of L-BFGS-B line-search:
// dcsrch (Moré-Thuente), dcstep (step update), dpmeps (machine epsilon).
#include <cmath>`n#include <limits>
#include <string>

double dpmeps() { return std::numeric_limits<double>::epsilon(); }

void dcstep(double& stx, double& fx, double& dx, double& sty, double& fy,
            double& dy, double& stp, double fp, double dp, bool& brackt,
            double stpmin, double stpmax) {
    const double zero = 0.0, p66 = 0.66, two = 2.0, three = 3.0;
    double sgnd = dp * (dx / std::fabs(dx));
    if (std::fabs(stp - stx) < 1e-5) stp = stp + 1e-5;
    double theta, s, gamma, p, q, r, stpc, stpq, stpf;
    if (fp > fx) {
        theta = three * (fx - fp) / (stp - stx) + dx + dp;
        s = std::max({std::fabs(theta), std::fabs(dx), std::fabs(dp)});
        gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if (stp < stx) gamma = -gamma;
        p = (gamma - dx) + theta;
        q = ((gamma - dx) + gamma) + dp;
        r = p / q;
        stpc = stx + r * (stp - stx);
        stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx);
        stpf = (std::fabs(stpc - stx) < std::fabs(stpq - stx))
                   ? stpc
                   : stpc + (stpq - stpc) / two;
        brackt = true;
    } else if (sgnd < zero) {
        theta = three * (fx - fp) / (stp - stx) + dx + dp;
        s = std::max({std::fabs(theta), std::fabs(dx), std::fabs(dp)});
        gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if (stp > stx) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dx;
        r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        stpf = (std::fabs(stpc - stp) > std::fabs(stpq - stp)) ? stpc : stpq;
        brackt = true;
    } else if (std::fabs(dp) < std::fabs(dx)) {
        theta = three * (fx - fp) / (stp - stx) + dx + dp;
        s = std::max({std::fabs(theta), std::fabs(dx), std::fabs(dp)});
        gamma = s * std::sqrt(std::max(0.0, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
        if (stp > stx) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = (gamma + (dx - dp)) + gamma;
        r = p / q;
        if (r < zero && gamma != zero) stpc = stp + r * (stx - stp);
        else if (stp > stx) stpc = stpmax;
        else stpc = stpmin;
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (brackt) {
            stpf = (std::fabs(stpc - stp) < std::fabs(stpq - stp)) ? stpc : stpq;
            stpf = (stp > stx) ? std::min(stp + p66 * (sty - stp), stpf)
                               : std::max(stp + p66 * (sty - stp), stpf);
        } else {
            stpf = (std::fabs(stpc - stp) > std::fabs(stpq - stp)) ? stpc : stpq;
            stpf = std::min(stpmax, stpf);
            stpf = std::max(stpmin, stpf);
        }
    } else if (brackt) {
        theta = three * (fp - fy) / (sty - stp) + dy + dp;
        s = std::max({std::fabs(theta), std::fabs(dy), std::fabs(dp)});
        gamma = s * std::sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
        if (stp > sty) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dy;
        r = p / q;
        stpc = stp + r * (sty - stp);
        stpf = stpc;
    } else if (stp > stx) {
        stpf = stpmax;
    } else {
        stpf = stpmin;
    }

    if (fp > fx) {
        sty = stp; fy = fp; dy = dp;
    } else {
        if (sgnd < zero) { sty = stx; fy = fx; dy = dx; }
        stx = stp; fx = fp; dx = dp;
    }
    stp = stpf;
}

void dcsrch(double f, double g, double& stp, double ftol, double gtol,
             double xtol, double stpmin, double stpmax, std::string& task,
             int* isave, double* dsave) {
    const double zero = 0.0, p5 = 0.5, p66 = 0.66, xtrapl = 1.1, xtrapu = 4.0;
    bool brackt;
    int stage;
    double finit, fm, ftest, fx, fxm, fy, fym, ginit, gm, gtest, gx, gxm, gy, gym;
    double stmax, stmin, stx, sty, width, width1;

    if (task.rfind("START", 0) == 0) {
        if (stp < stpmin) task = "ERROR: STP .LT. STPMIN";
        if (stp > stpmax) task = "ERROR: STP .GT. STPMAX";
        if (g >= zero) task = "ERROR: INITIAL G .GE. ZERO";
        if (ftol < zero) task = "ERROR: FTOL .LT. ZERO";
        if (gtol < zero) task = "ERROR: GTOL .LT. ZERO";
        if (xtol < zero) task = "ERROR: XTOL .LT. ZERO";
        if (stpmin < zero) task = "ERROR: STPMIN .LT. ZERO";
        if (stpmax < stpmin) task = "ERROR: STPMAX .LT. STPMIN";
        if (task.rfind("ERROR", 0) == 0) return;
        brackt = false; stage = 1;
        finit = f; ginit = g;
        gtest = ftol * ginit;
        width = stpmax - stpmin;
        width1 = width / p5;
        stx = 1e-5; fx = finit; gx = ginit;
        sty = 1e-4; fy = finit; gy = ginit;
        stmin = zero; stmax = stp + xtrapu * stp;
        task = "FG";
    } else {
        brackt = (isave[1] == 1);
        stage = isave[2];
        ginit = dsave[1]; gtest = dsave[2]; gx = dsave[3]; gy = dsave[4];
        finit = dsave[5]; fx = dsave[6]; fy = dsave[7];
        stx = dsave[8]; sty = dsave[9]; stmin = dsave[10]; stmax = dsave[11];
        width = dsave[12]; width1 = dsave[13];
        ftest = finit + stp * gtest;
        if (stage == 1 && f <= ftest && g >= zero) stage = 2;
        if (brackt && (stp <= stmin || stp >= stmax))
            task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS";
        if (brackt && stmax - stmin <= xtol * stmax)
            task = "WARNING: XTOL TEST SATISFIED";
        if (std::fabs(stp - stpmax) < 1e-20 && f <= ftest && g <= gtest)
            task = "WARNING: STP = STPMAX";
        if (std::fabs(stp - stpmin) < 1e-20 && (f > ftest || g >= gtest))
            task = "WARNING: STP = STPMIN";
        if (f <= ftest && std::fabs(g) <= gtol * (-ginit))
            task = "CONVERGENCE";
        if (task.rfind("WARN", 0) != 0 && task.rfind("CONV", 0) != 0) {
            if (stage == 1 && f <= fx && f > ftest) {
                fm = f - stp * gtest;
                fxm = fx - stx * gtest;
                fym = fy - sty * gtest;
                gm = g - gtest;
                gxm = gx - gtest;
                gym = gy - gtest;
                dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
                fx = fxm + stx * gtest;
                fy = fym + sty * gtest;
                gx = gxm + gtest;
                gy = gym + gtest;
            } else {
                dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
            }
            if (brackt) {
                if (std::fabs(sty - stx) >= p66 * width1)
                    stp = stx + p5 * (sty - stx);
                width1 = width;
                width = std::fabs(sty - stx);
            }
            if (brackt) {
                stmin = std::min(stx, sty);
                stmax = std::max(stx, sty);
            } else {
                stmin = stp + xtrapl * (stp - stx);
                stmax = stp + xtrapu * (stp - stx);
            }
            stp = std::max(stp, stpmin);
            stp = std::min(stp, stpmax);
            if ((brackt && (stp <= stmin || stp >= stmax)) ||
                (brackt && stmax - stmin <= xtol * stmax))
                stp = stx;
            task = "FG";
        }
    }
    isave[1] = brackt ? 1 : 0;
    isave[2] = stage;
    dsave[1] = ginit; dsave[2] = gtest; dsave[3] = gx; dsave[4] = gy;
    dsave[5] = finit; dsave[6] = fx; dsave[7] = fy;
    dsave[8] = stx; dsave[9] = sty; dsave[10] = stmin; dsave[11] = stmax;
    dsave[12] = width; dsave[13] = width1;
}

// dpofa: Cholesky factorization of symmetric positive-definite matrix.
// a is stored row-major [lda][n+1] (1-based); overwritten with L.
double ddot(int n, const double* dx, int incx, const double* dy, int incy);
void daxpy(int n, double sa, const double* sx, int incx, double* sy, int incy);

void dpofa(double* a, int lda, int n, int& info) {
    auto A = [&](int r, int c) -> double& { return a[(r) * lda + (c)]; };
    for (int j = 1; j <= n; ++j) {
        info = j;
        double s = 0.0;
        int jm1 = j - 1;
        for (int k = 1; k <= jm1; ++k) {
            double t = A(k, j) - ddot(k - 1, &A(1, k), lda, &A(1, j), lda);
            t /= A(k, k);
            A(k, j) = t;
            s += t * t;
        }
        s = A(j, j) - s;
        if (s <= 0.0) return;
        A(j, j) = std::sqrt(s);
    }
    info = 0;
}

// dtrsl: solve triangular system T*x = b (or transpose), job selects form.
void dtrsl(double* t, int ldt, int n, double* b, int job, int& info) {
    auto T = [&](int r, int c) -> double& { return t[(r) * ldt + (c)]; };
    for (info = 1; info <= n; ++info)
        if (T(info, info) == 0.0) return;
    info = 0;
    int cs = 1;
    if (job % 10 != 0) cs = 2;
    if ((job % 100) / 10 != 0) cs += 2;
    switch (cs) {
        case 2:  // solve T*x=b, upper triangular
            b[n] = b[n] / T(n, n);
            for (int jj = 2; jj <= n; ++jj) {
                int j = n - jj + 1;
                double temp = -b[j + 1];
                daxpy(j, temp, &T(1, j + 1), ldt, &b[1], 1);
                b[j] = b[j] / T(j, j);
            }
            break;
        case 3:  // solve T^T*x=b, upper triangular
            b[n] = b[n] / T(n, n);
            for (int jj = 2; jj <= n; ++jj) {
                int j = n - jj + 1;
                b[j] = b[j] - ddot(jj - 1, &T(j + 1, j), ldt, &b[j + 1], 1);
                b[j] = b[j] / T(j, j);
            }
            break;
        case 4:  // solve T^T*x=b, lower triangular
            b[1] = b[1] / T(1, 1);
            for (int j = 2; j <= n; ++j) {
                b[j] = b[j] - ddot(j - 1, &T(1, j), ldt, &b[1], 1);
                b[j] = b[j] / T(j, j);
            }
            break;
        default:  // solve T*x=b, lower triangular
            b[1] = b[1] / T(1, 1);
            for (int j = 2; j <= n; ++j) {
                double temp = -b[j - 1];
                daxpy(n - j + 1, temp, &T(j, j - 1), ldt, &b[j], 1);
                b[j] = b[j] / T(j, j);
            }
            break;
    }
}

// projgr: infinity-norm of the projected gradient (box constraints).
void projgr(int n, const double* l, const double* u, const int* nbd,
            const double* x, const double* g, double& sbgnrm) {
    sbgnrm = 0.0;
    for (int i = 1; i <= n; ++i) {
        double gi = g[i];
        if (nbd[i] != 0) {
            if (gi < 0.0) {
                if (nbd[i] >= 2) gi = std::max(x[i] - u[i], gi);
            } else if (nbd[i] <= 2) {
                gi = std::min(x[i] - l[i], gi);
            }
        }
        sbgnrm = std::max(sbgnrm, std::fabs(gi));
    }
}

// errclb: validate L-BFGS-B input arguments.
void errclb(int n, int m, double factr, const double* l, const double* u,
            const int* nbd, std::string& task, int& info, int& k) {
    if (n <= 0) task = "ERROR: N .LE. 0";
    if (m <= 0) task = "ERROR: M .LE. 0";
    if (factr < 0.0) task = "ERROR: FACTR .LT. 0";
    for (int i = 1; i <= n; ++i) {
        if (nbd[i] < 0 || nbd[i] > 3) { task = "ERROR: INVALID NBD"; info = -6; k = i; }
        if (nbd[i] == 2 && l[i] > u[i]) { task = "ERROR: NO FEASIBLE SOLUTION"; info = -7; k = i; }
    }
}

// active: project x to feasible set; classify variables (free/fixed/at bound).
void active(int n, const double* l, const double* u, const int* nbd, double* x,
            int* iwhere, int iprint, bool& prjctd, bool& cnstnd, bool& boxed) {
    int nbdd = 0;
    prjctd = false; cnstnd = false; boxed = true;
    for (int i = 1; i <= n; ++i) {
        if (nbd[i] > 0) {
            if (nbd[i] <= 2 && x[i] <= l[i]) {
                if (x[i] < l[i]) { prjctd = true; x[i] = l[i]; }
                ++nbdd;
            } else if (nbd[i] >= 2 && x[i] >= u[i]) {
                if (x[i] > u[i]) { prjctd = true; x[i] = u[i]; }
                ++nbdd;
            }
        }
    }
    for (int i = 1; i <= n; ++i) {
        if (nbd[i] != 2) boxed = false;
        if (nbd[i] == 0) {
            iwhere[i] = -1;
        } else {
            cnstnd = true;
            if (nbd[i] == 2 && u[i] - l[i] <= 0.0) iwhere[i] = 3;
            else iwhere[i] = 0;
        }
    }
    if (iprint >= 0) {
        if (prjctd) std::printf("The initial X is infeasible. Restart with its projection.\n");
        if (!cnstnd) std::printf("This problem is unconstrained.\n");
    }
    if (iprint > 0) std::printf("At X0 %d variables are exactly at the bounds\n", nbdd);
}

// hpsolb: one step of heap-sort on (t, iorder).
void hpsolb(int n, double* t, int* iorder, int iheap) {
    if (iheap == 0) {
        for (int k = 2; k <= n; ++k) {
            double ddum = t[k];
            int indxin = iorder[k];
            int i = k;
            while (i > 1) {
                int j = i / 2;
                if (ddum < t[j]) { t[i] = t[j]; iorder[i] = iorder[j]; i = j; }
                else break;
            }
            t[i] = ddum;
            iorder[i] = indxin;
        }
    }
    if (n <= 1) return;
    int i = 1;
    double out = t[1];
    int indxou = iorder[1];
    double ddum = t[n];
    int indxin = iorder[n];
    while (true) {
        int j = i + i;
        if (j > n - 1) break;
        if (t[j + 1] < t[j]) ++j;
        if (t[j] < ddum) { t[i] = t[j]; iorder[i] = iorder[j]; i = j; }
        else break;
    }
    t[i] = ddum;
    iorder[i] = indxin;
    t[n] = out;
    iorder[n] = indxou;
}

// freev: update the free/active variable sets; count entering/leaving vars.
void freev(int n, int& nfree, int* index, int& nenter, int& ileave, int* indx2,
           const int* iwhere, bool& wrk, bool updatd, bool cnstnd, int iprint, int iter) {
    nenter = 0;
    ileave = n + 1;
    if (iter > 0 && cnstnd) {
        for (int i = 1; i <= nfree; ++i) {
            int k = index[i];
            if (iwhere[k] > 0) {
                --ileave;
                indx2[ileave] = k;
                if (iprint >= 100) std::printf("Variable %d leaves the set of free variables\n", k);
            }
        }
        for (int i = nfree + 1; i <= n; ++i) {
            int k = index[i];
            if (iwhere[k] <= 0) {
                ++nenter;
                indx2[nenter] = k;
                if (iprint >= 100) std::printf("Variable %d enters the set of free variables\n", k);
            }
        }
        if (iprint >= 99)
            std::printf("%d variables leave; %d variables enter\n", n + 1 - ileave, nenter);
    }
    wrk = (ileave < n + 1) || (nenter > 0) || updatd;
    nfree = 0;
    int iact = n + 1;
    for (int i = 1; i <= n; ++i) {
        if (iwhere[i] <= 0) { ++nfree; index[nfree] = i; }
        else { --iact; index[iact] = i; }
    }
    if (iprint >= 99) std::printf("%d variables are free at GCP %d\n", nfree, iter + 1);
}

// formt: form T = theta*SS + L*D^-1*L' in upper triangle of wt, Cholesky factorize.
void formt(int m, double* wt, const double* sy, const double* ss, int col, double theta, int& info) {
    auto W = [&](int i, int j) -> double& { return wt[(i) * m + (j)]; };
    for (int j = 1; j <= col; ++j) W(1, j) = theta * ss[1 * m + j];
    for (int i = 2; i <= col; ++i) {
        for (int j = i; j <= col; ++j) {
            int k1 = std::min(i, j) - 1;
            double ddum = 0.0;
            for (int k = 1; k <= k1; ++k)
                ddum += sy[i * m + k] * sy[j * m + k] / sy[k * m + k];
            W(i, j) = ddum + theta * ss[i * m + j];
        }
    }
    dpofa(wt, m, col, info);
    if (info != 0) info = -3;
}

// matupd: append a new (s,y) BFGS pair; rotate the circular buffer (head/itail).
void matupd(int n, int m, double* ws, double* wy, double* sy, double* ss,
            double* d, double* r, int& itail, int iupdat, int& col, int& head,
            double& theta, double rr, double dr, double stp, double dtd) {
    auto SS = [&](int i, int j) -> double& { return ss[(i) * m + (j)]; };
    auto SY = [&](int i, int j) -> double& { return sy[(i) * m + (j)]; };
    auto WS = [&](int i, int j) -> double& { return ws[(i) * n + (j)]; };
    auto WY = [&](int i, int j) -> double& { return wy[(i) * n + (j)]; };
    if (iupdat <= m) { col = iupdat; itail = ((head + iupdat - 2) % m) + 1; }
    else { itail = (itail % m) + 1; head = (head % m) + 1; }
    for (int i = 1; i <= n; ++i) WS(i, itail) = d[i];
    for (int i = 1; i <= n; ++i) WY(i, itail) = r[i];
    theta = rr / dr;
    if (iupdat > m) {
        for (int j = 1; j <= col - 1; ++j) {
            for (int k = 1; k <= j; ++k) SS(k, j) = SS(k + 1, j + 1);
            for (int k = 1; k <= col - j; ++k) SY(j, j + k - 1) = SY(j + 1, j + k);
        }
    }
    int pointr = head;
    for (int j = 1; j <= col - 1; ++j) {
        double sd = 0.0, wd = 0.0;
        for (int i = 1; i <= n; ++i) { sd += d[i] * WY(i, pointr); wd += WS(i, pointr) * d[i]; }
        SY(col, j) = sd;
        SS(j, col) = wd;
        pointr = (pointr % m) + 1;
    }
    SS(col, col) = (std::fabs(stp - 1.0) < 1e-20) ? dtd : stp * stp * dtd;
    SY(col, col) = dr;
}

// bmv: product with middle matrix (L-BFGS-B).  Stub: real translation of
// lbfgs.F90 bmv pending in M07; no-op keeps mainlb/cmprlb linkable.
void bmv(int m, const double* sy, const double* wt, int col, double* v, double* wa, int& info) {
    info = 0;
    (void)m; (void)sy; (void)wt; (void)col; (void)v; (void)wa;
}

// cmprlb: compute the compressed search direction r.
void cmprlb(int n, int m, const double* x, const double* g, double* ws, double* wy,
            const double* sy, double* wt, const double* z, double* r, double* wa,
            const int* index, double theta, int col, int head, int nfree, bool cnstnd, int& info) {
    if (!cnstnd && col > 0) {
        for (int i = 1; i <= n; ++i) r[i] = -g[i];
    } else {
        for (int i = 1; i <= nfree; ++i) {
            int k = index[i];
            r[i] = -theta * (z[k] - x[k]) - g[k];
        }
        // copy r[1..nfree] -> wa[2*m+1 .. 2*m+nfree]
        for (int i = 1; i <= nfree; ++i) wa[2 * m + i] = r[i];
        bmv(m, sy, wt, col, wa + 2 * m, wa, info);
        if (info != 0) { info = -8; return; }
        int pointr = head;
        for (int j = 1; j <= col; ++j) {
            double a1 = wa[j];
            double a2 = theta * wa[col + j];
            for (int i = 1; i <= nfree; ++i) {
                int k = index[i];
                r[i] += wy[k * n + pointr] * a1 + ws[k * n + pointr] * a2;
            }
            pointr = (pointr % m) + 1;
        }
    }
}

// formk: build and factor the preconditioned K (Schur) matrix WN.
void formk(int n, int nsub, const int* ind, int nenter, int ileave, const int* indx2,
           int iupdat, bool updatd, double* wn, double* wn1, int m,
           double* ws, double* wy, const double* sy, double theta, int& col, int head, int& info) {
    int m2 = 2 * m;
    auto WN = [&](int i, int j) -> double& { return wn[(i) * m2 + (j)]; };
    auto WN1 = [&](int i, int j) -> double& { return wn1[(i) * m2 + (j)]; };
    int upcl, ipntr = head, jpntr = head;
    if (updatd) {
        if (iupdat > m) {
            for (int jy = 1; jy <= m - 1; ++jy) {
                int js = m + jy;
                for (int k = 0; k < m - jy; ++k) WN1(jy + 1 + k, jy + 1) = WN1(jy + k, jy);
                for (int k = 0; k < m - jy; ++k) WN1(js + 1 + k, js + 1) = WN1(js + k, js);
                for (int k = 0; k < m - 1; ++k) WN1(m + 2 + k, jy + 1) = WN1(m + 1 + k, jy);
            }
        }
        int pbegin = 1, pend = nsub, dbegin = nsub + 1, dend = n;
        int iy = col, is = m + col;
        int ipntr = head + col - 1; if (ipntr > m) ipntr -= m;
        int jpntr = head;
        for (int jy = 1; jy <= col; ++jy) {
            int js = m + jy;
            double t1 = 0, t2 = 0, t3 = 0;
            for (int k = pbegin; k <= pend; ++k) { int k1 = ind[k]; t1 += wy[k1 * n + ipntr] * wy[k1 * n + jpntr]; }
            for (int k = dbegin; k <= dend; ++k) { int k1 = ind[k]; t2 += ws[k1 * n + ipntr] * ws[k1 * n + jpntr]; t3 += ws[k1 * n + ipntr] * wy[k1 * n + jpntr]; }
            WN1(iy, jy) = t1; WN1(is, js) = t2; WN1(is, jy) = t3;
            jpntr = (jpntr % m) + 1;
        }
        {
            int jy = col;
            jpntr = head + col - 1; if (jpntr > m) jpntr -= m;
            ipntr = head;
            for (int i = 1; i <= col; ++i) {
                is = m + i; double t3 = 0;
                for (int k = pbegin; k <= pend; ++k) { int k1 = ind[k]; t3 += ws[k1 * n + ipntr] * wy[k1 * n + jpntr]; }
                ipntr = (ipntr % m) + 1;
                WN1(is, jy) = t3;
            }
        }
        upcl = col - 1;
    } else {
        upcl = col;
    }
    ipntr = head;
    for (int iy = 1; iy <= upcl; ++iy) {
        int is = m + iy; jpntr = head;
        for (int jy = 1; jy <= iy; ++jy) {
            int js = m + jy; double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
            for (int k = 1; k <= nenter; ++k) { int k1 = indx2[k]; t1 += wy[k1 * n + ipntr] * wy[k1 * n + jpntr]; t2 += ws[k1 * n + ipntr] * ws[k1 * n + jpntr]; }
            for (int k = ileave; k <= n; ++k) { int k1 = indx2[k]; t3 += wy[k1 * n + ipntr] * wy[k1 * n + jpntr]; t4 += ws[k1 * n + ipntr] * ws[k1 * n + jpntr]; }
            WN1(iy, jy) += t1 - t3;
            WN1(is, js) -= t2 - t4;
            jpntr = (jpntr % m) + 1;
        }
        ipntr = (ipntr % m) + 1;
    }
    ipntr = head;
    for (int is = m + 1; is <= m + upcl; ++is) {
        jpntr = head;
        for (int jy = 1; jy <= upcl; ++jy) {
            double t1 = 0, t3 = 0;
            for (int k = 1; k <= nenter; ++k) { int k1 = indx2[k]; t1 += ws[k1 * n + ipntr] * wy[k1 * n + jpntr]; }
            for (int k = ileave; k <= n; ++k) { int k1 = indx2[k]; t3 += ws[k1 * n + ipntr] * wy[k1 * n + jpntr]; }
            if (is <= jy + m) WN1(is, jy) += t1 - t3; else WN1(is, jy) -= t1 - t3;
            jpntr = (jpntr % m) + 1;
        }
        ipntr = (ipntr % m) + 1;
    }
    for (int iy = 1; iy <= col; ++iy) {
        int is = col + iy, is1 = m + iy;
        for (int jy = 1; jy <= iy; ++jy) {
            int js = col + jy, js1 = m + jy;
            WN(jy, iy) = WN1(iy, jy) / theta;
            WN(js, is) = WN1(is1, js1) * theta;
        }
        for (int jy = 1; jy <= iy - 1; ++jy) WN(jy, is) = -WN1(is1, jy);
        for (int jy = iy; jy <= col; ++jy) WN(jy, is) = WN1(is1, jy);
        WN(iy, iy) += sy[iy * m + iy];
    }
    dpofa(wn, m2, col, info);
    if (info != 0) { info = -1; return; }
    int col2 = 2 * col;
    for (int js = col + 1; js <= col2; ++js) dtrsl(wn, m2, col, &WN(1, js), 11, info);
    for (int is = col + 1; is <= col2; ++is)
        for (int js = is; js <= col2; ++js) {
            double acc = 0; for (int k = 1; k <= col; ++k) acc += WN(k, is) * WN(k, js);
            WN(is, js) += acc;
        }
    dpofa(&WN(col + 1, col + 1), m2, col, info);
    if (info != 0) info = -2;
}

// cauchy: find the Generalized Cauchy Point (GCP) xcp.
void cauchy(int n, double* x, const double* l, const double* u, const int* nbd,
            const double* g, int* iorder, int* iwhere, double* t, double* d, double* xcp,
            int m, double* wy, double* ws, const double* sy, double* wt, double theta,
            int col, int head, double* p, double* c, double* wbp, double* v, int& nint,
            int iprint, double sbgnrm, int& info, double epsmch) {
    if (sbgnrm <= 0.0) {
        if (iprint >= 0) std::printf("Subgnorm = 0.  GCP = X.\n");
        for (int i = 1; i <= n; ++i) xcp[i] = x[i];
        return;
    }
    bool bnded = true;
    int nfree = n + 1, nbreak = 0, ibkmin = 0;
    double bkmin = 0.0;
    int col2 = 2 * col;
    double f1 = 0.0;
    for (int i = 1; i <= col2; ++i) p[i] = 0.0;
    for (int i = 1; i <= n; ++i) {
        double neggi = -g[i];
        double tl = 0, tu = 0;
        if (iwhere[i] != 3 && iwhere[i] != -1) {
            if (nbd[i] <= 2) tl = x[i] - l[i];
            if (nbd[i] >= 2) tu = u[i] - x[i];
            bool xlower = (nbd[i] <= 2 && tl <= 0.0);
            bool xupper = (nbd[i] >= 2 && tu <= 0.0);
            iwhere[i] = 0;
            if (xlower) { if (neggi <= 0.0) iwhere[i] = 1; }
            else if (xupper) { if (neggi >= 0.0) iwhere[i] = 2; }
            else if (std::fabs(neggi) <= 0.0) iwhere[i] = -3;
        }
        int pointr = head;
        if (iwhere[i] != 0 && iwhere[i] != -1) {
            d[i] = 0.0;
        } else {
            d[i] = neggi;
            f1 -= neggi * neggi;
            for (int j = 1; j <= col; ++j) {
                p[j] += wy[i * n + pointr] * neggi;
                p[col + j] += ws[i * n + pointr] * neggi;
                pointr = (pointr % m) + 1;
            }
            if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.0) {
                ++nbreak; iorder[nbreak] = i; t[nbreak] = tl / (-neggi);
                if (nbreak == 1 || t[nbreak] < bkmin) { bkmin = t[nbreak]; ibkmin = nbreak; }
            } else if (nbd[i] >= 2 && neggi > 0.0) {
                ++nbreak; iorder[nbreak] = i; t[nbreak] = tu / neggi;
                if (nbreak == 1 || t[nbreak] < bkmin) { bkmin = t[nbreak]; ibkmin = nbreak; }
            } else {
                --nfree; iorder[nfree] = i;
                if (std::fabs(neggi) > 0.0) bnded = false;
            }
        }
    }
    if (std::fabs(theta - 1.0) > 1e-20)
        for (int i = 1; i <= col; ++i) p[col + i] *= theta;
    for (int i = 1; i <= n; ++i) xcp[i] = x[i];
    if (nbreak == 0 && nfree == n + 1) return;
    for (int j = 1; j <= col2; ++j) c[j] = 0.0;
    double f2 = -theta * f1, f2_org = f2;
    if (col > 0) {
        bmv(m, sy, wt, col, p, v, info);
        if (info != 0) return;
        double acc = 0; for (int i = 1; i <= col2; ++i) acc += v[i] * p[i];
        f2 -= acc;
    }
    double dtm = -f1 / f2, tsum = 0.0;
    nint = 1;
    if (nbreak == 0) return;
    int nleft = nbreak, iter = 1;
    double tj = 0.0, dt = 0.0;
    bool gcp_in_segment = false;
    while (true) {
        double tj0 = tj;
        int ibp;
        if (iter == 1) {
            tj = bkmin; ibp = iorder[ibkmin];
        } else {
            if (iter == 2 && ibkmin != nbreak) {
                t[ibkmin] = t[nbreak]; iorder[ibkmin] = iorder[nbreak];
            }
            hpsolb(nleft, t, iorder, iter - 2);
            tj = t[nleft]; ibp = iorder[nleft];
        }
        dt = tj - tj0;
        if (dtm < dt) { gcp_in_segment = true; break; }
        tsum += dt; --nleft; ++iter;
        double dibp = d[ibp]; d[ibp] = 0.0;
        double zibp;
        if (dibp > 0.0) { zibp = u[ibp] - x[ibp]; xcp[ibp] = u[ibp]; iwhere[ibp] = 2; }
        else { zibp = l[ibp] - x[ibp]; xcp[ibp] = l[ibp]; iwhere[ibp] = 1; }
        if (nleft == 0 && nbreak == n) break;
        ++nint;
        double dibp2 = dibp * dibp;
        f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
        f2 = f2 - theta * dibp2;
        if (col > 0) {
            for (int i = 1; i <= col2; ++i) c[i] += dt * p[i];
            int pointr = head;
            for (int j = 1; j <= col; ++j) {
                wbp[j] = wy[ibp * n + pointr];
                wbp[col + j] = theta * ws[ibp * n + pointr];
                pointr = (pointr % m) + 1;
            }
            bmv(m, sy, wt, col, wbp, v, info);
            if (info != 0) return;
            double wmc = 0, wmp = 0, wmw = 0;
            for (int i = 1; i <= col2; ++i) { wmc += c[i] * v[i]; wmp += p[i] * v[i]; wmw += wbp[i] * v[i]; }
            for (int i = 1; i <= col2; ++i) p[i] -= dibp * wbp[i];
            f1 += dibp * wmc;
            f2 += 2.0 * dibp * wmp - dibp2 * wmw;
        }
        f2 = std::max(epsmch * f2_org, f2);
        if (nleft > 0) { dtm = -f1 / f2; continue; }
        if (bnded) { f1 = 0; f2 = 0; dtm = 0; }
        else dtm = -f1 / f2;
    }
    if (!gcp_in_segment) dtm = dt;
    if (dtm < 0.0) dtm = 0.0;
    tsum += dtm;
    for (int i = 1; i <= n; ++i) xcp[i] += tsum * d[i];
    if (col > 0) for (int i = 1; i <= col2; ++i) c[i] += dtm * p[i];
}

// lnsrlb: bounded line search driver around the GCP direction d.
void lnsrlb(int n, const double* l, const double* u, const int* nbd, double* x,
            double f, double& fold, double& gd, double& gdold, const double* g,
            const double* d, double* r, double* t, const double* z, double& stp,
            double& dnorm, double& dtd, double& xstep, double& stpmx, int iter,
            int& ifun, int& iback, int& nfgv, int& info, std::string& task,
            bool boxed, bool cnstnd, std::string& csave, int* isave, double* dsave) {
    if (task.rfind("FG_LN", 0) != 0) {
        double acc = 0; for (int i = 1; i <= n; ++i) acc += d[i] * d[i];
        dtd = acc; dnorm = std::sqrt(dtd);
        stpmx = 1e-5;
        if (cnstnd) {
            if (iter == 0) stpmx = 1.0;
            else {
                for (int i = 1; i <= n; ++i) {
                    double a1 = d[i];
                    if (nbd[i] != 0) {
                        if (a1 < 0.0 && nbd[i] <= 2) {
                            double a2 = l[i] - x[i];
                            if (a2 >= 0.0) stpmx = 0.0;
                            else if (a1 * stpmx < a2) stpmx = a2 / a1;
                        } else if (a1 > 0.0 && nbd[i] >= 2) {
                            double a2 = u[i] - x[i];
                            if (a2 <= 0.0) stpmx = 0.0;
                            else if (a1 * stpmx > a2) stpmx = a2 / a1;
                        }
                    }
                }
            }
        }
        if (iter == 0 && !boxed) stp = std::min(1.0 / dnorm, stpmx);
        else stp = 1.0;
        for (int i = 1; i <= n; ++i) { t[i] = x[i]; r[i] = g[i]; }
        fold = f; ifun = 0; iback = 0; csave = "START";
    }
    { double acc = 0; for (int i = 1; i <= n; ++i) acc += g[i] * d[i]; gd = acc; }
    if (ifun == 0) {
        gdold = gd;
        if (gd >= 0.0) { info = -4; return; }
    }
    dcsrch(f, gd, stp, 1e-3, 0.9, 0.1, 0.0, stpmx, csave, isave, dsave);
    xstep = stp * dnorm;
    if (csave.rfind("CONV", 0) == 0 || csave.rfind("WARN", 0) == 0) { task = "NEW_X"; return; }
    task = "FG_LNSRCH";
    ++ifun; ++nfgv; iback = ifun - 1;
    if (std::fabs(stp - 1.0) < 1e-20) { for (int i = 1; i <= n; ++i) x[i] = z[i]; }
    else { for (int i = 1; i <= n; ++i) x[i] = stp * d[i] + t[i]; }
}

// bmv: apply the middle matrix inverse (two triangular solves with sy,wt).
void bmv(int m, const double* sy, double* wt, int col, const double* v, double* p, int& info) {
    if (col == 0) return;
    auto SY = [&](int i, int k) -> double { return sy[(i) * m + (k)]; };
    p[col + 1] = v[col + 1];
    for (int i = 2; i <= col; ++i) {
        int i2 = col + i; double sum = 0.0;
        for (int k = 1; k <= i - 1; ++k) sum += SY(i, k) * v[k] / SY(k, k);
        p[i2] = v[i2] + sum;
    }
    dtrsl(wt, m, col, &p[col + 1], 11, info);
    if (info != 0) return;
    for (int i = 1; i <= col; ++i) p[i] = v[i] / std::sqrt(SY(i, i));
    dtrsl(wt, m, col, &p[col + 1], 1, info);
    if (info != 0) return;
    for (int i = 1; i <= col; ++i) p[i] = -p[i] / std::sqrt(SY(i, i));
    for (int i = 1; i <= col; ++i) {
        double sum = 0.0;
        for (int k = i + 1; k <= col; ++k) sum += SY(k, i) * p[col + k] / SY(i, i);
        p[i] += sum;
    }
}

// setulb: lay out the work-space wa, then dispatch to mainlb.
void mainlb(int n, int m, double* x, const double* l, const double* u, const int* nbd,
            double& f, double* g, double factr, double pgtol, double* ws, double* wy,
            double* sy, double* ss, double* wt, double* wn, double* snd, double* z,
            double* r, double* d, double* t, double* wa, int* iwa1, int* iwa2, int* iwa3,
            std::string& task, int iprint, std::string& csave, bool* lsave, int* isave, double* dsave);

void setulb(int n, int m, double* x, const double* l, const double* u, const int* nbd,
            double& f, double* g, double factr, double pgtol, double* wa, int* iwa,
            std::string& task, int iprint, std::string& csave, bool* lsave, int* isave, double* dsave) {
    if (task == "START") {
        isave[1] = m * n; isave[2] = m * m; isave[3] = 4 * m * m; isave[4] = 1;
        isave[5] = isave[4] + isave[1]; isave[6] = isave[5] + isave[1];
        isave[7] = isave[6] + isave[2]; isave[8] = isave[7] + isave[2];
        isave[9] = isave[8]; isave[10] = isave[9] + isave[2];
        isave[11] = isave[10] + isave[3]; isave[12] = isave[11] + isave[3];
        isave[13] = isave[12] + n; isave[14] = isave[13] + n;
        isave[15] = isave[14] + n; isave[16] = isave[15] + n;
    }
    int lws = isave[4], lwy = isave[5], lsy = isave[6], lss = isave[7];
    int lwt = isave[9], lwn = isave[10], lsnd = isave[11], lz = isave[12];
    int lr = isave[13], ld = isave[14], lt = isave[15], lwa = isave[16];
    mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa + lws, wa + lwy, wa + lsy,
           wa + lss, wa + lwt, wa + lwn, wa + lsnd, wa + lz, wa + lr, wa + ld,
           wa + lt, wa + lwa, iwa + 1, iwa + (n + 1), iwa + (2 * n + 1),
           task, iprint, csave, lsave, isave + 22, dsave);
}

// External helpers used by mainlb.
void subsm(int n, int m, int nfree, const int* index, const double* l, const double* u,
           const int* nbd, double* z, double* r, double* ws, double* wy, double theta,
           int col, int head, int& iword, double* wa, double* wn, int iprint, int& info);
void prn1lb(int n, int m, const double* l, const double* u, double* x, int iprint,
            int itfile, double epsmch) {
    // L-BFGS-B iteration header printout.  Stub: real translation pending (M07).
    (void)n; (void)m; (void)l; (void)u; (void)x; (void)iprint; (void)itfile; (void)epsmch;
}
void prn2lb(int n, double* x, double f, const double* g, int iprint, int itfile, int iter,
            int nfgv, int nact, double sbgnrm, int nint, const char* word, int iword,
            int iback, double stp, double xstep) {
    // L-BFGS-B iteration line printout.  Stub: real translation pending (M07).
    (void)n; (void)x; (void)f; (void)g; (void)iprint; (void)itfile; (void)iter; (void)nfgv;
    (void)nact; (void)sbgnrm; (void)nint; (void)word; (void)iword; (void)iback; (void)stp; (void)xstep;
}
void prn3lb(int n, double* x, double f, std::string& task, int iprint, int info, int itfile,
            int iter, int nfgv, int nintol, int nskip, int nact, double sbgnrm, double time,
            int nint, const char* word, int iback, double stp, double xstep, int k,
            double cachyt, double sbtime, double lnscht) {
    // L-BFGS-B final printout.  Stub: real translation pending (M07).
    (void)n; (void)x; (void)f; (void)task; (void)iprint; (void)info; (void)itfile; (void)iter;
    (void)nfgv; (void)nintol; (void)nskip; (void)nact; (void)sbgnrm; (void)time; (void)nint;
    (void)word; (void)iback; (void)stp; (void)xstep; (void)k; (void)cachyt; (void)sbtime; (void)lnscht;
}

// mainlb: L-BFGS-B state machine (START / FG_ST / NEW_X / STOP / FG_LN).
void mainlb(int n, int m, double* x, const double* l, const double* u, const int* nbd,
            double& f, double* g, double factr, double pgtol, double* ws, double* wy,
            double* sy, double* ss, double* wt, double* wn, double* snd, double* z,
            double* r, double* d, double* t, double* wa, int* index, int* iwhere, int* indx2,
            std::string& task, int iprint, std::string& csave, bool* lsave, int* isave, double* dsave) {
    const double one = 1.0, zero = 0.0;
    std::string word = "---";
    bool boxed, cnstnd, prjctd, updatd, wrk;
    int col, head, i, iback, ifun, ileave, info, itail, iter, itfile, iupdat, iword, k = 0;
    int nact, nenter, nfgv, nfree, nint, nintol, nskip;
    double ddum, dnorm, dr, dtd, epsmch, fold, gd, gdold, rr, sbgnrm, sbtime, stp, stpmx;
    double theta, tol, xstep, cachyt = 0.0, lnscht = 0.0;

    if (task == "START") {
        epsmch = dpmeps();
        fold = 0; dnorm = 0; gd = 0; sbgnrm = 0; stp = 0; stpmx = 0; gdold = 0; dtd = 0;
        col = 0; head = 1; theta = one; iupdat = 0; updatd = false; iback = 0; itail = 0;
        ifun = 0; iword = 0; nact = 0; ileave = 0; nenter = 0;
        iter = 0; nfgv = 0; nint = 0; nintol = 0; nskip = 0; nfree = n;
        tol = factr * epsmch;
        cachyt = 0; sbtime = 0; lnscht = 0;
        word = "---"; info = 0; itfile = 0;
        errclb(n, m, factr, l, u, nbd, task, info, k);
        if (task.rfind("ERROR", 0) == 0) {
            xstep = 0; k = 0;
            prn3lb(n, x, f, task, iprint, info, itfile, iter, nfgv, nintol, nskip, nact,
                   sbgnrm, zero, nint, word.c_str(), iback, stp, xstep, k, cachyt, sbtime, lnscht);
            goto save;
        } else {
            prn1lb(n, m, l, u, x, iprint, itfile, epsmch);
            active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed);
        }
    } else {
        prjctd = lsave[1]; cnstnd = lsave[2]; boxed = lsave[3]; updatd = lsave[4];
        nintol = isave[1]; itfile = isave[3]; iback = isave[4]; nskip = isave[5];
        head = isave[6]; col = isave[7]; itail = isave[8]; iter = isave[9]; iupdat = isave[10];
        nint = isave[12]; nfgv = isave[13]; info = isave[14]; ifun = isave[15]; iword = isave[16];
        nfree = isave[17]; nact = isave[18]; ileave = isave[19]; nenter = isave[20];
        theta = dsave[1]; fold = dsave[2]; tol = dsave[3]; dnorm = dsave[4]; epsmch = dsave[5];
        cachyt = dsave[7]; sbtime = dsave[8]; lnscht = dsave[9];
        gd = dsave[11]; stpmx = dsave[12]; sbgnrm = dsave[13]; stp = dsave[14];
        gdold = dsave[15]; dtd = dsave[16];
        if (task.rfind("FG_LN", 0) == 0) goto save;

        if (task.rfind("NEW_X", 0) == 0) {
            if (sbgnrm <= pgtol) {
                task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"; goto done;
            }
            ddum = std::max({std::fabs(fold), std::fabs(f), one});
            if ((fold - f) <= tol * ddum) {
                task = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH";
                if (iback >= 10) info = -5;
                goto done;
            }
            for (i = 1; i <= n; ++i) r[i] = g[i] - r[i];
            { double acc = 0; for (i = 1; i <= n; ++i) acc += r[i] * r[i]; rr = acc; }
            if (std::fabs(stp - one) < 1e-20) { dr = gd - gdold; ddum = -gdold; }
            else { dr = (gd - gdold) * stp; for (i = 1; i <= n; ++i) d[i] *= stp; ddum = -gdold * stp; }
            if (dr <= epsmch * ddum) {
                ++nskip; updatd = false;
            } else {
                updatd = true; ++iupdat;
                matupd(n, m, ws, wy, sy, ss, d, r, itail, iupdat, col, head, theta, rr, dr, stp, dtd);
                formt(m, wt, sy, ss, col, theta, info);
                if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; }
            }
        } else if (task.rfind("FG_ST", 0) == 0) {
            nfgv = 1;
            projgr(n, l, u, nbd, x, g, sbgnrm);
            if (sbgnrm <= pgtol) { task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"; goto done; }
        } else if (task.rfind("STOP", 0) == 0) {
            for (i = 1; i <= n; ++i) { x[i] = t[i]; g[i] = r[i]; }
            f = fold;
            goto done;
        } else {
            goto fg_start;
        }

        // iteration / GCP loop
        while (true) {
            iword = -1;
            if (!cnstnd && col > 0) {
                for (i = 1; i <= n; ++i) z[i] = x[i];
                wrk = updatd; nint = 0;
            } else {
                cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, ws, sy, wt, theta,
                       col, head, wa + 1, wa + (2 * m + 1), wa + (4 * m + 1), wa + (6 * m + 1),
                       nint, iprint, sbgnrm, info, epsmch);
                if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; continue; }
                nintol += nint;
                freev(n, nfree, index, nenter, ileave, indx2, iwhere, wrk, updatd, cnstnd, iprint, iter);
                nact = n - nfree;
            }
            if (nfree == 0 || col == 0) break;
            bool sb_fail = false;
            if (wrk) formk(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd, m,
                           ws, wy, sy, theta, col, head, info);
            if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; sb_fail = true; }
            else {
                cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, col, head, nfree, cnstnd, info);
                if (info == 0) subsm(n, m, nfree, index, l, u, nbd, z, r, ws, wy, theta, col, head, iword, wa, wn, iprint, info);
                if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; sb_fail = true; }
                else break;
            }
            (void)sb_fail;
        }
        for (i = 1; i <= n; ++i) d[i] = z[i] - x[i];

        // line search loop
        while (true) {
            lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, stp, dnorm, dtd,
                   xstep, stpmx, iter, ifun, iback, nfgv, info, task, boxed, cnstnd, csave,
                   isave + 22, dsave + 17);
            if (info == 0 && iback < 20) goto after_ls;
            for (i = 1; i <= n; ++i) { x[i] = t[i]; g[i] = r[i]; }
            f = fold;
            if (col == 0) break;
            if (info == 0) --nfgv;
            info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false;
            // recompute GCP + subspace after bad direction
            while (true) {
                iword = -1;
                if (!cnstnd && col > 0) { for (i = 1; i <= n; ++i) z[i] = x[i]; wrk = updatd; nint = 0; }
                else {
                    cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, ws, sy, wt, theta,
                           col, head, wa + 1, wa + (2 * m + 1), wa + (4 * m + 1), wa + (6 * m + 1),
                           nint, iprint, sbgnrm, info, epsmch);
                    if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; continue; }
                    nintol += nint;
                    freev(n, nfree, index, nenter, ileave, indx2, iwhere, wrk, updatd, cnstnd, iprint, iter);
                    nact = n - nfree;
                }
                if (nfree == 0 || col == 0) break;
                if (wrk) formk(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd, m,
                               ws, wy, sy, theta, col, head, info);
                if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; }
                else {
                    cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, col, head, nfree, cnstnd, info);
                    if (info == 0) subsm(n, m, nfree, index, l, u, nbd, z, r, ws, wy, theta, col, head, iword, wa, wn, iprint, info);
                    if (info != 0) { info = 0; col = 0; head = 1; theta = one; iupdat = 0; updatd = false; }
                    else break;
                }
            }
            for (i = 1; i <= n; ++i) d[i] = z[i] - x[i];
        }
        if (info == 0) { info = -9; --nfgv; --ifun; --iback; }
        task = "ABNORMAL_TERMINATION_IN_LNSRCH";
        ++iter;
        goto done;

after_ls:
        if (task.rfind("FG_LN", 0) != 0) {
            ++iter;
            projgr(n, l, u, nbd, x, g, sbgnrm);
            prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, sbgnrm, nint,
                   word.c_str(), iword, iback, stp, xstep);
        }
        goto save;

done:
        prn3lb(n, x, f, task, iprint, info, itfile, iter, nfgv, nintol, nskip, nact,
               sbgnrm, zero, nint, word.c_str(), iback, stp, xstep, k, cachyt, sbtime, lnscht);
        goto save;

fg_start:
        task = "FG_START";
    }

save:
        lsave[1] = prjctd; lsave[2] = cnstnd; lsave[3] = boxed; lsave[4] = updatd;
        isave[1] = nintol; isave[3] = itfile; isave[4] = iback; isave[5] = nskip;
        isave[6] = head; isave[7] = col; isave[8] = itail; isave[9] = iter; isave[10] = iupdat;
        isave[12] = nint; isave[13] = nfgv; isave[14] = info; isave[15] = ifun; isave[16] = iword;
        isave[17] = nfree; isave[18] = nact; isave[19] = ileave; isave[20] = nenter;
        dsave[1] = theta; dsave[2] = fold; dsave[3] = tol; dsave[4] = dnorm; dsave[5] = epsmch;
        dsave[7] = cachyt; dsave[8] = sbtime; dsave[9] = lnscht;
        dsave[11] = gd; dsave[12] = stpmx; dsave[13] = sbgnrm; dsave[14] = stp;
        dsave[15] = gdold; dsave[16] = dtd;
}

// subsm: solve the preconditioned subspace system; backtrack to the box.
void subsm(int n, int m, int nsub, const int* ind, const double* l, const double* u,
           const int* nbd, double* x, double* d, double* ws, double* wy, double theta,
           int col, int head, int& iword, double* wv, double* wn, int iprint, int& info) {
    if (nsub <= 0) return;
    int pointr = head;
    for (int i = 1; i <= col; ++i) {
        double t1 = 0, t2 = 0;
        for (int j = 1; j <= nsub; ++j) { int k = ind[j]; t1 += wy[k * n + pointr] * d[j]; t2 += ws[k * n + pointr] * d[j]; }
        wv[i] = t1; wv[col + i] = theta * t2;
        pointr = (pointr % m) + 1;
    }
    int m2 = 2 * m, col2 = 2 * col;
    dtrsl(wn, m2, col2, wv, 11, info);
    if (info != 0) return;
    for (int i = 1; i <= col; ++i) wv[i] = -wv[i];
    dtrsl(wn, m2, col2, wv, 1, info);
    if (info != 0) return;
    pointr = head;
    for (int jy = 1; jy <= col; ++jy) {
        int js = col + jy;
        for (int i = 1; i <= nsub; ++i) {
            int k = ind[i];
            d[i] += wy[k * n + pointr] * wv[jy] / theta + ws[k * n + pointr] * wv[js];
        }
        pointr = (pointr % m) + 1;
    }
    for (int i = 1; i <= nsub; ++i) d[i] /= theta;
    double alpha = 1.0;
    double temp1 = alpha;
    int ibd = 0;
    for (int i = 1; i <= nsub; ++i) {
        int k = ind[i]; double dk = d[i];
        if (nbd[k] != 0) {
            if (dk < 0.0 && nbd[k] <= 2) {
                double temp2 = l[k] - x[k];
                if (temp2 >= 0.0) temp1 = 0.0;
                else if (dk * alpha < temp2) temp1 = temp2 / dk;
            } else if (dk > 0.0 && nbd[k] >= 2) {
                double temp2 = u[k] - x[k];
                if (temp2 <= 0.0) temp1 = 0.0;
                else if (dk * alpha > temp2) temp1 = temp2 / dk;
            }
            if (temp1 < alpha) { alpha = temp1; ibd = i; }
        }
    }
    if (alpha < 1.0) {
        double dk = d[ibd]; int k = ind[ibd];
        if (dk > 0.0) { x[k] = u[k]; d[ibd] = 0.0; }
        else if (dk < 0.0) { x[k] = l[k]; d[ibd] = 0.0; }
    }
    for (int i = 1; i <= nsub; ++i) { int k = ind[i]; x[k] = x[k] + alpha * d[i]; }
    iword = (alpha < 1.0) ? 1 : 0;
}