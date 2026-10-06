// geout.cpp — C++ translation of MOPAC 2016 "geout.F90".
// Prints the current geometry (MOPAC .out style or restart-readable style),
// including symmetry data and user-supplied pi bonds.
#include "geout.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "elemts_C.h"
#include "chanel_C.h"
#include "maps_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace parameters_C {
extern std::vector<double> ams;   // tore declared in parameters_C.h as double[N_PARAM+1]
}
namespace molkst_C {
extern int numat, natoms, ndep, maxtxt, ncomments;
extern bool gui, in_house_only;
extern std::string keywrd, line;
}
namespace symmetry_C {
extern std::vector<double> depmul;
}
namespace elemts_C {
extern std::vector<std::string> elemnt;
}
namespace chanel_C {
extern int iw;
}
namespace maps_C {
extern int lpara1, latom1, lpara2, latom2, lparam, latom;
}

extern void wrttxt(int);
extern void chrge(const std::vector<double>&, std::vector<double>&);
extern const char* elemnt(int);

namespace {
std::string fixed(double v, const char* fmt) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), fmt, v);
    return buf;
}
}  // namespace

void geout(int mode1) {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace symmetry_C;
    int mode = mode1;
    int store_maxtxt = maxtxt;
    bool charge = (keywrd.find(" PRTCHAR") != std::string::npos);
    bool lxyz = (keywrd.find(" COORD") != std::string::npos) ||
                  (keywrd.find("VELO") != std::string::npos);
    if (keywrd.find(" 0SCF") != std::string::npos && lxyz) {
        natoms = numat;
        lxyz = false;
    }
    std::string flag1, flag0, flagn;
    int iprt;
    if (mode == 1) {
        flag1 = " *"; flag0 = "  "; flagn = " +"; iprt = chanel_C::iw;
    } else {
        if (gui) { flag1 = " 1"; flag0 = " 0"; }
        else { flag1 = "+1"; flag0 = "+0"; }
        flagn = "-1";
        iprt = std::abs(mode);
    }
    double degree = 57.29577951308232;
    if (lxyz) degree = 1.0;
    std::string blank(80, ' ');
    if (mode != 1 && keywrd.find(" RESEQ") != std::string::npos)
        for (int i = 1; i <= numat; ++i) atmass[i] = parameters_C::ams[nat[i]];
    bool isotopes = false;
    if (mode != 1 && nat[1] != 0) {
        int i = 1;
        for (; i <= numat; ++i)
            if (std::fabs(parameters_C::ams[nat[i]] - atmass[i]) > 1.e-6) break;
        isotopes = (i <= numat);
    }
    bool cart = true;
    for (int i = 1; i <= natoms; ++i)
        if (na[i] > 0) cart = false;
    if (mode > 1 && ncomments > 0) {
        for (int i = 1; i <= ncomments; ++i) {
            std::string c = all_comments[i];
            while (!c.empty() && (c.back() == ' ' || c.back() == '\t')) c.pop_back();
            std::fprintf(stdout, "%s\n", c.c_str());
        }
    }
    std::string fmt1 = "%13.8f", fmt23 = "%13.7f";
    if (maxtxt != 0) maxtxt = maxtxt + 2;
    if (cart) {
        double x = 0.0;
        for (int i = 1; i <= natoms; ++i)
            x = std::min(x, std::min(geo[1][i], std::min(geo[2][i], geo[3][i])));
        if (x < -99.99) { fmt1 = "%14.8f"; fmt23 = "%14.8f"; }
        else fmt23 = "%13.8f";
        if (mode == 1) {
            int i = (maxtxt == 0) ? 6 : maxtxt / 2 + 1;
            std::string s1 = "   ATOM " + blank.substr(0, maxtxt / 2 + 1) + " CHEMICAL  " + blank.substr(0, i) + "  X               Y               Z";
            std::string s2 = "  NUMBER " + blank.substr(0, maxtxt / 2 + 1) + " SYMBOL" + blank.substr(0, i) + "(ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)";
            std::fprintf(stdout, "%s\n%s\n\n", s1.c_str(), s2.c_str());
        } else if (mode > 0) wrttxt(iprt);
    } else if (mode == 1) {
        int j = std::max(9, maxtxt + 2);
        if (maxtxt == 0) j = 8;
        std::fprintf(stdout, "%s\n%s\n%s\n",
            ("  ATOM" + blank.substr(0, j / 2) + "CHEMICAL " + blank.substr(0, (j + 1) / 2) + " BOND LENGTH      BOND ANGLE     TWIST ANGLE ").c_str(),
            (" NUMBER" + blank.substr(0, j / 2) + "SYMBOL   " + blank.substr(0, (j + 1) / 2) + "(ANGSTROMS)      (DEGREES)       (DEGREES) ").c_str(),
            ("   (I)       " + blank.substr(0, j) + "      NA:I           NB:NA:I    " + "   NC:NB:NA:I " + "      NA    NB    NC ").c_str());
    } else if (mode > 0) wrttxt(iprt);
    std::vector<double> q2(numat + 2, 0.0);
    if (mode != 1 && !p.empty()) {
        chrge(p, q2);
        for (int i = 1; i <= numat; ++i) q2[i] = parameters_C::tore[nat[i]] - q2[i];
        if (in_house_only)
            for (int i = 1; i <= numat; ++i)
                if (std::fabs(q2[i]) > 0.999) q2[i] = 0.0;
    } else {
        for (int i = 1; i <= numat; ++i) q2[i] = 0.0;
    }
    int n = 1, ia = loc[1][1], ii = 0;
    for (int i = 1; i <= natoms; ++i) {
        std::string q[4] = { "  ", "  ", "  ", "  " };
        for (int j = 1; j <= 3; ++j) {
            q[j] = flag0;
            if (ia != i) continue;
            if (j != loc[2][n]) continue;
            q[j] = flag1;
            n = n + 1;
            ia = loc[1][n];
        }
        double w, x;
        if (na[i] > 0) {
            w = geo[2][i] * degree;
            x = geo[3][i] * degree;
            w = w - std::floor(w / 360.0) * 360.0;
            if (w < -1.e-6) w = w + 360.0;
            if (w > 180.000001) { x = x + 180.0; w = 360.0 - w; }
            double sgn = (x >= 0 ? 1.0 : -1.0);
            x = x - std::floor(x / 360.0 + sgn * (0.5 - 1.e-9) - 1.e-9) * 360.0;
        } else {
            w = geo[2][i];
            x = geo[3][i];
        }
        if (maps_C::latom == i) q[maps_C::lparam] = flagn;
        if (!gui && maps_C::latom1 == i) q[maps_C::lpara1] = flagn;
        if (maps_C::latom2 == i) q[maps_C::lpara2] = flagn;
        std::string blank_s = elemnt(labels[i]);
        if (labels[i] != 99 && labels[i] != 107) ii++;
        if (labels[i] == 1 && ii > 0) {
            if (std::fabs(atmass[ii] - 2.014) < 1.e-3) blank_s = " D";
            if (std::fabs(atmass[ii] - 3.016) < 1.e-3) blank_s = " T";
        }
        int k = 4;
        int igui = -10;
        if (maxtxt > 0) {
            line = txtatm[i];
            if (gui && txtatm[i].size() >= 2 && txtatm[i].substr(2) != " ") {
                std::string tt = txtatm[i];
                while (!tt.empty() && tt.back() == ' ') tt.pop_back();
                blank_s = blank_s + "(" + tt.substr(2) + "  ";
                igui = -10;
            } else {
                if (line != " ") blank_s = blank_s + "(" + line.substr(0, store_maxtxt) + ")  ";
                igui = 0;
            }
        }
        k = 0;
        if (mode != 1 && ii > 0 && nat[1] != 0) {
            if (std::fabs(parameters_C::ams[nat[ii]] - atmass[ii]) > 1.e-6) {
                int jj = (int)std::log10(atmass[ii]);
                if (jj < 0) jj = 0;
                char fmtbuf[16];
                std::snprintf(fmtbuf, sizeof(fmtbuf), "(f8.%d)", 6 - jj);
                if (blank_s.substr(0, 2) != " D" && blank_s.substr(0, 2) != " T") {
                    char ln[80];
                    std::snprintf(ln, sizeof(ln), fmtbuf, atmass[ii]);
                    std::string ls = ln;
                    for (int m = (int)ls.size() - 1; m >= 0; --m) {
                        if (ls[m] != '0') break;
                        ls[m] = ' ';
                    }
                    ls = ls + "0";
                    // strip leading spaces for the isotope label
                    size_t nz = ls.find_first_not_of(' ');
                    ls = (nz == std::string::npos) ? ls : ls.substr(nz);
                    size_t jp = blank_s.find('(');
                    if (jp != std::string::npos)
                        blank_s = blank_s.substr(0, jp) + ls + blank_s.substr(jp);
                    else
                        blank_s = blank_s.substr(0, 2) + ls;
                }
            }
            int jj = std::max(4, maxtxt + 2);
            if (isotopes) jj = jj + 8;
            k = std::max(0, 8 - jj);
        } else {
            if (mode != 1) { /* j = max(4, maxtxt+2) */ }
            else { /* j = max(9, maxtxt+3) */ }
        }
        if (labels[i] == 0) continue;
        if (!l_atom.empty() && !l_atom[i]) continue;   // unset l_atom => all atoms printable
        std::string line_extra;
        char ch[128];
        if (na[i] == igui || cart) {
            if (mode != 1) {
                if (labels[i] != 99 && labels[i] != 107) {
                    if (charge) std::snprintf(ch, sizeof(ch), "%8.4f", q2[ii]);
                    else std::snprintf(ch, sizeof(ch), "%8.4f", 0.0);
                    line_extra = ch;
                    std::snprintf(ch, sizeof(ch), "%s %s %s %s %s %s %s %s\n",
                        blank_s.c_str(), fixed(geo[1][i], fmt1.c_str()).c_str(), q[1].c_str(),
                        fixed(w, fmt23.c_str()).c_str(), q[2].c_str(),
                        fixed(x, fmt23.c_str()).c_str(), q[3].c_str(), line_extra.c_str());
                    std::fprintf(stdout, "%s", ch);
                } else {
                    std::snprintf(ch, sizeof(ch), "%s %s %s %s %s %s %s %s\n",
                        blank_s.c_str(), fixed(geo[1][i], fmt1.c_str()).c_str(), q[1].c_str(),
                        fixed(w, fmt23.c_str()).c_str(), q[2].c_str(),
                        fixed(x, fmt23.c_str()).c_str(), q[3].c_str(), " ");
                    std::fprintf(stdout, "%s", ch);
                }
            } else {
                std::snprintf(ch, sizeof(ch), "%6d      %s %s %s %s %s %s %s\n", i,
                    blank_s.c_str(), fixed(geo[1][i], fmt1.c_str()).c_str(), q[1].c_str(),
                    fixed(w, fmt23.c_str()).c_str(), q[2].c_str(),
                    fixed(x, fmt23.c_str()).c_str(), q[3].c_str());
                std::fprintf(stdout, "%s", ch);
            }
        } else {
            if (mode != 1) {
                if (labels[i] != 99 && labels[i] != 107) {
                    if (charge) std::snprintf(ch, sizeof(ch), "%8.4f", q2[ii]);
                    else std::snprintf(ch, sizeof(ch), "%8.4f", 0.0);
                    line_extra = ch;
                    std::snprintf(ch, sizeof(ch), "%s %s %s %s %s %s %s %3d %3d %3d %s\n",
                        blank_s.c_str(), fixed(geo[1][i], fmt1.c_str()).c_str(), q[1].c_str(),
                        fixed(w, fmt23.c_str()).c_str(), q[2].c_str(),
                        fixed(x, fmt23.c_str()).c_str(), q[3].c_str(),
                        na[i], nb[i], nc[i], line_extra.c_str());
                    std::fprintf(stdout, "%s", ch);
                } else {
                    std::snprintf(ch, sizeof(ch), "%s %s %s %s %s %s %s %3d %3d %3d\n",
                        blank_s.c_str(), fixed(geo[1][i], fmt1.c_str()).c_str(), q[1].c_str(),
                        fixed(w, fmt23.c_str()).c_str(), q[2].c_str(),
                        fixed(x, fmt23.c_str()).c_str(), q[3].c_str(),
                        na[i], nb[i], nc[i]);
                    std::fprintf(stdout, "%s", ch);
                }
            } else {
                std::snprintf(ch, sizeof(ch), "%6d      %s %s %s %s %s %s %s %3d %3d %3d\n", i,
                    blank_s.c_str(), fixed(geo[1][i], fmt1.c_str()).c_str(), q[1].c_str(),
                    fixed(w, fmt23.c_str()).c_str(), q[2].c_str(),
                    fixed(x, fmt23.c_str()).c_str(), q[3].c_str(),
                    na[i], nb[i], nc[i]);
                std::fprintf(stdout, "%s", ch);
            }
        }
    }
    maxtxt = store_maxtxt;
    if (mode == 1) return;
    std::fprintf(stdout, "\n");
    if (ndep != 0) {
        int nn = 1;
        int i = 1;
        for (;;) {
            int j = i;
            for (;;) {
                if (j == ndep) break;
                if (locpar[j] != locpar[j + 1] || idepfn[j] != idepfn[j + 1] || j - i >= 9) break;
                if (idepfn[i] == 18 || idepfn[i] == 19) {
                    if (std::fabs(depmul[nn] - depmul[nn + 1]) > 1.e-10) break;
                    nn++;
                }
                j++;
            }
            if (idepfn[i] == 18 || idepfn[i] == 19) {
                std::fprintf(stdout, "%4d%3d%13.9f", locpar[i], idepfn[i], depmul[nn]);
                for (int k = i; k <= j; ++k) std::fprintf(stdout, "%5d", locdep[k]);
                std::fprintf(stdout, "\n");
                nn++;
            } else {
                std::fprintf(stdout, "%4d%3d", locpar[i], idepfn[i]);
                for (int k = i; k <= j; ++k) std::fprintf(stdout, "%5d", locdep[k]);
                std::fprintf(stdout, "\n");
            }
            i = j + 1;
            if (j == ndep) break;
            if (i > ndep) break;
        }
        if (i <= ndep) {
            if (idepfn[i] == 18 || idepfn[i] == 19) {
                std::fprintf(stdout, "%4d%3d%13.9f", locpar[i], idepfn[i], depmul[nn]);
                for (int k = i; k <= ndep; ++k) std::fprintf(stdout, "%5d", locdep[k]);
                std::fprintf(stdout, "\n");
            } else {
                std::fprintf(stdout, "%4d%3d", locpar[i], idepfn[i]);
                for (int k = i; k <= ndep; ++k) std::fprintf(stdout, "%5d", locdep[k]);
                std::fprintf(stdout, "\n");
            }
        }
        std::fprintf(stdout, "\n");
    }
    if (keywrd.find(" SETPI") != std::string::npos && !pibonds_txt.empty()) {
        for (int n = 1; n <= 10; ++n) {
            if (pibonds_txt[n] == " ") break;
            if (n == 1) {
                std::string s = pibonds_txt[n];
                while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) s.pop_back();
                std::fprintf(stdout, "%s  User-supplied pi bonds\n", s.c_str());
            } else {
                std::string s = pibonds_txt[n];
                while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) s.pop_back();
                std::fprintf(stdout, "%s\n", s.c_str());
            }
        }
    }
}
