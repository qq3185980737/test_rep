// getgeg.cpp — C++ translation of MOPAC 2016 "getgeg.F90".
// Reads a geometry dataset (MOPAC input format) with symbolic parameters.
#include "getgeg.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "symmetry_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace molkst_C {
extern int natoms, nvar, ndep, ltxt;
extern std::string line;
}
namespace parameters_C {
extern std::vector<double> ams;
}
namespace symmetry_C {
extern std::vector<int> locpar, idepfn, locdep;
}
namespace common_arrays_C {
extern std::vector<double> atmass, xparam;
extern std::vector<std::string> simbol, txtatm;
extern std::vector<std::vector<int>> loc;
}
namespace chanel_C {
extern int iw;
}

extern double reada(const std::string&, int);
extern void getval(const std::string&, double&, std::string&);
extern void mopend(const std::string&);

static const char* elemnt_syms[108] = {
    " H","HE","LI","BE"," B"," C"," N"," O"," F","NE","NA","MG","AL","SI"," P"," S","CL","AR",
    "K ","CA","SC","TI"," V","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR",
    "RB","SR"," Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE"," I","XE",
    "CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF",
    "TA"," W","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
    "PA"," U","NP","PU","AM","CM","BK","CF","XX","FM","MD","CB","++"," +","--"," -","TV",""
};

void getgeg(int iread, std::vector<int>& labels, std::vector<std::vector<double>>& geo,
            std::vector<std::vector<int>>& lopt, std::vector<int>& na,
            std::vector<int>& nb, std::vector<int>& nc) {
    int istart[21] = {0};
    std::vector<std::vector<int>> lgeo(4, std::vector<int>(10001, -1));
    std::vector<std::vector<std::string>> tgeo(4, std::vector<std::string>(10001, " "));
    int nerr = 0, numat = 0;
    na[1] = 0; nb[1] = 0; nc[1] = 0;
    nb[2] = 0; nc[2] = 0; nc[3] = 0;
    molkst_C::ltxt = 0;
    molkst_C::natoms = 1;
    for (; molkst_C::natoms <= 900; ++molkst_C::natoms) {
        geo[1][molkst_C::natoms] = 0.0; geo[2][molkst_C::natoms] = 0.0; geo[3][molkst_C::natoms] = 0.0;
        std::getline(std::cin, molkst_C::line);
        if (molkst_C::line == " ") break;
        // Text label "(...)"
        size_t ip = molkst_C::line.find('(');
        if (ip != std::string::npos) {
            size_t kp = molkst_C::line.find(')', ip);
            common_arrays_C::txtatm[molkst_C::natoms] = molkst_C::line.substr(ip, kp - ip + 1);
            molkst_C::ltxt = std::max(molkst_C::ltxt, (int)(kp - ip + 1));
            std::string s = molkst_C::line.substr(0, ip) + molkst_C::line.substr(kp + 1);
            molkst_C::line = s;
        } else {
            common_arrays_C::txtatm[molkst_C::natoms] = " ";
        }
        // Uppercase.
        for (char& ch : molkst_C::line)
            if (ch >= 'a' && ch <= 'z') ch = (char)(ch - 'a' + 'A');
        // Tokenise.
        int nvalue = 0;
        bool leadsp = true;
        for (int k = 1; k <= 80; ++k) {
            char ch = (k <= (int)molkst_C::line.size()) ? molkst_C::line[k - 1] : ' ';
            if (leadsp && ch != ' ') { nvalue++; istart[nvalue] = k; }
            leadsp = (ch == ' ');
        }
        int j = 107;
        std::string tok = " " + molkst_C::line.substr(istart[1] - 1, 3);
        bool found = false;
        for (j = 1; j <= 107; ++j) {
            if (tok.find(elemnt_syms[j - 1] + std::string(" ")) != std::string::npos ||
                tok == " " + std::string(elemnt_syms[j - 1])) { found = true; break; }
        }
        if (!found && tok.find(" X") != std::string::npos) { j = 99; }
        if (!found && j == 107 && tok.find(" X") == std::string::npos) {
            std::fprintf(stdout, " ELEMENT NOT RECOGNIZED: %s\n",
                         molkst_C::line.substr(istart[1] - 1, 3).c_str());
            nerr++;
        }
        labels[molkst_C::natoms] = j;
        if (j != 99 && j < 107) {
            numat++;
            common_arrays_C::atmass[numat] = parameters_C::ams[j];
        }
        tgeo[1][molkst_C::natoms] = " ";
        tgeo[2][molkst_C::natoms] = " ";
        tgeo[3][molkst_C::natoms] = " ";
        if (molkst_C::natoms == 1) continue;
        na[molkst_C::natoms] = (int)std::lround(reada(molkst_C::line, istart[2]));
        getval(molkst_C::line.substr(istart[3] - 1), geo[1][molkst_C::natoms], tgeo[1][molkst_C::natoms]);
        if (molkst_C::natoms == 2) continue;
        if (istart[4] == 0) continue;
        nb[molkst_C::natoms] = (int)std::lround(reada(molkst_C::line, istart[4]));
        getval(molkst_C::line.substr(istart[5] - 1), geo[2][molkst_C::natoms], tgeo[2][molkst_C::natoms]);
        if (molkst_C::natoms == 3) continue;
        nc[molkst_C::natoms] = (int)std::lround(reada(molkst_C::line, istart[6]));
        getval(molkst_C::line.substr(istart[7] - 1), geo[3][molkst_C::natoms], tgeo[3][molkst_C::natoms]);
    }
    molkst_C::natoms = molkst_C::natoms - 1;
    for (int k = 1; k <= molkst_C::natoms; ++k) { lgeo[1][k] = -1; lgeo[2][k] = -1; lgeo[3][k] = -1; }
    int ivar = -1;
    molkst_C::nvar = 0;
    molkst_C::ndep = 0;
    int kerr = 0;
    for (;;) {
        std::getline(std::cin, molkst_C::line);
        if (molkst_C::line == " ") {
            if (ivar == -1) {
                int merr = 0;
                for (int k = 1; k <= molkst_C::natoms; ++k)
                    for (int q = 1; q <= 3; ++q)
                        if (geo[q][k] < -998) merr++;
                if (merr == 0) break;
                ivar = molkst_C::nvar;
                continue;
            } else break;
        }
        for (char& ch : molkst_C::line)
            if (ch >= 'a' && ch <= 'z') ch = (char)(ch - 'a' + 'A');
        int p = 1;
        while (p <= 80 && (p > (int)molkst_C::line.size() || molkst_C::line[p - 1] == ' ')) p++;
        int l = p;
        while (l <= p + 12 && (l <= (int)molkst_C::line.size() && molkst_C::line[l - 1] != ' ')) l++;
        double sum = reada(molkst_C::line, l);
        int n = 0, lerr = 0;
        for (int j = 1; j <= molkst_C::natoms; ++j) {
            for (int k = 1; k <= 3; ++k) {
                bool match = (tgeo[k][j] == molkst_C::line.substr(p - 1, l - p));
                if (!match && tgeo[k][j].size() >= 2 &&
                    tgeo[k][j].substr(1) == molkst_C::line.substr(p - 1, l - p) &&
                    tgeo[k][j][0] == '-') match = true;
                if (!match) continue;
                if (lgeo[k][j] != -1) lerr = 1;
                lgeo[k][j] = lgeo[k][j] + 1;
                n++;
                geo[k][j] = sum;
                if (n == 1) {
                    molkst_C::nvar++;
                    common_arrays_C::loc[1][molkst_C::nvar] = j;
                    common_arrays_C::loc[2][molkst_C::nvar] = k;
                    common_arrays_C::xparam[molkst_C::nvar] = sum;
                    common_arrays_C::simbol[molkst_C::nvar] = tgeo[k][j];
                    if (common_arrays_C::simbol[molkst_C::nvar][0] == '-') {
                        std::fprintf(stdout, " NEGATIVE SYMBOLICS MUST BE PRECEEDED BY  THE POSITIVE EQUIVALENT\n");
                        std::fprintf(stdout, " FAULTY SYMBOLIC:  %s\n", common_arrays_C::simbol[molkst_C::nvar].c_str());
                        mopend("NEGATIVE SYMBOLICS MUST BE PRECEEDED BY THE POSITIVE EQUIVALENT");
                        return;
                    }
                }
                if (n <= 1) continue;
                molkst_C::ndep++;
                symmetry_C::locpar[molkst_C::ndep] = common_arrays_C::loc[1][molkst_C::nvar];
                symmetry_C::idepfn[molkst_C::ndep] = common_arrays_C::loc[2][molkst_C::nvar];
                if (tgeo[k][j][0] == '-') {
                    symmetry_C::idepfn[molkst_C::ndep] = 14;
                    if (common_arrays_C::loc[2][molkst_C::nvar] != 3) {
                        kerr++;
                        std::fprintf(stdout, " ONLY DIHEDRAL SYMBOLICS  CAN BE PRECEEDED BY A \"-\" SIGN\n");
                    }
                }
                symmetry_C::locdep[molkst_C::ndep] = j;
            }
        }
        kerr += lerr;
        if (lerr == 1) {
            std::fprintf(stdout, " THE FOLLOWING SYMBOL HAS BEEN DEFINED MORE THAN ONCE: %s\n",
                         molkst_C::line.substr(p - 1, l - p).c_str());
            nerr++;
        }
        if (n == 0) {
            std::fprintf(stdout, " THE FOLLOWING SYMBOLIC WAS NOT USED: %s\n",
                         molkst_C::line.substr(p - 1, l - p).c_str());
            nerr++;
        }
    }
    int merr = 0;
    for (int k = 1; k <= molkst_C::natoms; ++k)
        for (int q = 1; q <= 3; ++q)
            if (geo[q][k] < -998) merr++;
    if (merr != 0) std::fprintf(stdout, "%4d GEOMETRY VARIABLES WERE NOT DEFINED\n", merr);
    if (merr + kerr + nerr != 0) {
        std::fprintf(stdout, " THE GEOMETRY DATA-SET CONTAINED%3d ERRORS\n", merr + kerr + nerr);
        mopend("THE GEOMETRY DATA-SET CONTAINED ERRORS");
        return;
    }
    // Sort optimisation parameters into increasing atom order.
    if (ivar != -1) molkst_C::nvar = ivar;
    for (int i = 1; i <= molkst_C::nvar; ++i) {
        int j = 100000, k = i;
        for (int l = i; l <= molkst_C::nvar; ++l) {
            if (j <= common_arrays_C::loc[1][l] * 4 + common_arrays_C::loc[2][l]) continue;
            k = l;
            j = common_arrays_C::loc[1][l] * 4 + common_arrays_C::loc[2][l];
        }
        std::string st = common_arrays_C::simbol[i];
        common_arrays_C::simbol[i] = common_arrays_C::simbol[k];
        common_arrays_C::simbol[k] = st;
        double s = common_arrays_C::xparam[i];
        common_arrays_C::xparam[i] = common_arrays_C::xparam[k];
        common_arrays_C::xparam[k] = s;
        for (int q = 1; q <= 2; ++q) {
            int l = common_arrays_C::loc[q][i];
            common_arrays_C::loc[q][i] = common_arrays_C::loc[q][k];
            common_arrays_C::loc[q][k] = l;
        }
    }
    for (int l = molkst_C::nvar + 1; l <= (int)common_arrays_C::loc[1].size() - 1; ++l) {
        common_arrays_C::loc[1][l] = 0;
        common_arrays_C::loc[2][l] = 0;
    }
    // Convert to radians.
    double degree = 1.7453292519943e-2;
    for (int i = 1; i <= molkst_C::nvar; ++i)
        if (common_arrays_C::loc[2][i] != 1) common_arrays_C::xparam[i] *= degree;
    for (int k = 1; k <= molkst_C::natoms; ++k) {
        geo[2][k] *= degree;
        geo[3][k] *= degree;
    }
    molkst_C::ltxt = molkst_C::ltxt;
    for (int k = 1; k <= molkst_C::natoms; ++k) { lopt[1][k] = 0; lopt[2][k] = 0; lopt[3][k] = 0; }
    for (int i = 1; i <= molkst_C::nvar; ++i)
        lopt[common_arrays_C::loc[2][i]][common_arrays_C::loc[1][i]] = 1;
}
