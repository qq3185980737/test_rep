// getpdb.cpp — C++ translation of MOPAC 2016 "getpdb.F90".
// Reads a Brookhaven Protein Data Bank geometry from stdin: ATOM/HETATM
// records, element identification (77-78 / 13-14 fields), residue/chain
// bookkeeping (START_RES / CHAINS keywords), alternative location filters,
// and optional conversion to internal coordinates.
#include "getpdb.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "parameters_C.h"
#include "elemts_C.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace molkst_C {
extern int maxtxt, natoms, numat, na1, nbreaks, ncomments;
extern std::string keywrd, line;
}
namespace common_arrays_C {
extern std::vector<std::string> txtatm, txtatm1, all_comments;
extern std::vector<double> atmass;
extern std::vector<int> na, nb, nc, labels, breaks;
extern std::vector<std::vector<int>> lopt;
extern std::vector<std::vector<double>> coord, break_coords;
}
namespace chanel_C {
extern int iw, ir;
}
namespace parameters_C {
extern std::vector<double> ams;
}
namespace elemts_C {
extern std::vector<std::string> elemnt;
}

extern void upcase(std::string&, int);
extern double reada(const std::string&, int);
extern void mopend(const std::string&);
extern void web_message(int, const char*);
extern void xyzint(const std::vector<double>&, int, const std::vector<int>&,
                   const std::vector<int>&, const std::vector<int>&, double, std::vector<double>&);
extern const char* elemnt(int);

void getpdb(std::vector<std::vector<double>>& geo) {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace chanel_C;
    const int maxel = 100;
    std::vector<std::string> element(maxel + 2, "  ");
    std::vector<int> ielem(maxel + 2, 0);
    int i;
    for (i = 1; i <= 83; ++i) {
        std::string en = elemnt(i);
        if (en[0] == ' ')
            element[i] = std::string(1, en[1]) + std::string(1, en[0]);
        else {
            char c2 = en[1];
            if (c2 >= 'a' && c2 <= 'z') c2 = (char)(c2 - 'a' + 'A');
            element[i] = std::string(1, en[0]) + std::string(1, c2);
        }
        ielem[i] = i;
    }
    element[i] = "X "; ielem[i] = 1;
    i = i + 1;
    element[i] = "L "; ielem[i] = 0;
    i = i + 1;
    element[i] = "D "; ielem[i] = 0;
    int defined_elements = i;
    int new_elements = 0;
    std::string letter;
    size_t p0 = keywrd.find(" CHAINS");
    if (p0 != std::string::npos) {
        size_t p = keywrd.find('(', p0);
        size_t q = keywrd.find(") ", p0);
        if (p != std::string::npos && q != std::string::npos && q > p + 1)
            letter = keywrd.substr(p + 1, q - p - 2);
        else letter = " ";
    } else letter = " ";
    bool l_pdb = (keywrd.find(" PDB ") != std::string::npos || keywrd.find(" PDB(") != std::string::npos);
    if (keywrd.find(" PDB(") != std::string::npos) {
        // PDB(<symbol>:<Z>,...) user-defined element mapping.
        size_t ip = keywrd.find(" PDB(");
        ip = ip + 5;
        size_t jp = keywrd.find(')', ip);
        if (jp != std::string::npos) {
            int k = defined_elements + 1;
            new_elements = 1;
            for (;;) {
                element[k] = keywrd.substr(ip, jp - ip);
                int ii = (int)element[k].find(':');
                if (ii != 0) element[k] = element[k].substr(0, 2);
                size_t ic = keywrd.find(':', ip);
                if (ic == ip) break;
                ielem[k] = (int)(reada(keywrd, (int)ic) + 0.5);
                size_t icom = keywrd.find(',', ic);
                if (icom != std::string::npos && icom != ic) {
                    ip = icom + 1;
                    k++;
                    new_elements++;
                } else {
                    break;
                }
            }
        }
        std::fprintf(stdout, " KEYWORD PDB IS INCORRECT\n");
        mopend("KEYWORD PDB IS INCORRECT");
    }
    char typea = 'A', typer = 'A';
    size_t pa = keywrd.find(" ALT_A=");
    if (pa != std::string::npos) typea = keywrd[pa + 7];
    pa = keywrd.find(" ALT_R=");
    if (pa != std::string::npos) typer = keywrd[pa + 7];
    maxtxt = 0;
    natoms = 0;
    numat = 0;
    int npdb = 0;
    maxtxt = 26;
    int nline = 0;
    int previous_res = 2000;
    ncomments = 0;
    nbreaks = 0;
    std::string new_key_1 = " ", new_key_2 = " ";
    std::string ch = " ";
    bool lchain = true;
    std::vector<std::string> tmp_comments(10000);
    int n_water = 0, old_natoms = 0;
    bool last_atom = true, first = true;
    std::vector<std::string> txtpdb(20, "    ");
    std::vector<int> ntxt_loc(20, 2);
    int nvalue = 0;
    std::string line1;
    bool done = false;
    while (!done) {
        if (!std::getline(std::cin, line)) break;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        nline++;
        if (natoms > 0 && line == " ") break;
        bool is_atom = (line.size() >= 4 && line.substr(0, 4) == "ATOM");
        bool is_het = (line.size() >= 6 && line.substr(0, 6) == "HETATM");
        if (!is_atom && !is_het) {
            const char* tags[] = { "ATOM  ", "HETATM", "TITLE ", "HEADER", "COMPND", "SOURCE",
                "KEYWDS", "HELIX ", "SHEET ", "REMARK", "USER  ", "EXPDTA", "AUTHOR", "REVDAT",
                "JRNL  ", "DBREF ", "SEQRES", "HET   ", "HETNAM", "LINK  ", "CRYST1", "SCALE",
                "ORIGX", "FORMUL", "SEQRES" };
            bool found = false;
            for (const char* t : tags) {
                if (line.find(t) != std::string::npos) { found = true; break; }
            }
            if (found) {
                ncomments++;
                if ((int)tmp_comments.size() < ncomments) tmp_comments.resize(ncomments + 100);
                tmp_comments[ncomments] = "*" + line.substr(0, 80);
                continue;
            }
            if (line.size() >= 3 && line.substr(0, 3) == "TER" && natoms > 0) {
                nbreaks++;
                if (nbreaks > 400) {
                    std::fprintf(stdout, "Maximum number (400) of breaks (lines starting with TER) exceeded\n");
                    mopend("Maximum number (400) of breaks (lines starting with TER) exceeded");
                    return;
                }
                line1 = txtatm[natoms].substr(21);
                int j = natoms;
                for (i = natoms; i >= 1; --i) {
                    j = i;
                    if (labels[i] != 1) break;
                    if (txtatm[i].substr(21, 5) != line1.substr(0, 5)) { j = j + 1; break; }
                }
                for (int c = 1; c <= 3; ++c) break_coords[c - 1][nbreaks] = geo[c][j];
                breaks[nbreaks] = j;
                lchain = true;
            }
            if (line.size() >= 3 && line.substr(0, 3) == "END") break;
            if (l_pdb) {
                if (first) {
                    first = false;
                    std::fprintf(stdout, "\n          Keyword PDB present, but the following lines are not in PDB format\n");
                }
                std::fprintf(stdout, "%s\n", line.c_str());
            }
            continue;
        }
        if (letter != " " && line.size() >= 22 && letter.find(line[21]) == std::string::npos) continue;
        if (line.size() >= 17 && line[16] != ' ') {
            if (typea == ' ') {
                std::fprintf(stdout, " The data set contains alternative location indicators.  Keyword ALT_A must be used\n");
                web_message(iw, "ALT_A.html");
                std::fprintf(stdout, "%d lines have been read in\nFaulty line: %s\n", nline, line.c_str());
                mopend("The data set contains alternative location indicators.  Keyword ALT_A must be used");
                return;
            }
            if (line[16] != typea) continue;
        }
        if (line.size() >= 27 && line[26] != ' ') {
            if (typer == '9') {
                std::fprintf(stdout, " The data set contains alternative location indicators.  Keyword ALT_R must be used\n");
                web_message(iw, "ALT_A.html");
                std::fprintf(stdout, "%d lines have been read in\nFaulty line: %s\n", nline, line.c_str());
                mopend("The data set contains alternative location indicators.  Keyword ALT_R must be used");
                return;
            }
            if (line[26] != typer) continue;
        }
        natoms++;
        std::string up = line;
        upcase(up, 80);
        line = up;
        // Commas -> spaces.
        for (char& c : line) if (c == ',') c = ' ';
        if (line.size() >= 14 && (line[13] == 'O' || line[13] == 'H')) n_water++;
        if (line.size() < 14 || line[13] != 'H') {
            int current_res = (int)(reada(line, 23) + 0.5);
            bool het_chain = (line.size() >= 6 && line.substr(0, 6) == "HETATM" && line[21] == ' ' && ch == " ");
            if (lchain && (!het_chain || natoms < 100 || last_atom)) {
                lchain = false;
                i = 23;
                while (i <= 25 && i < (int)line.size() && line[i - 1] == ' ') ++i;
                if (new_key_1 != " ") {
                    new_key_1 = new_key_1 + " " + line.substr(i - 1, 4) + line.substr(21, 1);
                    new_key_2 = new_key_2 + line.substr(21, 1);
                    previous_res = 2000;
                } else {
                    new_key_1 = " START_RES=(" + line.substr(i - 1, 4) + line.substr(21, 1);
                    new_key_2 = " CHAINS=(" + line.substr(21, 1);
                }
            } else {
                if (current_res - previous_res > 1 || (natoms > 2 && line.size() >= 22 && line[21] != ch[0])) {
                    if (line.size() < 6 || line.substr(0, 6) != "HETATM" || line.size() < 22 || line[21] != ' ') {
                        i = 23;
                        while (i <= 25 && i < (int)line.size() && line[i - 1] == ' ') ++i;
                        if (line[21] != ch[0]) {
                            if (n_water != natoms - old_natoms) new_key_1 = new_key_1 + " " + line.substr(i - 1, 4) + line.substr(21, 1);
                            n_water = 0;
                            old_natoms = natoms;
                            new_key_2 = new_key_2 + line.substr(21, 1);
                        } else {
                            if (n_water != natoms - old_natoms) new_key_1 = new_key_1 + "-" + line.substr(i - 1, 4);
                            n_water = 0;
                            old_natoms = natoms;
                        }
                    }
                }
            }
            previous_res = current_res;
            ch = line.size() >= 22 ? std::string(1, line[21]) : " ";
        }
        last_atom = (line.size() >= 4 && line.substr(0, 4) == "ATOM");
        txtatm.resize(natoms + 2);
        txtatm1.resize(natoms + 2);
        txtatm[natoms] = line.substr(0, 26);
        txtatm1[natoms] = line.substr(0, 26);
        // Element identification from the 77-78 field, falling back to 13-14.
        std::string ele = "  ";
        if (line.size() >= 78) {
            ele = line.substr(76, 2);
            if (ele[0] != ' ' && ele[1] == ' ') ele = "  ";
            if (ele[0] == ' ' || (ele[0] >= '0' && ele[0] <= '9') || ele[0] == '.') {
                if (ele.size() >= 2) { ele[0] = ele[1]; ele[1] = ' '; }
                if (ele[0] == ' ' || (ele[0] >= '0' && ele[0] <= '9') || ele[0] == '.') {
                    ele = line.substr(12, 2);
                    if (line.substr(0, 4) == "ATOM" && ele[0] == 'H') {
                        ele[1] = ' ';
                    } else if (ele[0] >= '0' && ele[0] <= '9') {
                        ele[0] = ele[1];
                        ele[1] = ' ';
                    }
                    if (ele[1] >= '0' && ele[1] <= '9') ele[1] = ' ';
                }
            }
        } else {
            ele = line.substr(12, 2);
        }
        int label = 0;
        bool matched = false;
        for (i = defined_elements + 1; i <= defined_elements + new_elements; ++i) {
            if (ele == element[i]) { label = ielem[i]; matched = true; break; }
        }
        if (!matched) {
            for (i = 1; i <= defined_elements; ++i) {
                if (ele == element[i]) { label = ielem[i]; matched = true; break; }
            }
        }
        if (!matched) {
            natoms--;
            bool seen = false;
            for (i = 1; i <= npdb; ++i)
                if (ele == txtpdb[i]) { seen = true; break; }
            if (seen) continue;
            npdb++;
            txtpdb[npdb] = ele;
            if (ele[1] == ' ') ntxt_loc[npdb] = 1;
            std::fprintf(stdout, "  UNRECOGNIZED SPECIES: ('%s', on line %d\n",
                ele.substr(0, ntxt_loc[npdb]).c_str(), nline);
            continue;
        }
        if (label == 0) {
            natoms--;
        } else {
            if (label == 99) label = 1;
            numat++;
            labels.resize(natoms + 2);
            labels[natoms] = label;
            if (label != 99) {
                atmass.resize(numat + 2);
                atmass[numat] = parameters_C::ams[label];
            }
            geo[1][natoms] = reada(line.substr(30, 8), 1);
            geo[2][natoms] = reada(line.substr(38, 8), 1);
            geo[3][natoms] = reada(line.substr(46, 8), 1);
            lopt[1][natoms] = 1;
            lopt[2][natoms] = 1;
            lopt[3][natoms] = 1;
            na.resize(natoms + 2);
            nb.resize(natoms + 2);
            nc.resize(natoms + 2);
            na[natoms] = 0;
            nb[natoms] = 0;
            nc[natoms] = 0;
        }
    }
    if (npdb != 0) {
        std::fprintf(stdout, " THE SPECIES THAT WERE NOT RECOGNIZED CAN BE RECOGNIZED BY USING THE FOLLOWING KEYWORD\n");
        std::string kw = " PDB(";
        for (i = 1; i <= npdb; ++i) {
            kw += txtpdb[i].substr(0, ntxt_loc[i]) + ":Z";
            kw += (i == npdb ? ")" : ",");
        }
        std::fprintf(stdout, "%s\n", kw.c_str());
        std::fprintf(stdout, " WHERE 'Z' IS/ARE ATOMIC NUMBERS\n");
        mopend(" UNRECOGNIZED ELEMENT ");
    }
    for (i = 1; i <= numat; ++i) {
        if ((int)txtatm[i].size() >= 17) txtatm[i][16] = ' ';
        if ((int)txtatm1[i].size() >= 17) txtatm1[i][16] = ' ';
    }
    all_comments.assign(ncomments + 101, "");
    for (i = 1; i <= ncomments; ++i) all_comments[i] = tmp_comments[i];
    if (keywrd.find("START_RES") == std::string::npos && new_key_1.size() >= 12 && new_key_1[11] != '(') {
        // Compress repeated spaces out of the residue list.
        std::string nk;
        std::string src = new_key_1;
        int sl = (int)src.size();
        for (i = 2; i < sl; ++i) {
            if (i + 1 <= sl && src[i - 1] == ' ' && src[i] == ' ') continue;
            if (i - 1 >= 0 && src[i - 2] == '-' && src[i - 1] == ' ') continue;
            nk += src[i - 1];
        }
        if (nk.size() >= 12 && nk[11] == ' ') nk = nk.substr(0, 11) + nk.substr(12);
        keywrd = " " + nk + ") " + keywrd;
    }
    if (keywrd.find("CHAINS") == std::string::npos && (new_key_2.size() < 9 || new_key_2[8] != '('))
        keywrd = " " + new_key_2 + ") " + keywrd;
    bool lxyz = (keywrd.find(" INT") == std::string::npos);
    if (lxyz) {
        na1 = 99;
        return;
    }
    for (i = 1; i <= natoms; ++i)
        for (int j = 1; j <= 3; ++j)
            coord[j-1][i] = geo[j][i];
    double degree = 57.29577951308232;
    std::vector<double> out(3 * (natoms + 2), 0.0);
    std::vector<double> flat(3 * (natoms + 2), 0.0);
    for (i = 1; i <= natoms; ++i)
        for (int j = 1; j <= 3; ++j) flat[(i - 1) * 3 + (j - 1)] = geo[j][i];
    xyzint(flat, natoms, na, nb, nc, degree, out);
    for (i = 1; i <= natoms; ++i)
        for (int j = 1; j <= 3; ++j) geo[j][i] = out[(i - 1) * 3 + (j - 1)];
    for (int j = 1; j <= 3; ++j)
        for (i = j; i <= 3; ++i) lopt[j][i] = 0;
    if (natoms > 2) {
        if (std::fabs(geo[2][3] - 180.0) < 1.e-4 || std::fabs(geo[2][3]) < 1.e-4) {
            std::fprintf(stdout, " DUE TO PROGRAM BUG, THE FIRST THREE ATOMS MUST NOT LIE IN A STRAIGHT LINE.\n");
            mopend("FIRST 3 ATOMS MUST NOT LIE IN A LINE.");
            return;
        }
    }
    degree = 1.7453292519943e-2;
    for (int j = 2; j <= 3; ++j)
        for (i = 1; i <= natoms; ++i) geo[j][i] = geo[j][i] * degree;
    na1 = 0;
}
