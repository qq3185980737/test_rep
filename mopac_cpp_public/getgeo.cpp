// getgeo.cpp — C++ translation of MOPAC 2016 "getgeo.F90".
// Reads the geometry from a MOPAC input (or .ARC) file: element symbols or
// atomic numbers, internal or Cartesian coordinates, optimization flags,
// connectivity (na/nb/nc), atom labels, velocity vectors, and conversion
// between coordinate systems. Also implements txt_to_atom_no (lewis.F90),
// which resolves quoted atom labels.
#include "getgeo.h"
#include "getdat_lines_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "maps_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace molkst_C {
extern int natoms, numat, maxtxt, numcal, id;
extern bool moperr, units, Angstroms;
extern std::string keywrd, line;
extern double arc_hof_1, arc_hof_2;
}
namespace chanel_C {
extern int iw, ir;
extern std::string input_fn;
}
namespace common_arrays_C {
extern std::vector<double> atmass;
extern std::vector<std::string> simbol, txtatm;
extern std::vector<int> na_store, nat;
extern std::vector<char> l_atom;
}
namespace parameters_C {
extern std::vector<double> ams;
}
namespace maps_C {
extern std::vector<double> react;
}
namespace funcon_C {
extern double a0;
}

extern void upcase(std::string&, int);
extern double reada(const std::string&, int);
extern void mopend(const std::string&);
extern void geout(int);
extern void nuchar(char*, int, double*, int&);
extern void gmetry(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
extern void xyzint(double*, int, int*, int*, int*, double, double*);
extern void web_message(int, const char*);
extern void txt_to_atom_no(std::string&, int, bool, int&);

static const char* kElem[108] = {
    "", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F",
    "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC",
    "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS",
    "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH",
    "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA",
    "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM",
    "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL",
    "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP",
    "PU", "AM", "CM", "BK", "MI", "XX", "+T", "-T", "CB", "++", "+", "--",
    "-", "TV"};

static std::string trimr(const std::string& s) {
    size_t b = s.find_last_not_of(' ');
    return b == std::string::npos ? "" : s.substr(0, b + 1);
}
static int len_trim(const std::string& s) {
    size_t b = s.find_last_not_of(' ');
    return (b == std::string::npos) ? 0 : (int)(b + 1);
}

void getgeo(int iread, std::vector<int>& labels,
            std::vector<std::vector<double>>& geo,
            std::vector<std::vector<double>>& xyz,
            std::vector<std::vector<int>>& lopt,
            std::vector<int>& na, std::vector<int>& nb, std::vector<int>& nc,
            bool int_) {
    using namespace molkst_C;
    using namespace chanel_C;
    using namespace common_arrays_C;
    using namespace parameters_C;
    using namespace maps_C;
    using namespace funcon_C;
    (void)iread;
    bool ircdrc = ((keywrd.find(" IRC") != std::string::npos) +
                   (keywrd.find(" DRC") != std::string::npos)) != 0;
    if (keywrd.find(" FORCETS") != std::string::npos) ircdrc = false;
    int icapa = 'A', icapz = 'Z';
    bool lturn = true;
    bool saddle = (keywrd.find("SADDLE") != std::string::npos);
    bool lxyz = (keywrd.find(" XYZ") != std::string::npos || saddle);
    int_ = (keywrd.find(" INT ") != std::string::npos);
    bool velo = (keywrd.find(" VELO") != std::string::npos);
    bool lmop = (keywrd.find(" MOPAC") != std::string::npos);
    int max_atoms = (int)txtatm.size();
    maxtxt = 0;
    if ((int)simbol.size() < natoms * 3) {
        natoms = 0;
        numat = 0;
        return;
    }
    for (int i = 0; i < natoms * 3 && i < (int)simbol.size(); ++i) simbol[i] = "---------";
    natoms = 0;
    numat = 0;
    int iserr = 0;
    bool mini = (keywrd.find(" MINI") != std::string::npos);
    units = false;
    Angstroms = true;
    if (keywrd.find(" A0 ") != std::string::npos) {
        units = true;
        Angstroms = false;
    } else if (keywrd.find(" Ang") != std::string::npos) {
        units = true;
    }
    std::vector<int> istart(41, 0);
    int nvalue = 0;
    double weight = 0.0;
    std::string ele, string(120, ' ');
    int label = 0;
    double sum = 0.0;
    int ndmy = 0;
    bool solid = false;
    bool leadsp = true;
    std::string turn;
    bool done = false;
    while (!done) {
        if (getdat_line_idx >= getdat_lines.size()) { done = true; break; }
        if (natoms > (int)geo[1].size() - 2 || natoms > 100000) {
            done = true; break;
        }
        line = getdat_lines[getdat_line_idx++];
        if (!line.empty() && line.back() == '\r') line.pop_back();
        line.resize(120, ' ');  // Fortran reads a fixed 120-column record
        if (trimr(line) == "$coord" || trimr(line) == "$end") continue;   // Fortran char compare pads with blanks
        if (line.find_first_not_of(' ') == std::string::npos) {   // blank line (Fortran: LINE == ' '  auto-pads)
            if (natoms == 0 && numcal == 1) {
                // Is this an ARC file?
                numcal = 2;
                sum = 0.0;
                std::string l;
                while (getdat_line_idx < getdat_lines.size()) {
                    l = getdat_lines[getdat_line_idx++];
                    if (!l.empty() && l.back() == '\r') l.pop_back();
                    if (l.find("HEAT OF FORMATION") != std::string::npos) sum = reada(l, 20);
                    if (l.find("FINAL GEOMETRY OBTAINED") != std::string::npos) break;
                    if (l.find("GEOMETRY IN CARTESIAN COORDINATE") != std::string::npos) break;
                    if (l.find("GEOMETRY IN MOPAC Z-MATRIX") != std::string::npos) break;
                }
                std::fprintf(stderr, "The data set was defined as a MOPAC input data file\n");
                std::fprintf(stderr, "but it appears to an archive file\n");
                std::fprintf(stderr, "An attempt will be made to read it as an .ARC file\n");
                if (std::fabs(sum) > 1.e-4) {
                    if (std::fabs(arc_hof_1) < 1.e-4) arc_hof_1 = sum;
                    else arc_hof_2 = sum;
                }
                natoms = -3;
                return;
            }
            done = true;
            break;
        }
        int ltl = len_trim(line);
        for (int i = 0; i < ltl; ++i) {
            char kh = line[i];
            if (kh == ',' || kh == '\t') line[i] = ' ';
        }
        if (natoms == max_atoms) {
            std::fprintf(stdout, "\n          Maximum number of atoms exceeded\n");
            std::fprintf(stdout, "          Maximum allowed:%6d\n", max_atoms);
            mopend("Maximum number of atoms exceeded");
            return;
        }
        if (natoms < 10) {
            const char* tags[] = { "ATOM", "HETATM", "TITLE", "HEADER", "ANISOU", "COMPND",
                "SOURCE", "KEYWDS", "USER ", "HELIX", "SHEET", "REMARK", "AUTHOR", "REVDAT",
                "JRNL  ", "SEQRES" };
            bool pdblike = false;
            for (const char* t : tags)
                if (line.find(t) != std::string::npos) { pdblike = true; break; }
            if (pdblike) {
                natoms = -2;
                return;
            }
        }
        // Element labelled with text in parentheses.
        size_t ip = line.find('(');
        if (ip != std::string::npos) {
            size_t kp = line.find(')');
            if (kp == std::string::npos) {
                std::fprintf(stdout, " Atom%5d has an opening parenthesis but no closing parenthesis\n", natoms);
                std::fprintf(stdout, " Line :'%s'\n", line.c_str());
                std::fprintf(stdout, "\n GEOMETRY IS FAULTY.  GEOMETRY READ IN IS\n");
                atmass[1] = -1.0;
                nat[1] = 0;
                std::fill(l_atom.begin(), l_atom.end(), true);
                geout(iw);
                mopend("GEOMETRY IS FAULTY");
                return;
            }
            int jj = (line[kp] == ')') ? (int)kp : (int)(kp + 1);
            txtatm[natoms + 1] = line.substr(ip + 1, jj - ip - 1);
            if (maxtxt == 26 && (int)(kp - ip - 1) != 26) {
                if (keywrd.find(" GEO-OK") == std::string::npos) {
                    int i = natoms + 1;
                    char no[4];
                    std::snprintf(no, sizeof(no), "%d", (int)std::log10(i + 1.1) + 1);
                    std::fprintf(stdout, "\n          Atom label length of 26 detected, but atom %*d has a label of%3d characters.\n",
                                 no[0] - '0', natoms + 1, (int)(kp - ip - 1));
                    std::fprintf(stdout, "\nFaulty line = '%s'\n", trimr(line).c_str());
                    std::fprintf(stdout, "\n          If this is intended, add keyword \"GEO-OK\"\n");
                    char buf[80];
                    std::snprintf(buf, sizeof(buf), "Fault detected in atom number %*d", no[0] - '0', natoms + 1);
                    mopend(buf);
                    std::exit(1);
                }
            }
            maxtxt = std::max(maxtxt, (int)(kp - ip - 1));
            if (maxtxt > 38) {
                const char* tags[] = { "ATOM", "HETATM", "TITLE", "HEADER", "ANISOU", "COMPND",
                    "SOURCE", "KEYWDS", "USER ", "HELIX", "SHEET" };
                bool hit = false;
                for (const char* t : tags)
                    if (line.find(t) != std::string::npos) { hit = true; break; }
                if (hit) goto label70;
                std::fprintf(stdout, " Atom labels must not exceed 38 characters\n");
                std::fprintf(stdout, "                1         2         3         4\n");
                std::fprintf(stdout, "       1234567890123456789012345678901234567890\n");
                std::fprintf(stdout, "Line :\"%s\"\n", line.c_str());
                std::fprintf(stdout, "\n GEOMETRY IS FAULTY.  GEOMETRY READ IN IS\n");
                atmass[1] = -1.0;
                nat[1] = 0;
                geout(iw);
                mopend("GEOMETRY IS FAULTY");
                return;
            }
            line = line.substr(0, ip) + line.substr(kp + 1);
        } else {
            txtatm[natoms + 1] = " ";
        }
        upcase(line, ltl);
        for (int k = 0; k < 10; ++k) istart[k + 1] = len_trim(line) + 1;
        leadsp = true;
        nvalue = 0;
        int ii = 0;
        for (int jj = 1; jj <= ltl; ++jj) {
            ii = ii + 1;
            if (leadsp && line[ii - 1] != ' ') {
                nvalue = nvalue + 1;
                istart[nvalue] = ii;
                if (ii > 2) {
                    if (line[ii - 2] == '-' && line[ii - 3] != ' ') istart[nvalue] = ii - 1;
                }
                if (line[ii - 1] == '"') {
                    for (int j = 1; j <= 40; ++j) {
                        ii = ii + 1;
                        if (line[ii - 1] == '"') break;
                    }
                }
            }
            leadsp = (line[ii - 1] == ' ' ||
                      (ii != istart[std::max(nvalue, 1)] && line[ii - 1] == '-'));
            if (ii > 1) {
                if ((line[ii - 2] == 'D' && line[ii - 1] == '-') ||
                    (line[ii - 2] == 'E' && line[ii - 1] == '-')) leadsp = false;
            }
        }
        if (nvalue == 4) {
            // Four data: Cartesian coordinates without optimization flags?
            // Fortran index(line(istart(4):), " ") is 1-based; std::string::find is 0-based, so add 1.
            int i = (int)line.substr(istart[4] - 1).find(' ') + istart[4] - 1;
            int k = 0;
            if (i > istart[4] + 1) k = 1;
            for (int j = istart[4]; j <= i; ++j) {
                if (line[j - 1] < 'A' || line[j - 1] > 'Z') k = 1;
            }
            i = (int)line.substr(istart[1] - 1).find(' ') + istart[1] - 1;
            int j = istart[1];
            for (; j <= i; ++j) {
                if (line[j - 1] >= 'A' && line[j - 1] <= 'Z') break;
            }
            if (j > i && k == 0) {
                natoms = natoms + 1;
                geo[1][natoms] = reada(line, istart[1]);
                geo[2][natoms] = reada(line, istart[2]);
                geo[3][natoms] = reada(line, istart[3]);
                lopt[1][natoms] = 1;
                lopt[2][natoms] = 1;
                lopt[3][natoms] = 1;
                na[natoms] = 0;
                nb[natoms] = 0;
                nc[natoms] = 0;
                std::string el2 = line.substr(istart[4] - 1, 2);
                for (int ie = 1; ie <= 107; ++ie) {
                    std::string kk = kElem[ie];
                    if ((int)kk.size() < 2) kk += ' ';
                    if (el2 != kk) continue;
                    labels[natoms] = ie;
                    atmass[natoms] = ams[ie];
                    break;
                }
                if (!units) {
                    Angstroms = false;
                    units = true;
                }
                continue;
            }
        }
        // Establish element name and isotope.
        weight = 0.0;
        string = line.substr(istart[1] - 1, istart[2] - istart[1]);
        if (string == "+3") string = "+T";
        if (string == "-3") string = "-T";
        if (string[0] >= '0' && string[0] <= '9') {
            // Atomic number used: no isotope allowed.
            label = (int)(reada(string, 1) + 0.5);
            if (label == 0) { done = true; break; }
            if (label < 0 || label > 107) {
                std::fprintf(stdout, "  ILLEGAL ATOMIC NUMBER\n");
                goto label210;
            }
            goto label70;
        }
        // Atomic symbol used.
        double real = std::fabs(reada(string, 1));
        if (real < 1.e-15) {
            ele = string.substr(0, 2);
        } else {
            weight = real;
            if (string.size() >= 2 && string[1] >= '0' && string[1] <= '9') ele = string.substr(0, 1);
            else ele = string.substr(0, 2);
        }
        if (ele.size() >= 2 && ele[0] == '-' && ele[1] != '-') ele[1] = ' ';
        for (int ie = 1; ie <= 107; ++ie) {
            std::string kk = kElem[ie];
            if ((int)kk.size() < 2) kk += ' ';
            if (ele != kk) continue;
            label = ie;
            goto label70;
        }
        if (ele[0] == 'X') { label = 99; goto label70; }
        if (ele == "D ") { label = 1; weight = 2.014; goto label70; }
        else if (ele == "T ") { label = 1; weight = 3.016; goto label70; }
        else {
            const char* tags[] = { "ATOM", "HETATM", "TITLE", "HEADER", "ANISOU", "COMPND",
                "SOURCE", "KEYWDS", "HELIX", "SHEET", "REMARK", "USER ", "SEQRES" };
            bool hit = false;
            for (const char* t : tags)
                if (line.find(t) != std::string::npos) { hit = true; break; }
            if (hit) goto label70;
            if (trimr(line) == " *                    *") return;
            std::fprintf(stdout, "  UNRECOGNIZED ELEMENT NAME: ('%s')\n", ele.c_str());
            std::fprintf(stdout, "\n  Faulty line: \"%s\"\n\n", trimr(line).c_str());
            goto label210;
        }

    label70:
        natoms = natoms + 1;
        nb[natoms] = 0;
        nc[natoms] = 0;
        if (label != 99) numat = numat + 1;
        if (weight != 0.0) atmass[numat] = weight;
        else if (label != 99) atmass[numat] = ams[label];
        labels[natoms] = label;
        if (nvalue == 4) {
            geo[1][natoms] = reada(line, istart[2]);
            geo[2][natoms] = reada(line, istart[3]);
            geo[3][natoms] = reada(line, istart[4]);
            lopt[1][natoms] = 1;
            lopt[2][natoms] = 1;
            lopt[3][natoms] = 1;
        } else {
            geo[1][natoms] = reada(line, istart[2]);
            geo[2][natoms] = reada(line, istart[4]);
            geo[3][natoms] = reada(line, istart[6]);
            if (!mini && ircdrc) {
                turn = line.substr(istart[3] - 1, 1);
                if (turn == "T") {
                    lopt[1][natoms] = 1;
                    if (lturn) {
                        std::fprintf(stdout, " IN DRC MONITOR POTENTIAL ENERGY TURNING POINTS\n");
                        lturn = false;
                    }
                } else {
                    lopt[1][natoms] = 0;
                }
                turn = line.substr(istart[5] - 1, 1);
                lopt[2][natoms] = (turn == "T") ? 1 : 0;
                turn = line.substr(istart[7] - 1, 1);
                lopt[3][natoms] = (turn == "T") ? 1 : 0;
            } else {
                lopt[1][natoms] = (int)(reada(line, istart[3]) + 0.5);
                lopt[2][natoms] = (int)(reada(line, istart[5]) + 0.5);
                lopt[3][natoms] = (int)(reada(line, istart[7]) + 0.5);
                for (int i = 3; i <= 7; i += 2) {
                    char ch1 = line[istart[i] - 1];
                    if (!(ch1 >= icapa && ch1 <= icapz && natoms > 1)) continue;
                    iserr = 1;
                }
            }
        }
        if (line[istart[10] - 1] == '"') txt_to_atom_no(line, istart[10], false, ndmy);
        if (line[istart[9] - 1] == '"') txt_to_atom_no(line, istart[9], false, ndmy);
        if (line[istart[8] - 1] == '"') txt_to_atom_no(line, istart[8], false, ndmy);
        if (moperr) {
            std::fprintf(stdout, "          Error detected in definition of atom no.:%5d\n", natoms);
            std::fprintf(stdout, "\nText of atom: '%s'\n", line.c_str());
            return;
        }
        int j = len_trim(line);
        if (ltl != j) {
            int i = (int)line.substr(istart[8] - 1).find(' ') + istart[8];
            nvalue = 8;
            leadsp = true;
            for (; i <= j; ++i) {
                if (leadsp && line[i - 1] != ' ') {
                    nvalue = nvalue + 1;
                    istart[nvalue] = i;
                }
                leadsp = (line[i - 1] == ' ');
            }
        }
        sum = reada(line, istart[8]);
        int i = (int)line.substr(istart[8] - 1).find(' ') + istart[8];
        if (line.substr(istart[8] - 1, i - istart[8] + 1).find('.') != std::string::npos) sum = 0.0;
        na[natoms] = (int)(sum + 0.5);
        if (natoms == 1) na[1] = 0;
        if (lmop && natoms == 2) {
            na[2] = 1;
        } else if (lmop && natoms == 3 && na[3] == 0) {
            na[3] = 2;
            nb[3] = 1;
            geo[2][3] = geo[2][3] * 1.7453292519943e-2;
            geo[3][3] = 0.0;
            lopt[3][3] = 0;
        } else if (na[natoms] > 0) {
            nb[natoms] = (int)(reada(line, istart[9]) + 0.5);
            nc[natoms] = (int)(reada(line, istart[10]) + 0.5);
            geo[2][natoms] = geo[2][natoms] * 1.7453292519943e-2;
            geo[3][natoms] = geo[3][natoms] * 1.7453292519943e-2;
        }
        if (natoms == 3) {
            if (lopt[1][3] != 2 && lopt[2][3] != 2 && lopt[3][3] == 2) {
                na[3] = 1;
                nb[3] = 2;
                geo[3][3] = 0.0;
                lopt[3][3] = 0;
            } else if (lopt[3][3] == 1 && line.substr(istart[6] - 1, 2) == "2 ") {
                na[3] = 2;
                nb[3] = 1;
                geo[3][3] = 0.0;
                geo[2][3] = geo[2][3] * 1.7453292519943e-2;
                lopt[3][3] = 0;
            }
        }
        if (!mini) {
            if ((lopt[1][natoms] > 1 || lopt[2][natoms] > 1 || lopt[3][natoms] > 1) && natoms > 1) iserr = 1;
        } else {
            if (std::abs(lopt[1][natoms]) != std::abs(lopt[2][natoms]) ||
                std::abs(lopt[2][natoms]) != std::abs(lopt[3][natoms])) {
                std::fprintf(stdout, "When MINI is used, all coordinates for each atom must have the same flags\n");
                iserr = 1;
            }
        }
        if (iserr == 1) {
            const char* tags[] = { "ATOM", "HETATM", "TITLE", "HEADER", "ANISOU", "COMPND",
                "SOURCE", "KEYWDS", "HELIX", "SHEET", "REMARK", "USER ", "SEQRES" };
            bool hit = false;
            for (const char* t : tags)
                if (line.find(t) != std::string::npos) { hit = true; break; }
            if (hit) {
                natoms = -2;
                return;
            }
            // Must be Gaussian geometry input.
            for (int i = 2; i <= natoms; ++i) {
                for (int k = 1; k <= 3; ++k) {
                    int jj = (int)(geo[k][i] + 0.5);
                    if (std::fabs(geo[k][i] - jj) <= 1.e-5) continue;
                    std::fprintf(stdout, " GEOMETRY IS FAULTY.  GEOMETRY READ IN IS\n");
                    atmass[1] = -1.0;
                    nat[1] = 0;
                    geout(iw);
                    mopend("GEOMETRY IS FAULTY");
                    return;
                }
            }
            natoms = -1;
            return;
        }
    }
    // 120: all data read in.
    if (natoms == 0) {
        if (numcal == 1) mopend("Error in GETGEO");
        return;
    }
    if (!Angstroms) {
        for (int i = 1; i <= natoms; ++i) {
            geo[1][i] = geo[1][i] * a0;
            if (na[i] == 0) {
                geo[2][i] = geo[2][i] * a0;
                geo[3][i] = geo[3][i] * a0;
            }
        }
    }
    if (maxtxt > 0) {
        for (int i = 1; i <= natoms; ++i) {
            int j = std::max(1, len_trim(txtatm[i]));
            if (j != maxtxt) {
                txtatm[i].resize(std::max(j, maxtxt), ' ');
            }
        }
    }
    na[1] = 0; nb[1] = 0; nc[1] = 0;
    nb[2] = 0; nc[2] = 0;
    nc[3] = 0;
    if (velo) {
        int j = 0;
        for (int i = 1; i <= natoms; ++i) if (na[i] != 0) j = 1;
        if ((j != 0 || (int_ && keywrd.find(" LET") == std::string::npos)) &&
            keywrd.find(" 0SCF") == std::string::npos) {
            std::fprintf(stdout, " COORDINATES MUST BE CARTESIAN WHEN VELOCITY VECTOR IS USED.\n");
            mopend("COORDINATES MUST BE CARTESIAN WHEN VELOCITY VECTOR IS USED.");
            return;
        }
        react.assign(3 * natoms + 1, 0.0);
        for (int i = 1; i <= natoms; ++i) {
            std::getline(std::cin, line);
            double v0 = 0.0;
            nuchar(const_cast<char*>(line.data()), len_trim(line), &v0, ndmy);
            react[(i - 1) * 3 + 1] = v0;
            if (ndmy == 3) continue;
            std::fprintf(stdout, "  THERE MUST BE EXACTLY THREE VELOCITY DATA PER LINE\n");
            mopend("THERE MUST BE EXACTLY THREE VELOCITY DATA PER LINE.");
            return;
        }
    }
    if (ircdrc) {
        if (numat != natoms) {
            std::fprintf(stdout, " Only real atoms are allowed in IRC and DRC calculations.\n");
            mopend("Only real atoms are allowed in IRC and DRC calculations.");
            return;
        }
        // gmetry() consumes the module-level NA/NB/NC arrays; mirror our
        // parameters into them and copy any renumbering back afterwards.
        common_arrays_C::labels.assign(labels.begin(), labels.end());
        common_arrays_C::na.assign(na.begin(), na.end());
        common_arrays_C::nb.assign(nb.begin(), nb.end());
        common_arrays_C::nc.assign(nc.begin(), nc.end());
        gmetry(geo, xyz);
        na.assign(common_arrays_C::na.begin(), common_arrays_C::na.end());
        nb.assign(common_arrays_C::nb.begin(), common_arrays_C::nb.end());
        nc.assign(common_arrays_C::nc.begin(), common_arrays_C::nc.end());
        std::vector<double> flat(3 * (numat + 2), 0.0), out(3 * (numat + 2), 0.0);
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= 3; ++j) flat[(i - 1) * 3 + (j - 1)] = xyz[j - 1][i];
        xyzint(&flat[0], numat, &na[0], &nb[0], &nc[0], 1.0, &out[0]);
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= 3; ++j) geo[j][i] = out[(i - 1) * 3 + (j - 1)];
        for (int i = 1; i <= numat; ++i) {
            geo[1][i] = xyz[0][i];
            geo[2][i] = xyz[1][i];
            geo[3][i] = xyz[2][i];
        }
        na_store.resize(numat + 1);
        for (int i = 1; i <= numat; ++i) na_store[i] = na[i];
        std::fill(na.begin(), na.begin() + numat + 1, 0);
    }
    solid = false;
    for (int i = 1; i <= natoms; ++i) {
        if (labels[i] <= 0) {
            std::fprintf(stdout, " ATOMIC NUMBER OF %3d ?\n", labels[i]);
            if (i == 1) {
                std::fprintf(stdout, " THIS WAS THE FIRST ATOM\n");
            } else {
                std::fprintf(stdout, "    GEOMETRY UP TO, BUT NOT INCLUDING, THE FAULTY ATOM\n");
                natoms = i - 1;
                geout(iw);
            }
            mopend("Error in READMO");
            return;
        }
        if (labels[i] == 107) solid = true;
        if (solid && labels[i] < 99) {
            std::fprintf(stdout, " TRANSLATION VECTORS MUST BE AT THE END OF THE DATA SET\n");
            mopend("TRANSLATION VECTORS MUST BE AT THE END OF THE DATA SET");
            return;
        }
        if (i == 1 || na[i] == 0) continue;
        int j = 0;
        if ((na[i] == nb[i] && i > 1) ||
            ((na[i] == nc[i] || nb[i] == nc[i]) && i > 2) ||
            (nb[i] * nc[i] == 0 && i > 3)) j = 1;
        if (na[i] >= i || nb[i] >= i || nc[i] >= i) {
            j = j + 1;
            if ((na[i] < i || na[std::max(1, na[i])] == 0) &&
                (nb[i] < i || na[std::max(1, nb[i])] == 0) &&
                (nc[i] < i || na[std::max(1, nc[i])] == 0)) j = j - 1;
        }
        if (j == 0) continue;
        j = std::max(i, std::max(na[i], std::max(nb[i], nc[i])));
        char no[4];
        std::snprintf(no, sizeof(no), "%d", (int)std::log10(j + 1.1) + 1);
        char buf[96];
        std::snprintf(buf, sizeof(buf), " ATOM NUMBER %*d IS ILL-DEFINED", no[0] - '0', i);
        std::fprintf(stdout, "\n          %s\n", buf);
        std::fprintf(stdout, "\n          Connectivity of atom %*d: NA=%*d, NB=%*d, NC=%*d\n",
                     no[0] - '0', i, no[0] - '0', na[i], no[0] - '0', nb[i], no[0] - '0', nc[i]);
        if (na[i] > i && na[na[i]] != 0)
            std::fprintf(stdout, "\n          NA is defined using internal coordinates\n");
        if (nb[i] > i && na[nb[i]] != 0)
            std::fprintf(stdout, "\n          NB is defined using internal coordinates\n");
        if (nc[i] > i && na[nc[i]] != 0)
            std::fprintf(stdout, "\n          NC is defined using internal coordinates\n");
        web_message(0, "geometry_specification.html");
        if (i == 1) return;
        std::fprintf(stderr, "\n          %s\n", buf);
        std::fprintf(stdout, "\n  GEOMETRY READ IN\n\n");
        std::fill(nat.begin(), nat.end(), 0);
        geout(iw);
        mopend(buf);
        return;
    }
    if (natoms > 0) {
        common_arrays_C::labels.assign(labels.begin(), labels.end());
        common_arrays_C::na.assign(na.begin(), na.end());
        common_arrays_C::nb.assign(nb.begin(), nb.end());
        common_arrays_C::nc.assign(nc.begin(), nc.end());
        if ((int)xyz.size() < 3 || (int)xyz[0].size() < natoms + 1)
            xyz.assign(3, std::vector<double>((size_t)natoms + 1, 0.0));
        gmetry(geo, xyz);
        na.assign(common_arrays_C::na.begin(), common_arrays_C::na.end());
        nb.assign(common_arrays_C::nb.begin(), common_arrays_C::nb.end());
        nc.assign(common_arrays_C::nc.begin(), common_arrays_C::nc.end());
    } else {
        mopend("No atoms!");
    }
    if (moperr) return;
    if (lxyz || int_) {
        int k = 0;
        for (int i = 1; i <= natoms; ++i)
            for (int j = 1; j <= 3; ++j) k = k + lopt[j][i];
        numat = 0;
        for (int i = 1; i <= natoms; ++i) {
            if (labels[i] != 99) {
                numat = numat + 1;
                labels[numat] = labels[i];
                txtatm[numat] = txtatm[i];
                lopt[1][numat] = lopt[1][i];
                lopt[2][numat] = lopt[2][i];
                lopt[3][numat] = lopt[3][i];
                na[numat] = 0;
            }
            geo[1][i] = xyz[0][i];
            geo[2][i] = xyz[1][i];
            geo[3][i] = xyz[2][i];
        }
        if (k >= 3 * numat - 6) {
            int nm = std::min(3, numat);
            for (int i = 1; i <= nm; ++i) { lopt[1][i] = 1; lopt[2][i] = 1; lopt[3][i] = 1; }
        }
        natoms = numat;
        if (saddle)
            for (int i = 1; i <= numat; ++i) { lopt[1][i] = 1; lopt[2][i] = 1; lopt[3][i] = 1; }
    }
    if (int_) {
        if (id > 0) {
            for (int i = numat - id + 1; i <= numat; ++i) {
                xyz[0][i] = xyz[0][i] + xyz[0][1];
                xyz[1][i] = xyz[1][i] + xyz[1][1];
                xyz[2][i] = xyz[2][i] + xyz[2][1];
            }
        }
        std::vector<double> flat(3 * (numat + 2), 0.0), out(3 * (numat + 2), 0.0);
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= 3; ++j) flat[(i - 1) * 3 + (j - 1)] = xyz[j - 1][i];
        xyzint(&flat[0], numat, &na[0], &nb[0], &nc[0], 1.0, &out[0]);
        for (int i = 1; i <= numat; ++i)
            for (int j = 1; j <= 3; ++j) geo[j][i] = out[(i - 1) * 3 + (j - 1)];
        for (int i = 1; i <= 3; ++i)
            for (int j = i; j <= 3; ++j) lopt[j][i] = 0;
    }
    return;

label210:
    // Error conditions.
    {
        int j = natoms - 1;
        std::fprintf(stdout, " DATA CURRENTLY READ IN ARE: \n\n");
        for (int k = 1; k <= j; ++k) {
            if (na[k] > 0) {
                geo[2][k] = geo[2][k] / 1.7453292519943e-2;
                geo[3][k] = geo[3][k] / 1.7453292519943e-2;
            }
            std::fprintf(stdout, "%3s%10.5f%2d %10.5f%2d %10.5f%2d %2d %2d %2d\n",
                         kElem[labels[k]],
                         geo[1][k], lopt[1][k], geo[2][k], lopt[2][k], geo[3][k], lopt[3][k],
                         na[k], nb[k], nc[k]);
        }
        natoms = 0;
        return;
    }
}

void txt_to_atom_no(std::string& text, int j_in, bool let, int& m) {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace chanel_C;
    if (maxtxt < 26) {
        mopend("Labeled atoms can only be used when atoms have labels");
        m = 0;
        return;
    }
    int j = j_in + 1;
    size_t l = text.find('"', j - 1);
    l = (l == std::string::npos) ? text.size() : l + 1;
    std::string line = text.substr(j - 1, (l > (size_t)j) ? (l - j) : 0);
    std::string txt2 = trimr(line);
    if (!line.empty() && line[0] == '[') {
        int k = len_trim(line);
        if (k > 14) {
            char b[120];
            std::snprintf(b, sizeof(b), " CVB Atom defined by '%s' is not in JSmol format", trimr(line).c_str());
            mopend(b);
            return;
        }
        // Atom defined using Jmol format.
        size_t n = text.find('.', j - 1);
        n = (n == std::string::npos) ? 0 : n + 1 + j;
        size_t kk = text.find(':', j - 1);
        kk = (kk == std::string::npos) ? 0 : kk + 1 + j;
        size_t mm = text.find(']', j - 1);
        mm = (mm == std::string::npos) ? 0 : mm + 1 + j;
        if (n == j || kk == j || mm == j) {
            char b[120];
            std::snprintf(b, sizeof(b), " CVB Atom defined by '%s' is not in JSmol format", trimr(line).c_str());
            mopend(b);
            if (n == j) mopend("(The dot separator, \".\", is missing.)");
            if (kk == j) mopend("(The colon separator, \":\", is missing.)");
            if (mm == j) mopend("(The close square bracket, \"]\", is missing.)");
            return;
        }
        line = text.substr(n, l - 1 - n) + text.substr(j, mm - j - 1) + text.substr(kk, 1) +
               text.substr(mm, kk - 1 - mm);
    }
    m = 0;
    for (int k = 1; k <= len_trim(line); ++k) {
        if (line[k - 1] != ' ') {
            m = m + 1;
            line[m - 1] = line[k - 1];
        }
    }
    line.resize(m + 1, ' ');
    if (m > 12) {
        char b[120];
        std::snprintf(b, sizeof(b), " CVB Atom defined by '%s' is not in PDB format", trimr(line).c_str());
        mopend(b);
        return;
    }
    line.resize(m + 5, ' ');
    int matched = 0;
    for (int i = 1; i <= numat; ++i) {
        std::string txt = txtatm[i];
        int n = 0;
        for (int k = 13; k <= 26 && k <= (int)txt.size(); ++k) {
            if (txt[k - 1] != ' ') {
                n = n + 1;
                if (n <= (int)txt.size()) txt[n - 1] = txt[k - 1];
            }
        }
        txt.resize(std::max(n, m), '?');
        int n2 = std::max(n, m);
        int k = 1;
        for (; k <= n2; ++k) {
            if (txt[k - 1] != line[k - 1] && line[k - 1] != '*') break;
        }
        if (k > n2) {
            char b[16];
            std::snprintf(b, sizeof(b), "%d", i);
            std::string repl = b;
            text = text.substr(0, j - 1) + repl + text.substr(l - 1);
            matched = i;
            break;
        }
    }
    if (matched == 0) {
        char b[160];
        if (txt2.find('[') != std::string::npos) {
            std::snprintf(b, sizeof(b), " CVB Atom defined by JSmol label '%s' was not found in the data set", trimr(txt2).c_str());
        } else {
            std::snprintf(b, sizeof(b), " CVB Atom defined by PDB label '%s' was not found in the data set", trimr(txt2).c_str());
        }
        if (!let) {
            mopend(b);
            std::fprintf(stdout, "          (Hint: Check that the atom-label in the data set matches the atom-label in the PDB file.)\n");
            return;
        }
    }
    m = matched;
}
