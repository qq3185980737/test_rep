// gettxt.cpp — C++ translation of MOPAC 2016 "gettxt.F90".
// Reads the first three lines of a MOPAC input file: the keyword line
// (with optional continuation via '+' or '&'), the comment line and the
// title line. Also merges a SETUP file (keywords + parameters), removes
// duplicates (' -' handling), and converts PM3/PM5 keywords (GUI only).
#include "gettxt.h"
#include "getdat_lines_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace molkst_C {
extern std::string keywrd, koment, title, line, verson;
extern std::vector<std::string> refkey;
extern bool gui, in_house_only, moperr;
extern int numcal;
}
namespace chanel_C {
extern int ir, iw, isetup;
extern std::string input_fn;
}
extern void upcase(std::string&, int);
extern void mopend(const std::string&);
extern void web_message(int, const char*);
extern void add_path(std::string&);
extern std::string get_text(std::string&, int, int);

static std::string trimr(const std::string& s) {
    size_t b = s.find_last_not_of(' ');
    return b == std::string::npos ? "" : s.substr(0, b + 1);
}
static int len_trim(const std::string& s) {
    size_t b = s.find_last_not_of(' ');
    return (b == std::string::npos) ? 0 : (int)(b + 1);
}
// Fortran 1-based index(): returns position (1-based) or 0 when absent.
static int f_index(const std::string& s, const std::string& sub, int start1 = 1) {
    size_t p = s.find(sub, start1 - 1);
    return (p == std::string::npos) ? 0 : (int)p + 1;
}
// Fortran slice s(i:j), 1-based inclusive; j <= 0 means open-ended.
static std::string f_slice(const std::string& s, int i, int j) {
    if (i < 1) i = 1;
    if (j <= 0 || j > (int)s.size()) j = (int)s.size();
    if (i > j) return "";
    return s.substr(i - 1, j - i + 1);
}
// s(i:j) = repl (Fortran overwrite semantics; pad/truncate to slice width).
static void f_set(std::string& s, int i, int j, const std::string& repl) {
    if (i < 1) i = 1;
    if (j <= 0 || j > (int)s.size()) j = (int)s.size();
    if (i > j) return;
    std::string r = repl;
    if ((int)r.size() > j - i + 1) r.resize(j - i + 1);
    s.replace(i - 1, j - i + 1, r);
}
static bool read_stdin_line(std::string& l) {
    // Translated driver: records come from the in-memory queue filled by getdat(),
    // not from std::cin (Fortran reads unit "ir").
    if (getdat_line_idx >= getdat_lines.size()) return false;
    l = getdat_lines[getdat_line_idx++];
    if (!l.empty() && l.back() == '\r') l.pop_back();
    return true;
}

void gettxt() {
    using namespace molkst_C;
    using namespace chanel_C;
    koment = "    NULL  ";
    title = "    NULL  ";
    refkey.assign(6, "    NULL  ");
    bool aux = (keywrd.find("AUX") != std::string::npos);
    std::string oldkey;
    std::string filen, path;
    std::string l;
    bool zero_scf = false;

    // label60 body: final keyword handling.
    auto finish = [&]() {
        upcase(keywrd, len_trim(keywrd));
        if (gui) {
            int i = f_index(keywrd, "PM3");
            if (i != 0) {
                f_set(keywrd, i, i + 2, "PM6");
                std::fprintf(stdout, " Keyword PM3 was supplied. PM3 is not supported, so keyword converted to PM6\n");
            }
            i = f_index(keywrd, "PM5");
            if (i != 0) {
                f_set(keywrd, i, i + 2, "PM7");
                std::fprintf(stdout, " Keyword PM5 was supplied. PM5 is not supported, so keyword converted to PM7\n");
            }
        }
    };
    // label100 body: end-of-input handling.
    auto eof_handling = [&]() {
        if (numcal > 1) {
            if (f_index(keywrd, "OLDGEO") != 0) return;  // User forgot title/comment lines.
            if (aux) keywrd = " AUX";
            line = "JOB ENDED NORMALLY";
        } else if (f_index(keywrd, "GEO_DAT") == 0) {
            line = " ERROR IN READ OF FIRST THREE LINES";
        }
        if (trimr(line) != " ") mopend(trimr(line));
    };
    // label50 body: SETUP file missing/empty/corrupt.
    auto label50_body = [&]() {
        if (zero_scf) { finish(); return; }
        numcal = 2;
        std::fprintf(stdout, " SETUP FILE MISSING, EMPTY OR CORRUPT\n");
        std::fprintf(stdout, " (Setup file name: '%s')\n", trimr(filen).c_str());
        mopend("SETUP FILE MISSING, EMPTY OR CORRUPT");
    };

    if (!read_stdin_line(l)) { eof_handling(); return; }
    refkey[0] = l;
    keywrd = refkey[0];
    oldkey = keywrd;
    upcase(keywrd, len_trim(keywrd));
    zero_scf = (keywrd.find("0SCF") != std::string::npos);
    int ipath = 0;
    for (int i = len_trim(input_fn); i >= 2; --i) {
        if (input_fn[i - 1] == '\\' || input_fn[i - 1] == '/') { ipath = i; break; }
    }
    path = (ipath > 2) ? input_fn.substr(0, ipath) : "";
    filen = (ipath > 2) ? path + "SETUP" : "SETUP";
    if (in_house_only) {
        if (f_index(keywrd, "PM6") == 0 && verson.size() >= 7 && verson[6] == 'M')
            filen = "/Users/jstewart/SETUP.txt";
    }
    std::ifstream f(filen.c_str());
    bool exists = f.good();
    f.close();
    int i = len_trim(keywrd);
    bool setup_present = (f_index(keywrd, "SETUP") != 0);
    if (setup_present || (in_house_only && exists)) {
        i = f_index(keywrd, "SETUP=");
        if (i != 0) {
            filen = get_text(oldkey, i + 6, 1);
        } else {
            if (in_house_only && verson.size() >= 7 && verson[6] == 'M') {
                filen = "SETUP";
            } else {
                filen = (ipath > 2) ? path + "SETUP" : "SETUP";
            }
        }
        add_path(filen);
        std::ifstream f2(filen.c_str());
        exists = f2.good();
        f2.close();
        if (!exists) {
            std::ifstream f3((trimr(filen) + ".txt").c_str());
            exists = f3.good();
            f3.close();
            if (exists) filen = trimr(filen) + ".txt";
        }
        if (!exists) {
            if (setup_present && !zero_scf) {
                std::fprintf(stdout, " SETUP FILE MISSING\n");
                std::fprintf(stdout, " (Setup file name: '%s')\n", trimr(filen).c_str());
                numcal = 2;
                if (!gui) std::fprintf(stderr, "\n%30s\n", "SETUP FILE MISSING, EMPTY OR CORRUPT");
                if (!gui) std::fprintf(stderr, "\n%2s\n\n", ("(An attempt was made to open the Setup file named: '" + trimr(filen) + "')").c_str());
                mopend("SETUP FILE MISSING");
                return;
            }
        } else {
            std::ifstream fs(filen.c_str(), std::ios::in);
            if (!fs) {
                if (!zero_scf) {
                    mopend("COULD NOT OPEN SETUP FILE: " + trimr(filen));
                    if (zero_scf) moperr = false;
                    return;
                }
            } else {
                if (!std::getline(fs, l)) { fs.close(); label50_body(); return; }
                if (!l.empty() && l.back() == '\r') l.pop_back();
                refkey[1] = l;
                fs.close();
                upcase(refkey[1], len_trim(refkey[1]));
                // Check for " -" signs in the SETUP file.
                if (refkey[1].empty() || refkey[1][0] != ' ') refkey[1] = " " + trimr(refkey[1]);
                for (;;) {
                    int ii = f_index(refkey[1], " -");
                    if (ii == 0) break;
                    int jj = f_index(f_slice(refkey[1], ii + 2, 0), " ") + ii + 1;
                    for (;;) {
                        int kk = f_index(" " + keywrd, " " + f_slice(refkey[1], ii + 2, jj - 1));
                        if (kk == 0) break;
                        int ll = f_index(f_slice(keywrd, kk + 1, 0), " ") + kk + 1;
                        keywrd = keywrd.substr(0, kk - 1) + keywrd.substr(ll - 1);
                    }
                    refkey[1] = refkey[1].substr(0, ii - 1) + refkey[1].substr(jj - 1);
                }
                // Keywords on the keyword line take precedence: delete any in SETUP.
                i = 1;
                for (;;) {
                    i = i + 1;
                    if (i > len_trim(refkey[1])) break;
                    if (refkey[1][i - 2] == ' ' && refkey[1][i - 1] != ' ') {
                        int jj = i + 1;
                        for (; jj <= len_trim(refkey[1]); ++jj)
                            if (refkey[1][jj - 1] == ' ') break;
                        std::string line_ = f_slice(refkey[1], i, jj);
                        int k = std::min(6, jj - i);
                        if (f_index(keywrd, " " + line_.substr(0, k)) > 0) {
                            refkey[1] = refkey[1].substr(0, i - 1) + refkey[1].substr(jj);
                            i = i - 1;
                        }
                    }
                }
                i = len_trim(keywrd);
                keywrd = keywrd.substr(0, i) + " " + refkey[1].substr(0, std::max(0, 999 - i));
                // Check for " -" signs in the keyword line.
                for (;;) {
                    int ii = f_index(keywrd, " -");
                    if (ii == 0) break;
                    int jj = f_index(f_slice(keywrd, ii + 2, 0), " ") + ii + 1;
                    for (;;) {
                        int kk = f_index(" " + keywrd, " " + f_slice(keywrd, ii + 2, jj - 1));
                        if (kk == 0) break;
                        int ll = f_index(f_slice(keywrd, kk + 1, 0), " ") + kk + 1;
                        keywrd = keywrd.substr(0, kk - 1) + keywrd.substr(ll - 1);
                    }
                    ii = f_index(keywrd, " -");
                    jj = f_index(f_slice(keywrd, ii + 2, 0), " ") + ii + 1;
                    keywrd = keywrd.substr(0, ii - 1) + keywrd.substr(jj - 1);
                }
                // Remove excess EXTERNALS (in-house builds only).
                if (in_house_only) {
                    for (;;) {
                        int ii = f_index(keywrd, " EXTER");
                        if (ii == 0) break;
                        int jj = f_index(f_slice(keywrd, ii + 3, 0), " EXTER");
                        if (jj != 0) {
                            jj = f_index(f_slice(keywrd, ii + 1, 0), " ") + ii;
                            keywrd = keywrd.substr(0, ii - 1) + " " + keywrd.substr(jj - 1);
                        } else {
                            break;
                        }
                    }
                }
                refkey[0] = trimr(keywrd);
                refkey[1] = refkey[2];
                if (trimr(keywrd.substr(i)) == "") { label50_body(); return; }
            }
        }
        if (!read_stdin_line(l)) { eof_handling(); return; }
        koment = l;
        if (!read_stdin_line(l)) { eof_handling(); return; }
        title = l;
        finish();
        return;
    } else if (f_index(keywrd.substr(0, i), " +") != 0) {
        // READ SECOND KEYWORD LINE.
        int ii = f_index(keywrd.substr(0, i), " +");
        f_set(keywrd, ii, ii + 1, "  ");
        i = len_trim(keywrd);
        if (!read_stdin_line(l)) { eof_handling(); return; }
        refkey[1] = l;
        keywrd = keywrd.substr(0, i) + refkey[1].substr(0, std::max(0, 999 - i));
        oldkey = keywrd;
        upcase(keywrd, len_trim(keywrd));
        if (f_index(keywrd, "SETUP") != 0) {
            int j;
            ii = f_index(keywrd, "SETUP=");
            if (ii != 0) {
                j = f_index(keywrd.substr(ii - 1), " ");
                filen = oldkey.substr(ii + 5, j - 2);
                f_set(keywrd, ii, ii + j, "  ");
            } else {
                filen = "SETUP";
                ii = f_index(keywrd, "SETUP");
                f_set(keywrd, ii, ii + 5, "      ");
            }
            f_set(keywrd, ii, ii + 6, "       ");
            add_path(filen);
            std::ifstream fs(filen.c_str(), std::ios::in);
            if (fs) {
                if (std::getline(fs, l)) {
                    if (!l.empty() && l.back() == '\r') l.pop_back();
                    refkey[1] = l;
                }
                fs.close();
                i = len_trim(keywrd) + 1;
                keywrd = keywrd.substr(0, i - 1) + refkey[1].substr(0, std::max(0, 1001 - i));
                upcase(keywrd, len_trim(keywrd));
            }
            // label30: continue below
        } else if (f_index(keywrd.substr(i), " +") != 0) {
            // READ THIRD KEYWORD LINE.
            if (!read_stdin_line(l)) { eof_handling(); return; }
            refkey[2] = l;
            int ii = f_index(refkey[2], " + ");
            if (ii != 0) {
                std::fprintf(stdout, " A maximum of three lines of keywords are allowed.\n");
                std::fprintf(stdout, " On the third line of keywords is a '+' sign, implying more lines of keywords.\n");
                std::fprintf(stdout, " Remove the '+' sign from the third line of keywords, and re-run.\n");
                web_message(iw, "plus.html");
                mopend("A maximum of three lines of keywords are allowed.");
                return;
            }
            ii = f_index(keywrd.substr(0, len_trim(keywrd)), " +");
            f_set(keywrd, ii, ii + 1, "  ");
            i = len_trim(keywrd);
            keywrd = keywrd.substr(0, i) + refkey[2].substr(0, std::max(0, 1001 - i));
            upcase(keywrd, len_trim(keywrd));
        }
        // READ TITLE LINE.
        if (!read_stdin_line(l)) { eof_handling(); return; }
        koment = l;
        if (!read_stdin_line(l)) { eof_handling(); return; }
        title = l;
        finish();
        return;
    } else if (f_index(keywrd, "&") != 0) {
        int ii = f_index(keywrd, "&");
        f_set(keywrd, ii, ii, " ");
        i = len_trim(keywrd);
        if (!read_stdin_line(l)) { eof_handling(); return; }
        refkey[1] = l;
        keywrd = keywrd.substr(0, i) + " " + refkey[1].substr(0, std::max(0, 1001 - i));
        oldkey = keywrd;
        upcase(keywrd, len_trim(keywrd));
        if (f_index(keywrd, "SETUP") != 0) {
            int j;
            ii = f_index(keywrd, "SETUP=");
            if (ii != 0) {
                j = f_index(keywrd.substr(ii - 1), " ");
                filen = oldkey.substr(ii + 5, j);
                f_set(keywrd, ii, ii + j, "  ");
            } else {
                filen = "SETUP";
                ii = f_index(keywrd, "SETUP");
                f_set(keywrd, ii, ii + 6, "       ");
            }
            add_path(filen);
            std::ifstream fs(filen.c_str(), std::ios::in);
            if (fs) {
                if (std::getline(fs, l)) {
                    if (!l.empty() && l.back() == '\r') l.pop_back();
                    keywrd = keywrd.substr(0, len_trim(keywrd) + 1) + l;
                }
                fs.close();
                upcase(keywrd, len_trim(keywrd));
                if (!read_stdin_line(l)) { eof_handling(); return; }
                title = l;
            }
            // label40: continue below
        } else if (f_index(keywrd.substr(ii - 1), "&") != 0) {
            if (!read_stdin_line(l)) { eof_handling(); return; }
            keywrd = keywrd.substr(0, len_trim(keywrd) - 1) + l;
            upcase(keywrd, len_trim(keywrd));
        } else {
            if (!read_stdin_line(l)) { eof_handling(); return; }
            title = l;
        }
        finish();
        return;
    } else {
        if (!read_stdin_line(l)) { eof_handling(); return; }
        koment = l;
        if (!read_stdin_line(l)) { eof_handling(); return; }
        title = l;
        finish();
        return;
    }
}

// get_text: return text between character i_start and the next space.
// If character i_start is '"' or '\'', return text between i_start+1 and
// the closing quote. When zero==0 the consumed span is blanked out.
std::string get_text(std::string& line, int i_start, int zero) {
    // Fortran 1-based.
    const char limits[2] = { '"', '\'' };
    int i = 1;
    for (; i <= 2; ++i)
        if (i_start <= (int)line.size() && line[i_start - 1] == limits[i - 1]) break;
    int j = i_start;
    char ch;
    if (i > 2) {
        ch = ' ';
        i = i_start;
    } else {
        j = j + 1;
        ch = limits[i - 1];
        i = i_start + 1;
    }
    while (i + 1 <= (int)line.size() && line[i] != ch) ++i;
    std::string out = (j <= (int)line.size()) ? line.substr(j - 1, i - j + 1) : "";
    if (zero == 0) {
        int upto = std::min(i + 1, (int)line.size());
        for (int k = i_start; k <= upto; ++k) line[k - 1] = ' ';
    }
    return out;
}
