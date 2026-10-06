// getdat.cpp — C++ translation of MOPAC 2016 "getdat.F90".
// Reads the data-set (command-line argument), resolves <file>.mop/.dat/.arc,
// copies cleaned lines to a scratch buffer, and returns the line count.
#include "getdat_lines_C.h"
#include "getdat.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#if defined(_MSC_VER)
extern "C" { int __argc; char** __argv; }
#endif

namespace molkst_C {
extern int natoms, ncomments, ijulian, run;
extern bool gui, is_PARAM;
extern std::string keywrd, line, jobnam, verson;
extern double arc_hof_1, arc_hof_2;
}
namespace chanel_C {
extern int iw0, iw;
extern std::string job_fn, input_fn;
}
namespace common_arrays_C {
extern std::vector<std::string> all_comments;
}

extern void upcase(std::string&, int);
extern void init_filenames();
extern void mopend(const std::string&);
extern void to_screen(const std::string&);
extern double reada(const std::string&, int);
extern void web_message(int, const char*);

namespace {
std::string trimmed(const std::string& s) {
    size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}
}  // namespace

// The cleaned data-set lines are placed here (analogous to the Fortran
// scratch unit "input"); a consumer (getgeo/getgeg) reads from this.
std::vector<std::string> getdat_lines;
std::size_t getdat_line_idx = 0;  // defined here; shared via getdat_lines_C.h

void getdat(int input, int output) {
    using namespace molkst_C;
    using namespace chanel_C;
    using namespace common_arrays_C;
    getdat_lines.clear();
    natoms = 0;
    to_screen("To_file: getdat");
    bool exists = false, arc_file = false;
    if (gui) {
        jobnam = "MOPAC input";
        natoms = 1;
    } else {
        if (run != 2 || jobnam == " ") {
            int i = __argc - 1;
            if (i >= run) {
                jobnam = __argv[run];
                natoms = 1;
                // Strip trailing unprintable characters.
                size_t e = jobnam.size();
                while (e > 0) {
                    unsigned char c = (unsigned char)jobnam[e - 1];
                    if (c > 39 && c < 126) break;
                    --e;
                }
                jobnam = jobnam.substr(0, e);
            } else if (i == 0) {
                if (is_PARAM) {
                    line = " PARAM is the parameter optimization program for use with MOPAC";
                    std::fprintf(stderr, "\n          %s\n", trimmed(line).c_str());
                    mopend(trimmed(line));
                    std::fprintf(stderr, "          It uses a single argument, the PARAM data-set\n");
                    std::fprintf(stderr, "          The command to run PARAM is 'PARAM.exe <data-set>.dat'\n");
                    std::fprintf(stderr, "          Press (return) to continue\n");
                    std::string dummy;
                    std::getline(std::cin, dummy);
                    return;
                } else {
                    line = " MOPAC is a semiempirical quantum chemistry program";
                    std::fprintf(stderr, "\n          %s\n", trimmed(line).c_str());
                    std::fprintf(stderr, "          It uses a single argument, the MOPAC data-set\n");
                    std::fprintf(stderr, "          The command to run MOPAC is 'MOPAC2016.exe <data-set>.mop'\n");
                    web_message(0, "running_MOPAC.html");
                    std::fprintf(stderr, "          Press (return) to continue\n");
                    std::string dummy;
                    std::getline(std::cin, dummy);
                    return;
                }
            }
        } else {
            natoms = 1;
        }
    }
    if (natoms == 0) return;
    line = jobnam;
    std::string up = line;
    upcase(up, (int)trimmed(up).size());
    int idx = (int)(up.find(".MOP") != std::string::npos) +
              (int)(up.find(".DAT") != std::string::npos) +
              (int)(up.find(".ARC") != std::string::npos);
    arc_file = (up.find(".ARC") != std::string::npos);
    if (idx > 0) {
        // User supplied a suffix — use it directly.
        line = jobnam;
        std::ifstream f(line.c_str());
        exists = f.good();
        f.close();
        size_t m = jobnam.rfind('.');
        if (m != std::string::npos) jobnam = jobnam.substr(0, m);
    } else {
        line = jobnam;
        std::ifstream f(line.c_str());
        exists = f.good();
        f.close();
        if (!exists) {
            line = jobnam + ".mop";
            std::ifstream f2(line.c_str());
            exists = f2.good();
            f2.close();
        }
        if (!exists) {
            line = jobnam + ".dat";
            std::ifstream f3(line.c_str());
            exists = f3.good();
            f3.close();
        }
    }
    if (exists) {
        if (iw0 > -1) {
            if (ijulian < 365) {
                to_screen("********************************************************************************");
                char txt[128];
                std::snprintf(txt, sizeof(txt), "**                          MOPAC (%s)        Days remaining%4d         **", verson.c_str(), ijulian);
                to_screen(txt);
                to_screen("********************************************************************************");
            }
            to_screen("Preparing to read the following MOPAC file: ");
            int i = (int)trimmed(line).size();
            i = i < 240 ? i : 240;
            to_screen(line.substr(0, i < 120 ? i : 120));
            if (i > 120) to_screen(line.substr(120, i - 120 > 120 ? 120 : i - 120));
        }
        job_fn = line;
        std::ifstream src(job_fn.c_str());
        if (!src.is_open()) {
            line = " Data file: '" + trimmed(job_fn) + "' exists, but it cannot be opened.";
            std::fprintf(stderr, "\n          %s\n\n", trimmed(line).c_str());
            to_screen(" File '" + job_fn + "' cannot be opened");
            mopend("File '" + job_fn + "' cannot be opened");
            return;
        }
        init_filenames();
        arc_hof_1 = 0.0;
        arc_hof_2 = 0.0;
        if (arc_file) {
            // Scan for a geometry section header in the archive file.
            std::string l2;
            bool keep_arc = false;
            while (std::getline(src, l2)) {
                if (l2.find("HEAT OF FORMATION") != std::string::npos) arc_hof_1 = reada(l2, 20);
                if (l2.find("FINAL GEOMETRY OBTAINED") != std::string::npos ||
                    l2.find("GEOMETRY IN CARTESIAN COORDINATE") != std::string::npos ||
                    l2.find("GEOMETRY IN MOPAC Z-MATRIX") != std::string::npos) {
                    keep_arc = true;
                    break;
                }
            }
            if (!keep_arc) {
                std::fprintf(stderr, "The data set was defined as an ARC file\n");
                std::fprintf(stderr, "but it does not appear to an archive file\n");
                std::fprintf(stderr, "An attempt will be made to read it as a normal data set\n");
                src.clear();
                src.seekg(0);
            }
        }
        int nlines = 0;
        ncomments = 0;
        std::vector<std::string> tmp_comments;
        keywrd = " ";
        std::string data_line, line1;
        bool first_loop = true, stop = false;
        while (std::getline(src, data_line)) {
            if (!data_line.empty() && data_line.back() == '\r') data_line.pop_back();
            nlines++;
            if (first_loop) {
                first_loop = false;
                line1 = data_line;
                std::string up1 = line1;
                upcase(up1, (int)trimmed(up1).size());
                if (up1.find("DATA=") != std::string::npos) break;
            }
            if (data_line.empty() || data_line[0] != '*') {
                // Not a comment: expand tabs to spaces.
                std::string out;
                int col = 0;
                for (char ch : data_line) {
                    if (ch == '\t') {
                        out += ' ';
                        col++;
                        int l = col % 8;
                        if (l != 0) { for (int q = 0; q < 8 - l; ++q) { out += ' '; col++; } }
                    } else {
                        out += ch;
                        col++;
                    }
                }
                getdat_lines.push_back(trimmed(out));
                if (keywrd == " ") keywrd = trimmed(out);
            } else {
                ncomments++;
                tmp_comments.push_back(trimmed(data_line));
            }
        }
        getdat_line_idx = 0;   // rewind record queue for gettxt/getgeo consumers
        // Uppercase the keyword line.
        std::string k = keywrd;
        upcase(k, (int)k.size());
        keywrd = k;
        if (!keywrd.empty() && keywrd[0] != ' ') keywrd = " " + trimmed(keywrd);
        if (keywrd.find(" GEO_DAT") != std::string::npos) {
            nlines += 3;
        }
        keywrd = " ";
        all_comments.resize(ncomments + 101);
        for (int i = 1; i <= ncomments; ++i) all_comments[i] = tmp_comments[i - 1];
        natoms = nlines;
        src.close();
        return;
    }
    // Input file missing.
    {
        std::string jn = jobnam;
        while (!jn.empty() && jn.back() == ' ') jn.pop_back();
        size_t e = jn.size();
        while (e > 0 && (jn[e - 1] >= '0' && jn[e - 1] <= '9')) --e;
        std::string prefix = jn.substr(0, e);
        bool pwd = !prefix.empty() && (prefix.back() == 'a' || prefix.back() == 'A') && jn.size() - e == 8;
        if (pwd) {
            line = " Password detected, but password has already been correctly installed.";
            mopend(trimmed(line));
            line = " MOPAC is ready to run data sets.";
            mopend(trimmed(line));
            return;
        }
        line = " The input data file \"" + trimmed(jobnam) + "\" does not exist.";
        std::string out_fn = trimmed(jobnam) + ".out";
        std::ofstream out(out_fn.c_str());
        if (!out.is_open()) {
            line = "Cannot write to output file '" + trimmed(jobnam) + ".out'. Reason: \"Permission to access file denied\"";
            std::fprintf(stderr, "\n          %s\n", trimmed(line).c_str());
            return;
        }
        std::fprintf(stderr, "\n          %s\n", trimmed(line).c_str());
        to_screen(trimmed(line));
        mopend(trimmed(line));
        return;
    }
}
