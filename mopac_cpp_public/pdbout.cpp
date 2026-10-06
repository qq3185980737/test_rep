// pdbout.cpp — C++ translation of MOPAC 2016 "pdbout.F90".
// Writes the current geometry in Brookhaven Protein Data Bank format.
// Includes the HTML/JSmol viewer page (write_html) and the data table
// (write_data_to_html) triggered by the HTML keyword.
#include "pdbout.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "elemts_C.h"
#include "MOZYME_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace molkst_C {
extern int numat, natoms, ncomments, nbreaks, maxtxt, numcal, nelecs;
extern std::string verson, line, keywrd, koment, title, formula, geo_ref_name, geo_dat_name;
extern double escf, arc_hof_1, arc_hof_2;
}
namespace chanel_C {
extern int iw;
extern std::string input_fn;
extern bool log;
}
namespace common_arrays_C {
extern std::vector<std::string> txtatm, txtatm1, all_comments;
extern std::vector<std::vector<double>> coord;
extern std::vector<int> nat, breaks;
extern std::vector<double> p;
}
namespace parameters_C {
extern double tore[];
}
namespace elemts_C {
extern std::vector<std::string> elemnt;
}
namespace MOZYME_C {
extern std::vector<std::string> tyres, tyr;
}

extern void chrge(const std::vector<double>&, std::vector<double>&);
extern void empiri();
extern const char* elemnt(int);
extern double reada(const std::string&, int);
void write_html();
void write_data_to_html(int iprt);
void add_path(std::string& line);

static std::string fdate() {
    std::time_t t = std::time(nullptr);
    std::tm tm;
    localtime_s(&tm, &t);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%a %b %d %H:%M:%S %Y", &tm);
    return std::string(buf);
}
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(' ');
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(' ');
    return s.substr(a, b - a + 1);
}
static int len_trim(const std::string& s) {
    size_t b = s.find_last_not_of(' ');
    return (b == std::string::npos) ? 0 : (int)(b + 1);
}

void pdbout(int mode1) {
    using namespace molkst_C;
    using namespace chanel_C;
    using namespace common_arrays_C;
    using namespace parameters_C;
    bool html = (keywrd.find(" HTML") != std::string::npos);
    std::vector<double> q2(numat + 2, 0.0);
    if (html && !p.empty()) {
        std::vector<double> cq(numat + 2, 0.0);
        chrge(p, cq);
        for (int i = 1; i <= numat; ++i) q2[i] = tore[nat[i]] - cq[i];
    } else {
        for (int i = 1; i <= numat; ++i) q2[i] = 0.0;
    }
    int nline = 0;
    int iprt = (mode1 == 1) ? iw : std::abs(mode1);
    std::string idate = fdate();
    int i = len_trim(input_fn);
    auto header = [&]() {
        line = "HEADER  data-set: " + input_fn.substr(0, (i > 5) ? i - 5 : 0);
        if (len_trim(line) > 80) {
            for (;;) {
                size_t k = line.find('\\');
                if (k == std::string::npos) break;
                size_t ii = line.find("data-set:") + 9;
                line = line.substr(0, ii) + line.substr(k + 1);
            }
        }
        line.resize(80, ' ');
        std::fprintf(stdout, "%s\n", trim(line).c_str());
        line = "REMARK  MOPAC 2012, Version: " + verson + " Date: " +
               idate.substr(4, 7) + idate.substr(20) + idate.substr(10, 6);
        std::fprintf(stdout, "%s\n", trim(line).c_str());
    };
    if (ncomments > 0) {
        if (all_comments[1].find("HEADER") == std::string::npos) {
            header();
        }
        for (i = 1; i <= ncomments; ++i) {
            std::string l7 = all_comments[i].substr(0, 7);
            const char* tags[] = { "ATOM  ", "HETATM", "TITLE ", "HEADER", "ANISOU", "COMPND",
                "SOURCE", "KEYWDS", "HELIX ", "SHEET ", "REMARK", "USER  ", "EXPDTA", "AUTHOR",
                "REVDAT", "JRNL  ", "DBREF ", "SEQRES", "HET   ", "HETNAM", "LINK  ", "CRYST1",
                "SCALE", "ORIGX", "FORMUL", "SEQRES" };
            bool hit = false;
            for (const char* t : tags)
                if (l7.find(t) != std::string::npos) { hit = true; break; }
            if (!hit) continue;
            std::fprintf(stdout, "%s\n", all_comments[i].substr(1, len_trim(all_comments[i]) - 1).c_str());
        }
    } else {
        header();
    }
    if (maxtxt == 0 && txtatm[1] != " ") {
        for (i = 1; i <= natoms; ++i) txtatm1[i] = txtatm[i];
    }
    int ii = 0, iii = 0;
    nbreaks = 1;
    for (i = 1; i <= numat; ++i) {
        ii = ii + 1;
        iii = iii + 1;
        if ((int)txtatm[iii].size() >= 14 && txtatm[iii][13] == 'X') iii = iii + 1;
        nline = nline + 1;
        bool ter = (i == breaks[nbreaks]);
        if (ter) nbreaks = nbreaks + 1;
        std::string en = elemnt(nat[i]);
        std::string ele_pdb(2, ' ');
        if (en[0] == ' ') {
            ele_pdb[0] = ' ';
            ele_pdb[1] = en[1];
        } else {
            ele_pdb[0] = en[0];
            char c2 = en[1];
            if (c2 >= 'a' && c2 <= 'z') c2 = (char)(c2 - 'a' + 'A');
            ele_pdb[1] = c2;
        }
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%s%5d%s%12.3f%8.3f%8.3f%s%5.2f%s%2s%s",
                      txtatm[iii].substr(0, 6).c_str(), ii,
                      txtatm[iii].substr(11, 15).c_str(),
                      coord[0][i], coord[1][i], coord[2][i],
                      "  1.00 ", q2[i], "      PROT", ele_pdb.c_str(), " ");
        std::fprintf(stdout, "%s\n", buf);
        if (ter) {
            ii = ii + 1;
            std::fprintf(stdout, "TER%8d      %s\n", ii,
                         txtatm[iii].substr(17, 9).c_str());
        }
    }
    std::fprintf(stdout, "END\n");
    if (html) write_html();
}

void write_html() {
    using namespace molkst_C;
    using namespace chanel_C;
    using namespace common_arrays_C;
    using namespace MOZYME_C;
    int nres = 0, i, j, k, nprt, ncol, biggest_res, iprt = 27, it;
    static int icalcn = -1;
    if (icalcn == numcal) return;
    icalcn = numcal;
    size_t i0 = keywrd.find(" HTML");
    size_t j0 = keywrd.find(' ', i0 + 1);
    bool l_prt_res = (keywrd.substr(i0 + 1, j0 - i0 - 1).find("NORES") == std::string::npos);
    line = input_fn.substr(0, len_trim(input_fn) - 4) + "html";
    std::ofstream fout(line, std::ios::out | std::ios::trunc);
    if (!fout) return;
    std::vector<std::string> res_txt(4001, "          ");
    std::string l_res(10, ' '), n_res(11, ' '), wrt_res(10, ' '), num(1, '1'), line_1(400, ' ');
    const int limres = 260;
    if (l_prt_res) {
        for (i = 1; i <= numat; ++i) {
            if (txtatm[i].substr(17, 3) == "HOH") continue;
            if (txtatm[i].substr(17, 3) == "SO4") continue;
            j = (int)(reada(txtatm[i], 23) + 0.5);
            l_res = std::string(10, ' ');
            char b[8];
            std::snprintf(b, sizeof(b), "%d", j);
            std::string ns = txtatm[i].substr(17, 3) + b + ":" + txtatm[i].substr(21, 1);
            l_res = ns + std::string(10 - (int)ns.size(), ' ');
            for (k = 1; k <= 6; ++k)
                if (l_res[k - 1] != ' ') break;
            if (l_res[k - 1] >= '0' && l_res[k - 1] <= '9') l_res[k - 1] = 'Q';
            for (k = 1; k <= nres; ++k)
                if (res_txt[k] == l_res) break;
            if (k > nres) {
                nres = nres + 1;
                res_txt[nres] = l_res;
            }
        }
        for (i = 1; i <= nres; ++i) {
            for (j = 1; j <= 20; ++j)
                if (res_txt[i].substr(0, 3) == tyres[j]) break;
            if (j < 21) break;
        }
        num = txtatm[1].substr(21, 1);
        for (i = 1; i <= nres; ++i)
            if (res_txt[i].substr(8, 1) != num) break;
        if (i > nres) {
            for (i = 1; i <= nres; ++i) {
                if (res_txt[i].size() >= 9) {
                    res_txt[i][7] = ' ';
                    res_txt[i][8] = ' ';
                }
            }
        }
    }
    // Heading.
    if (len_trim(koment) == 0 || koment.substr(0, 8).find(" NULL") != std::string::npos) {
        fout << "<HTML><HEAD><TITLE>" << input_fn.substr(0, len_trim(input_fn) - 5) << "</TITLE></HEAD>\n";
    } else {
        for (it = 1; it <= len_trim(koment); ++it)
            if (koment[it - 1] != ' ') break;
        fout << "<HTML><HEAD><TITLE>" << trim(koment.substr(it - 1)) << "</TITLE></HEAD>\n";
    }
    fout << "<style type=\"text/css\">\n";
    fout << ".auto-style4 {\n";
    fout << "margin-top: 0px;\n";
    fout << "text-align: center;\n";
    fout << "line-height: 85%;\n";
    fout << "}\n";
    fout << ".auto-style5 {\n";
    fout << "text-decoration: none;\n";
    fout << "}\n";
    fout << "</style>\n";
    fout << "<!--   Start of JSmol script    -->\n";
    fout << "<meta charset=\"utf-8\"> <script type=\"text/javascript\" src=\"../jsmol/JSmol.min.js\"></script>  \n";
    fout << "<script type=\"text/javascript\">\n";
    fout << "\n$(document).ready(function() {Info = {\n";
    fout << "          width: 1500,\n";
    fout << "          height: 1000,\n";
    fout << "          color: \"0xB0B0B0\",\n";
    fout << "          disableInitialConsole: true, \n";
    fout << "          addSelectionOptions: false,\n";
    fout << "          j2sPath: \"../jsmol/j2s\",\n";
    fout << "          jarPath: \"../jsmol/java\",\n";
    fout << "          use: \"HTML5\", script:  \n";
    fout << "// Data set to be loaded\n";
    line = input_fn.substr(0, len_trim(input_fn) - 4) + "pdb";
    for (i = len_trim(line); i >= 1; --i)
        if (line[i - 1] == '/' || line[i - 1] == '\\') break;
    fout << " \"load \\\\\n";
    bool l_geo_ref = (keywrd.find(" 0SCF") != std::string::npos && keywrd.find(" GEO_REF") != std::string::npos);
    if (l_geo_ref) {
        line = geo_dat_name.substr(0, len_trim(geo_dat_name) - 3) + "pdb" + "' '" +
               geo_ref_name.substr(0, len_trim(geo_ref_name) - 3) + "pdb" + "'; \\\\\n";
        fout << "          FILES '" << trim(line) << "\n";
    } else {
        fout << "          '" << line.substr(i) << "'; \\\\\n";
    }
    fout << "          set measurementUnits ANGSTROMS; \\\\\n";
    if (l_geo_ref) {
        fout << "          set bondRadiusMilliAngstroms (25); \\\\\n          spacefill 10%; \\\\\n";
    } else {
        fout << "          set bondRadiusMilliAngstroms (50); \\\\\n          spacefill 15%; \\\\\n";
    }
    fout << "          set display selected; \\\\\n          hBonds calculate; \\\\\n";
    fout << "          set defaultDistanceLabel '%0.3VALUE %UNITS'; \\\\\n";
    if (keywrd.find(" 0SCF") != std::string::npos && keywrd.find(" GEO_REF") != std::string::npos) {
        fout << "          select */2.1; color bonds green;  select off;\\\n";
    } else {
        fout << "          select off; \\\\\n";
    }
    fout << "          set perspectivedepth off; \\\\\n";
    fout << "          connect 0.8  1.5 (hydrogen) (phosphorus) create; \\\\\n";
    if (l_geo_ref) {
        fout << "          set zoomLarge false; frame 0;\" \n";
    } else {
        fout << "          set zoomLarge false;\" \n";
    }
    fout << "} \n";
    fout << "$(\"#mydiv\").html(Jmol.getAppletHtml(\"jmolApplet0\",Info))}); \n";
    fout << "</script>\n";
    fout << "<!--   End of JSmol script    -->\n";
    if (len_trim(koment) == 0 || koment.substr(0, 8).find(" NULL") != std::string::npos) {
        fout << "<h1 align=\"center\">" << input_fn.substr(0, len_trim(input_fn) - 5) << "</h1>\n";
    } else {
        fout << "<h1 align=\"center\">" << trim(koment.substr(it - 1)) << "</h1>\n";
    }
    if (len_trim(title) != 0 && title.substr(0, 8).find(" NULL") == std::string::npos)
        fout << "<h2 align=\"center\">" << trim(title) << "</h2>\n";
    if (keywrd.find(" 0SCF") != std::string::npos && keywrd.find(" GEO_REF") != std::string::npos) {
        fout << "<h2 align=\"center\">Compare \"" << trim(geo_dat_name)
             << "\" and \"<span style=\"color:green\">" << trim(geo_ref_name) << "</span>\"</h2>\n";
    }
    fout << "<BODY  BGCOLOR=\"#ffffff\">\n";
    fout << "<TABLE>\n<TD>\n<TABLE>\n<TR>\n";
    if (nres < limres) {
        fout << "<TD colspan=\"2\">\n";
        write_data_to_html(iprt);
        fout << "</TD></TR><TR>\n";
    }
    if (keywrd.find(" 0SCF") != std::string::npos && keywrd.find(" GEO_REF") != std::string::npos) {
        fout << "<TD>Toggle display<br><a href=\"javascript:Jmol.script(jmolApplet0,'\n";
    } else {
        fout << "<TD><a href=\"javascript:Jmol.script(jmolApplet0,'\n";
    }
    fout << "if (isOK1);  Display *; zoom 0; isOK1 = FALSE; else hide *; \n";
    line = " ";
    biggest_res = 0;
    for (i = 1; i <= nres; ++i) {
        n_res = std::string(11, ' ');
        n_res.replace(0, 5, "[" + res_txt[i].substr(0, 3) + "]");
        j = (int)(reada(res_txt[i], 4) + 0.5);
        biggest_res = std::max(biggest_res, j);
        if (j < 0) {
            num = std::string(1, (char)('1' + (int)std::log10(-j + 1.0)));
            char b[24];
            std::snprintf(b, sizeof(b), "%*d%2s", num[0] - '0', -j, res_txt[i].substr(7, 2).c_str());
            n_res.replace(5, (int)std::strlen(b), b);
        } else {
            num = std::string(1, (char)('2' + (int)std::log10(j + 1.0)));
            char b[24];
            std::snprintf(b, sizeof(b), "%*d%2s", num[0] - '0', j, res_txt[i].substr(7, 2).c_str());
            n_res.replace(5, (int)std::strlen(b), b);
        }
        wrt_res = n_res.substr(1, 3) + n_res.substr(5);
        if (wrt_res[1] >= 'A' && wrt_res[1] <= 'Z') wrt_res[1] = (char)(wrt_res[1] - 'A' + 'a');
        if (wrt_res[2] >= 'A' && wrt_res[2] <= 'Z') wrt_res[2] = (char)(wrt_res[2] - 'A' + 'a');
        j = len_trim(line);
        if (j > 120) {
            j = 0;
            fout << trim(line) << "\n";
            line = " ";
        }
        l_res = res_txt[i].substr(0, 7) + res_txt[i].substr(8, 1);
        if (l_res[0] >= '0' && l_res[0] <= '9')
            l_res[0] = (char)('A' + (l_res[0] - '0'));
        k = (int)l_res.find('-');
        if (k > 0) l_res[k] = '_';
        line += std::string(std::max(0, j + 2 - (int)line.size()), ' ') + l_res + " = FALSE;";
    }
    fout << trim(line) << "\n";
    fout << "isOK1 = TRUE; isOK2 = FALSE; lzoom = TRUE; lcenter = TRUE;\n";
    if (keywrd.find(" 0SCF") != std::string::npos && keywrd.find(" GEO_REF") != std::string::npos) {
        fout << "endif;')\"> 1&2</a>\n";
        fout << "<a href=\"javascript:Jmol.script(jmolApplet0,'if (isOKone);  Display add */1.1; zoom 0; isOKone = FALSE;"
                "else hide add */1.1; zoom 0; isOKone = TRUE;  lzoom = TRUE; lcenter = TRUE; endif;')\"> 1</a>\n";
        fout << "<a href=\"javascript:Jmol.script(jmolApplet0,'if (isOKtwo);  Display add */2.1; zoom 0; isOKtwo = FALSE;"
                "else hide add */2.1; zoom 0; isOKtwo = TRUE;  lzoom = TRUE; lcenter = TRUE; endif;')\"> 2</a>\n";
    } else {
        fout << "endif;')\">Toggle display all</a>\n";
    }
    fout << "</TD>\n";
    fout << "<TD><a href=\"javascript:Jmol.script(jmolApplet0,'\n";
    fout << "if (lcenter);  lcenter = FALSE; else lcenter = TRUE; center {visible}; end if;\n";
    fout << "')\">Toggle center picture</a> </TD>\n</TR> <TR>\n";
    fout << "<TD><a href=\"javascript:Jmol.script(jmolApplet0,'console;')\">Console</a> </TD> \n";
    fout << " \n";
    fout << "<TD> <a href=\"javascript:Jmol.script(jmolApplet0,'if (lzoom); zoom 0; "
            "select */2.1; color bonds green; select off; lzoom = FALSE; else lzoom = TRUE; "
            "end if ')\">Fit to screen</a> </TD></TR><TR>\n";
    fout << "<TD><a href=\"javascript:Jmol.script(jmolApplet0,'display within(3,visible);"
            "select */2.1; color bonds green; select off; if (lzoom); zoom 0; end if ')\">Near Neighbors</a> </TD> \n";
    fout << "<TD><a href=\"javascript:Jmol.script(jmolApplet0,'"
            "connect (all) (all) delete; connect; display add connected(visible); select *; hbonds calculate; select "
            " */2.1; color bonds green; select off; if (lzoom); delay 0.001; zoom 0; end if ')\">Connected Neighbors</a> </TD> \n";
    if (!p.empty()) {
        fout << "<TR><TD>\n";
        fout << "<a href=\"javascript:Jmol.script(jmolApplet0,'if (!lcharge_x); \n";
        fout << "var use = {visible}; var sel = {selected};\n";
        fout << "var z = 0; for (var i IN @sel){z = 3}\n";
        fout << "if (z = 3); use = sel; end if;\n";
        fout << "for (var x IN @use){select @x; var txt =  (x.temperature > 0 ? \\'+\\':\\'\\')"
                "+format(\\'%1.2f\\',x.temperature ); label @txt; color label black;\n";
        fout << "set labelOffset 0 0;}  select @sel; lcharge_x= TRUE;\n";
        fout << "else lcharge_x= FALSE; var use = {visible}; var sel = {selected};\n";
        fout << "var z = 0; for (var i IN @sel){z = 3}\n";
        fout << "if (z = 3); use = sel; end if;\n";
        fout << "select @use; label OFF; select @sel; end if ')\">Charges as Nos.</a>\n";
        fout << "</TD><TD>\n";
        fout << "<a href=\"javascript:Jmol.script(jmolApplet0,'if (!lcharge_s); \n";
        fout << "var use = {visible}; var sel = {selected};\n";
        fout << "var z = 0; for (var i IN @sel){z = 3}\n";
        fout << "if (z = 3); use = sel; end if;\n";
        fout << "for (var x IN @use)\n";
        fout << "{select @x; var txt =  @x.temperature*0.5;\n";
        fout << "if (@txt > 0){spacefill @txt; color atom deepskyblue;}\n";
        fout << "if (!@txt > 0){txt = -txt; spacefill @txt; color atom deeppink;}}\n";
        fout << "select @sel; lcharge_s= TRUE;\n";
        fout << "else lcharge_s= FALSE; var use = {visible}; var sel = {selected};\n";
        fout << "var z = 0; for (var i IN @sel){z = 3}\n";
        fout << "if (z = 3); use = sel; end if;\n";
        fout << "select @use; spacefill 15%; color cpk; select @sel; end if ')\">Charges as Sizes</a>\n";
        fout << "</TD></TR>\n";
    }
    fout << "<TR>\n";
    line = input_fn.substr(0, len_trim(input_fn) - 4) + "txt";
    line_1 = trim(line);
    for (char& cc : line_1) if (cc >= 'a' && cc <= 'z') cc = (char)(cc - 'a' + 'A');
    i = 0;
    for (j = 1; j <= len_trim(line) - 4; ++j)
        if (line[j - 1] == '/' || line[j - 1] == '\\') i = j;
    if (i != 0) line = line.substr(i);
    fout << "<TD colspan=\"2\"><a href=\"javascript:Jmol.script(jmolApplet0,'script common.txt;')\">"
            "<strong style=\"font-size:20px\">Common Script</strong>\n";
    fout << "</a>&nbsp; (Read file from:<br> \"<a href=\"common.txt\"  target=\"_blank\">common.txt</a>\")</TD>\n";
    fout << "</TR> <TR>\n";
    fout << "<TD colspan=\"2\"><a href=\"javascript:Jmol.script(jmolApplet0,'script \\'" << trim(line)
         << "\\';')\">"
            "<strong style=\"font-size:20px\">Specific Script</strong>\n";
    fout << "</a>&nbsp; (Read file from:<br> \"<a href=\"" << trim(line)
         << "\"  target=\"_blank\">" << trim(line) << "</a>\")</TD>\n";
    fout << "</TR></TABLE>\n";
    if (nres > 0) {
        fout << "<p align=\"center\">Toggle Individual Residues</p>\n";
        fout << "<TABLE>\n<TR>\n";
        ncol = std::max(7, nres / 16) + 1;
        nprt = 1;
        for (i = 1; i <= nres; ++i) {
            n_res = std::string(11, ' ');
            n_res.replace(0, 5, "[" + res_txt[i].substr(0, 3) + "]");
            j = (int)(reada(res_txt[i], 4) + 0.5);
            if (j < 0) {
                num = std::string(1, (char)('2' + (int)std::log10(-j + 1.0)));
                char b[24];
                std::snprintf(b, sizeof(b), "%*d%2s", num[0] - '0', j, res_txt[i].substr(7, 2).c_str());
                n_res.replace(5, (int)std::strlen(b), b);
            } else {
                num = std::string(1, (char)('1' + (int)std::log10(j + 1.0)));
                char b[24];
                std::snprintf(b, sizeof(b), "%*d%2s", num[0] - '0', j, res_txt[i].substr(7, 2).c_str());
                n_res.replace(5, (int)std::strlen(b), b);
            }
            wrt_res = n_res.substr(1, 3) + n_res.substr(5);
            l_res = res_txt[i].substr(0, 7) + res_txt[i].substr(8, 1);
            if (l_res[0] >= '0' && l_res[0] <= '9')
                l_res[0] = (char)('A' + (l_res[0] - '0'));
            for (j = 1; j <= 23; ++j)
                if (n_res.substr(1, 3) == tyres[j]) break;
            if (j > 23) {
                n_res.replace(0, 5, "     ");
                wrt_res[0] = 'X';
            } else {
                wrt_res[0] = tyr[j][0];
                if (wrt_res[0] == '?') wrt_res[0] = 'X';
            }
            k = (int)l_res.find('-');
            if (k > 0) l_res[k] = '_';
            fout << "<TD> <a href=\"javascript:Jmol.script(jmolApplet0,'if (!" << l_res
                 << ");   display ADD " << trim(n_res) << ";  " << l_res << " = TRUE; else \n";
            fout << " hide ADD " << trim(n_res) << ";  " << l_res << " = FALSE; end if; "
                 << "if (lcenter); center {visible}; end if; if (lzoom); zoom 0; end if;')\" class=\"auto-style5\">  \n";
            if (wrt_res.find(':') != std::string::npos) {
                j = (int)wrt_res.find(':') + 1;
                fout << "<p class=\"auto-style4\">" << wrt_res[0] << "<br>" << wrt_res.substr(3, j - 2 - 3) << "</p> </a></TD>\n";
            } else {
                fout << "<p class=\"auto-style4\">" << wrt_res[0] << "<br>" << wrt_res.substr(3) << "</p> </a></TD>\n";
            }
            if (nprt == ncol) {
                nprt = 1;
                fout << "</TR> <TR>\n";
            } else {
                nprt = nprt + 1;
            }
        }
        fout << "</TR>\n</TABLE><br> \n";
    }
    fout << " </TD><TD>\n";
    if (nres >= limres) {
        write_data_to_html(iprt);
    }
    fout << "<span id=mydiv></span><a href=\"javascript:Jmol.script(jmolApplet0)\"></a></TD></TABLE>\n";
    fout << "</BODY>\n</HTML>\n";
    fout.close();
    line = line.substr(0, len_trim(line) - 3);
    add_path(line);
    if (!std::ifstream(line + "txt").good()) {
        std::ofstream ftxt(line + "txt", std::ios::out | std::ios::trunc);
        if (ftxt) {
            i = 0;
            for (j = 1; j <= len_trim(line); ++j)
                if (line[j - 1] == '/' || line[j - 1] == '\\') i = j;
            if (i != 0) line = line.substr(i);
            ftxt << "#\n# Script for use with the HTML file \"" << trim(line) << "html\"\n#\n";
        }
    }
}

void write_data_to_html(int iprt) {
    using namespace molkst_C;
    using namespace common_arrays_C;
    using namespace parameters_C;
    (void)iprt;
    std::string idate = fdate();
    std::fprintf(stdout, "<TABLE>\n");
    std::fprintf(stdout, "<TR><TD> Date:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%s</TD></TR>\n",
                 (idate.substr(4, 7) + idate.substr(20) + idate.substr(10, 6)).c_str());
    std::fprintf(stdout, "<TR><TD> No. atoms:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%5d</TD></TR>\n", numat);
    bool store_log = chanel_C::log;
    chanel_C::log = false;
    empiri();
    chanel_C::log = store_log;
    size_t colon = formula.find(':') + 1;
    size_t eqls = formula.find('=');
    std::string line = " ";
    for (size_t ic = colon; ic < eqls; ++ic) {
        size_t jc = ic + 1;
        if ((formula[ic] < '0' || formula[ic] > '9') && formula[jc] >= '0' && formula[jc] <= '9') {
            line = trim(line) + formula[ic] + "<sub>";
        } else if ((formula[jc] < '0' || formula[jc] > '9') && formula[ic] >= '0' && formula[ic] <= '9') {
            line = trim(line) + formula[ic] + "</sub>";
        } else {
            line = trim(line) + formula[ic];
        }
    }
    std::fprintf(stdout, "<TR><TD> Formula:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%s</TD></TR>\n", trim(line).c_str());
    if (keywrd.find(" 0SCF") != std::string::npos && keywrd.find(" GEO_REF") != std::string::npos) {
        if (std::fabs(arc_hof_1) > 1.e-4)
            std::fprintf(stdout, "<TR><TD> Dataset:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%12.3f kcal/mol</TD></TR>\n", arc_hof_1);
        if (std::fabs(arc_hof_2) > 1.e-4)
            std::fprintf(stdout, "<TR><TD> GEO_REF:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%12.3f kcal/mol</TD></TR>\n", arc_hof_2);
        size_t irms = keywrd.find(" RMS_DIFF");
        if (irms != std::string::npos) {
            size_t if_ = keywrd.find("F=", irms);
            if (if_ != std::string::npos) {
                size_t sp = keywrd.find(' ', if_);
                std::fprintf(stdout, "<TR><TD> RMS Diff.:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%s&Aring;</TD></TR>\n",
                             keywrd.substr(if_ + 2, sp - if_ - 2).c_str());
            }
        }
    }
    if (line.find('H') != std::string::npos && nelecs > 0) {
        double sum = -nelecs;
        for (int i = 1; i <= numat; ++i) sum += tore[nat[i]];
        int ii = (int)(sum + 0.5);
        if (ii != 0) {
            std::fprintf(stdout, "<TR><TD> Net charge:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%+d</TD></TR>\n", ii);
        } else {
            std::fprintf(stdout, "<TR><TD> Net charge:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>Zero</TD></TR>\n");
        }
    }
    if (std::fabs(escf) > 1.e-10)
        std::fprintf(stdout, "<TR><TD> Heat of Formation:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>%12.3f Kcal/mol</TD></TR>\n", escf);
    std::fprintf(stdout, "</TABLE>\n");
}

void add_path(std::string& line) {
    // Fortran add_path prepends the working directory when the name is relative.
    // The C++ translation keeps the name as-is (files are opened in the CWD).
    (void)line;
}
