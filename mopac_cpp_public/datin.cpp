// datin.cpp — C++ translation of MOPAC 2016 "datin.F90".
#include "datin.h"

#include <string>
#include <vector>

#include "chanel_C.h"
#include "molkst_C.h"

using namespace molkst_C;

namespace {
const char* elemnt[108] = {
    "",
    "H ","HE","LI","BE","B ","C ","N ","O ","F ","NE",
    "NA","MG","AL","SI","P ","S ","CL","AR","K ","CA",
    "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN",
    "GA","GE","AS","SE","BR","KR","RB","SR","Y ","ZR",
    "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN",
    "SB","TE","I ","XE","CS","BA","LA","CE","PR","ND",
    "PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
    "LU","HF","TA","W ","RE","OS","IR","PT","AU","HG",
    "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
    "PA","U ","NP","PU","AM","CM","BK","MI","XX","FM",
    "MD","CB","++","+","--","-","TV"
};
}

std::string elemnt_sym(int z) {
    if (z < 1 || z > 107) return "??";
    return elemnt[z];
}

std::vector<std::string> split_param_files(const std::string& kw) {
    std::vector<std::string> out;
    auto pos = kw.find("EXTERNAL=");
    if (pos == std::string::npos) pos = kw.find("PARAMS=");
    if (pos == std::string::npos) return out;
    pos = kw.find('=', pos);
    if (pos == std::string::npos) return out;
    // collect until next whitespace
    size_t end = kw.find(' ', pos);
    std::string seg = (end == std::string::npos) ? kw.substr(pos + 1) : kw.substr(pos + 1, end - pos - 1);
    size_t start = 0;
    while (true) {
        size_t semi = seg.find(';', start);
        if (semi == std::string::npos) { out.push_back(seg.substr(start)); break; }
        out.push_back(seg.substr(start, semi - start));
        start = semi + 1;
    }
    return out;
}

void datin(int) {
    // File-I/O body: open parameter files, parse FN*/PAR lines, call update().
    // Stubbed pending upcase/reada/update/add_path/write_params porting.
}
