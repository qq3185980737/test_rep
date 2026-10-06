// myword.cpp — C++ translation of MOPAC 2016 "myword.F90".
// F90: j = index(keywrd,testwd); do while keywrd(j:j)==' ': j=j+1  (check-then-advance).
#include "myword.h"
bool myword(std::string& keywrd, const std::string& testwd) {
    bool res=false;
    size_t len_key=keywrd.find_last_not_of(' ');
    if (len_key==std::string::npos) return false;
    bool quote=false;
    size_t j;
    while ((j=keywrd.find(testwd))!=std::string::npos) {
        while (j<keywrd.size() && keywrd[j]==' ') j++;
        res=true;
        size_t k;
        for (k=j;k<=len_key;++k) {
            if (keywrd[k]=='"') quote=!quote;
            if (!quote) {
                if (keywrd[k]=='=' || keywrd[k]==' ') {
                    if (keywrd[k]=='=') { keywrd[k]=' '; continue; }
                    size_t j2;
                    for (j2=k+1;j2<=len_key;++j2) {
                        if (keywrd[j2]=='=') { keywrd[j2]=' '; break; }
                        if (keywrd[j2]!=' ') break;
                    }
                    break;
                }
            }
            keywrd[k]=' ';
        }
    }
    return res;
}
