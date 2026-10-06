// upcase.cpp
#include "upcase.h"
#include <cctype>
void upcase(std::string& s, int n) {
    std::string buf=s;
    for (int i=0;i<n;++i) {
        char c=s[i];
        if (c>='a'&&c<='z') s[i]=char(c-'a'+'A');
        else if (c=='\t') s[i]=' ';
    }
    size_t p=s.find("EXTERNAL");
    if (p!=std::string::npos) {
        size_t q=s.find(' ',p+1);
        if (q==std::string::npos) q=s.size();
        for (size_t k=p+9;k<q;++k) s[k]=buf[k];
    }
}
