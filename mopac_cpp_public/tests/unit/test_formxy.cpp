// test_formxy.cpp
#include <cstdio>
#include <vector>
#include "formxy.h"
int main() {
    int kr=0, na=1, nb=1;
    std::vector<double> w={0,1.0}, wca={0,0}, wcb={0,0};
    std::vector<double> ca={0,2.0}, cb={0,3.0};
    formxy(w,kr,wca,wcb,ca,cb,na,nb);
    bool ok=(wca[1]==0.75 && wcb[1]==0.5 && kr==1);
    std::printf("formxy wca=%g wcb=%g kr=%d %s\n",wca[1],wcb[1],kr, ok?"PASS":"FAIL");
    return ok?0:1;
}
