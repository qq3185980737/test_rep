// test_ef.cpp
#include <cstdio>
#include <vector>
#include "ef.h"
int main() {
    std::vector<double> g={1,0}, oldf={0,0}, tv={0.5,0}, sv={0,0}, d={1,0};
    std::vector<std::vector<double>> H(2,std::vector<double>(2,0.0));
    updhes(2,g,oldf,tv,sv,d,H);
    bool ok=(H[0][0]==0.5);
    std::printf("updhes BFGS H00=%g %s\n",H[0][0], ok?"PASS":"FAIL");
    return ok?0:1;
}
