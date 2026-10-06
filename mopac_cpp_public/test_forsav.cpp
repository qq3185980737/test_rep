// test_forsav.cpp
#include <cstdio>
#include "forsav.h"
int main() {
    std::vector<std::vector<double>> dd(3,std::vector<double>(1,0));
    std::vector<double> fm(1,0), co(1,0), ev(1,0), fc(1,0);
    double t=0, refh=0; int ipt=0,js=0;
    forsav(t,dd,ipt,fm,co,1,refh,ev,js,fc);
    std::printf("forsav links OK PASS\n");
    return 0;
}
