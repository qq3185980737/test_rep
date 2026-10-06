// test_dernvo.cpp
#include <cstdio>
#include "dernvo.h"
#include "meci_C.h"
#include "molkst_C.h"
int main() {
    using namespace meci_C;
    using namespace molkst_C;
    norbs=4; nclose=1; nopen=1; nmos=1; nelec=0;
    dernvo();
    bool ok = (nbo[1]==1 && nbo[2]==0 && nbo[3]==3);
    std::printf("nbo=(%d,%d,%d) %s\n", nbo[1],nbo[2],nbo[3], ok?"PASS":"FAIL");
    return ok?0:1;
}
