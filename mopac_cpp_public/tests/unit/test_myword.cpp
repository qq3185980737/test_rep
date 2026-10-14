// test_myword.cpp
#include <cstdio>
#include "myword.h"
int main() {
    std::string kw="PM6 EPS=78.4  FORCE";
    bool r1=myword(kw,"EPS");
    bool r2=myword(kw,"NO");
    std::printf("r1=%d r2=%d kw=\"%s\"\n",(int)r1,(int)r2,kw.c_str());
    return (r1&&!r2)?0:1;
}
