// web_message.cpp — print a pointer to the online manual page.
#include "web_message.h"

#include <cstdio>

// 1-arg form used by getgeo.cpp (defaults to the standard output unit).
void web_message(const char* page) {
    web_message(6, page);
}

void web_message(int iw, const char* page) {
    std::printf("\n%10sFor more information, see: HTTP://OpenMOPAC.net/Manual/%s\n\n",
                "", page ? page : "");
}
