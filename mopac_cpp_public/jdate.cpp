// jdate.cpp
#define _CRT_SECURE_NO_WARNINGS
#include "jdate.h"
#include <ctime>
#include <cstdio>

std::string jdate() {
    static const int dim[12] = {0,31,28,31,30,31,30,31,31,30,31,30};
    std::time_t t = std::time(nullptr);
    std::tm* tm = std::localtime(&t);
    int iyear = tm->tm_year % 100;
    int imonth = tm->tm_mon + 1;
    int iday = tm->tm_mday;
    int ijulian = iday;
    for (int i = 1; i <= imonth; ++i) ijulian += dim[i-1];
    char buf[9];
    std::snprintf(buf, sizeof(buf), "%02d%03d", iyear, ijulian);
    return std::string(buf) + "   ";
}
