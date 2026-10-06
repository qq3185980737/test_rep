// fdate.cpp — Fortran fdate intrinsic: current date/time string (writmo dep).
#include "fdate.h"
#include <ctime>
#include <string>

void fdate(std::string& date_str) {
  std::time_t t = std::time(nullptr);
  std::tm tm{};
  localtime_s(&tm, &t);
  char buf[40];
  std::strftime(buf, sizeof(buf), "%a %b %e %H:%M:%S %Y", &tm);
  date_str = buf;
}
