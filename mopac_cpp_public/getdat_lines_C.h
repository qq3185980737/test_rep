// getdat_lines_C.h — shared input line queue for the translated MOPAC driver.
// getdat() fills this queue from the .mop/.dat data-set; gettxt()/getgeo() and
// other record readers consume it in order (replacing Fortran's unit "ir" reads).
#pragma once
#include <cstddef>
#include <string>
#include <vector>

extern std::vector<std::string> getdat_lines;  // one record per data-set line
extern std::size_t getdat_line_idx;            // next unconsumed record index
