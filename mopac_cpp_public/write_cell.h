// write_cell.h — C++ translation of "write_cell" (moldat.F90).
#pragma once
// write_cell(iprt): print unit-cell data and formula-unit properties.
// iprt < 0 or GUI: no output; iprt == 0: only print if nstep changed.
void write_cell(int iprt);
