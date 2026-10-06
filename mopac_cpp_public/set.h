// set.h — C++ translation.
#pragma once
namespace overlaps_C {
    extern thread_local int isp, ips;
    extern thread_local double sa, sb;
    extern thread_local double a[64], b[64];
    extern double fact[18];
}
void set_fn(double s1, double s2, int na, int nb, double rab, int ii);
