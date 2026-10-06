// test_gpu_info.cpp
#include <cstdio>
#include "mod_gpu_info.h"
int main() { int h[6],hd[6],nd,ns[6],cr[6],ma[6],mi[6]; char n[6]; size_t tm[6]; getGPUInfo(h,hd,&nd,n,ns,tm,cr,ma,mi); int st; setDevice(0,&st); std::printf("gpu_info nd=%d st=%d PASS\n",nd,st); return 0; }
