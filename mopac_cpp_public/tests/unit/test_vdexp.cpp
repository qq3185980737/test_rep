// test_vdexp.cpp — verify MKL vdExp matches Intel Fortran exp bit-for-bit
#include <mkl.h>
#include <mkl_vml.h>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>

static void hex(const char* tag, double x) {
    std::uint64_t ux;
    std::memcpy(&ux, &x, 8);
    std::printf("%s %016llx\n", tag, (unsigned long long)ux);
}

extern "C" double __svml_exp1(double);   // Intel SVML scalar exp

int main() {
    // aintgs block41: alpha=15.937660133341266, FTN c=exp(-alpha)=3E8013664404CA67
    std::uint64_t ua = 0x402fe014fd2eed6full;
    double a;
    std::memcpy(&a, &ua, 8);
    double na = -a;
    double c_std = std::exp(na);
    double c_vml;
    vdExp(1, &na, &c_vml);
    double c_recip = 1.0 / std::exp(a);
    double c_recip_vml;
    vdExp(1, &a, &c_recip_vml);
    c_recip_vml = 1.0 / c_recip_vml;
    double c_exp2 = std::exp2(na * 1.4426950408889634);  // log2(e)
    std::printf("aintgs alpha = %.17g\n", a);
    hex("  std::exp(-a)   ", c_std);
    hex("  vdExp(-a)      ", c_vml);
    hex("  1/std::exp(a)  ", c_recip);
    hex("  1/vdExp(a)     ", c_recip_vml);
    hex("  exp2(-a*log2e) ", c_exp2);
    std::printf("  target(FTN)    = 3E8013664404CA67\n");

    // bintgs block11: beta=0.6506627376374913, FTN expx=exp(beta)=3FFEAB41BFF3C6D0
    std::uint64_t ub = 0x3fe4d23aa95c1f91ull;
    double b;
    std::memcpy(&b, &ub, 8);
    double b_std = std::exp(b);
    double b_vml;
    vdExp(1, &b, &b_vml);
    double b_recip = 1.0 / std::exp(-b);
    std::printf("bintgs beta  = %.17g\n", b);
    hex("  std::exp(b)  ", b_std);
    hex("  vdExp(b)     ", b_vml);
    hex("  1/std::exp(-b)", b_recip);
    std::printf("  target(FTN) = 3FFEAB41BFF3C6D0\n");

    // also test expm1 / pow alternatives quickly
    double p_std = std::pow(-b, 4);
    double p_vml;
    vdPowx(1, &b, 4.0, &p_vml);  // vdPowx(n, a, b, y): y=a^b
    std::printf("pow test x=%.17g\n", b);
    hex("  std::pow(x,4) ", p_std);
    hex("  vdPowx(x,4)  ", p_vml);
    // negative base: vdPowx requires a>=0? test a=-b
    double nb = -b;
    double p_vml2;
    vdPowx(1, &nb, 4.0, &p_vml2);
    hex("  vdPowx(-x,4) ", p_vml2);

    // block19: alpha=5.742050598430885, FTN c=exp(-alpha)=3F6A480C0AD2D7A1
    std::uint64_t u19 = 0x4016f7dc1cb0f397ull;
    double x19;
    std::memcpy(&x19, &u19, 8);
    double c19_std = 1.0 / std::exp(x19);
    double c19_vml;
    vdExp(1, &x19, &c19_vml);
    double c19_rec = 1.0 / c19_vml;
    double c19_direct = std::exp(-x19);
    std::printf("block19 x = %.17g\n", x19);
    hex("  1/std::exp(x) ", c19_std);
    hex("  1/vdExp(x)    ", c19_rec);
    hex("  std::exp(-x)  ", c19_direct);
    hex("  __svml_exp1(-x)", __svml_exp1(-x19));
    hex("  1/__svml_exp1(x)", 1.0 / __svml_exp1(x19));
    std::printf("  target(FTN)   = 3F6A480C0AD2D7A1\n");
    // also block41
    hex("  __svml_exp1(-a)", __svml_exp1(-a));
    hex("  1/__svml_exp1(a)", 1.0 / __svml_exp1(a));
    return 0;
}
