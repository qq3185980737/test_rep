// libm_exp.cpp — bit-exact reimplementation of Intel LIBM exp
// (from OpenJDK stubGenerator_x86_64_exp.cpp, GPLv2, Intel copyright),
// translated instruction-by-instruction to scalar double.
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <cstdlib>

typedef std::uint32_t juint;

// Tables verbatim from the OpenJDK source.
static const juint _cv[] = {
    0x652b82feUL, 0x40571547UL, 0x652b82feUL, 0x40571547UL, 0xfefa0000UL,
    0x3f862e42UL, 0xfefa0000UL, 0x3f862e42UL, 0xbc9e3b3aUL, 0x3d1cf79aUL,
    0xbc9e3b3aUL, 0x3d1cf79aUL, 0xfffffffeUL, 0x3fdfffffUL, 0xfffffffeUL,
    0x3fdfffffUL, 0xe3289860UL, 0x3f56c15cUL, 0x555b9e25UL, 0x3fa55555UL,
    0xc090cf0fUL, 0x3f811115UL, 0x55548ba1UL, 0x3fc55555UL
};
static const juint _mmask[] = { 0xffffffc0UL, 0x00000000UL, 0xffffffc0UL, 0x00000000UL };
static const juint _bias[]  = { 0x0000ffc0UL, 0x00000000UL, 0x0000ffc0UL, 0x00000000UL };
static const juint _Tbl_addr[] = {
    0x00000000UL, 0x00000000UL, 0x00000000UL, 0x00000000UL, 0x0e03754dUL,
    0x3cad7bbfUL, 0x3e778060UL, 0x00002c9aUL, 0x3567f613UL, 0x3c8cd252UL,
    0xd3158574UL, 0x000059b0UL, 0x61e6c861UL, 0x3c60f74eUL, 0x18759bc8UL,
    0x00008745UL, 0x5d837b6cUL, 0x3c979aa6UL, 0x6cf9890fUL, 0x0000b558UL,
    0x702f9cd1UL, 0x3c3ebe3dUL, 0x32d3d1a2UL, 0x0000e3ecUL, 0x1e63bcd8UL,
    0x3ca3516eUL, 0xd0125b50UL, 0x00011301UL, 0x26f0387bUL, 0x3ca4c554UL,
    0xaea92ddfUL, 0x0001429aUL, 0x62523fb6UL, 0x3ca95153UL, 0x3c7d517aUL,
    0x000172b8UL, 0x3f1353bfUL, 0x3c8b898cUL, 0xeb6fcb75UL, 0x0001a35bUL,
    0x3e3a2f5fUL, 0x3c9aecf7UL, 0x3168b9aaUL, 0x0001d487UL, 0x44a6c38dUL,
    0x3c8a6f41UL, 0x88628cd6UL, 0x0002063bUL, 0xe3a8a894UL, 0x3c968efdUL,
    0x6e756238UL, 0x0002387aUL, 0x981fe7f2UL, 0x3c80472bUL, 0x65e27cddUL,
    0x00026b45UL, 0x6d09ab31UL, 0x3c82f7e1UL, 0xf51fdee1UL, 0x00029e9dUL,
    0x720c0ab3UL, 0x3c8b3782UL, 0xa6e4030bUL, 0x0002d285UL, 0x4db0abb6UL,
    0x3c834d75UL, 0x0a31b715UL, 0x000306feUL, 0x5dd3f84aUL, 0x3c8fdd39UL,
    0xb26416ffUL, 0x00033c08UL, 0xcc187d29UL, 0x3ca12f8cUL, 0x373aa9caUL,
    0x000371a7UL, 0x738b5e8bUL, 0x3ca7d229UL, 0x34e59ff6UL, 0x0003a7dbUL,
    0xa72a4c6dUL, 0x3c859f48UL, 0x4c123422UL, 0x0003dea6UL, 0x259d9205UL,
    0x3ca8b846UL, 0x21f72e29UL, 0x0004160aUL, 0x60c2ac12UL, 0x3c4363edUL,
    0x6061892dUL, 0x00044e08UL, 0xdaa10379UL, 0x3c6ecce1UL, 0xb5c13cd0UL,
    0x000486a2UL, 0xbb7aafb0UL, 0x3c7690ceUL, 0xd5362a27UL, 0x0004bfdaUL,
    0x9b282a09UL, 0x3ca083ccUL, 0x769d2ca6UL, 0x0004f9b2UL, 0xc1aae707UL,
    0x3ca509b0UL, 0x569d4f81UL, 0x0005342bUL, 0x18fdd78eUL, 0x3c933505UL,
    0x36b527daUL, 0x00056f47UL, 0xe21c5409UL, 0x3c9063e1UL, 0xdd485429UL,
    0x0005ab07UL, 0x2b64c035UL, 0x3c9432e6UL, 0x15ad2148UL, 0x0005e76fUL,
    0x99f08c0aUL, 0x3ca01284UL, 0xb03a5584UL, 0x0006247eUL, 0x0073dc06UL,
    0x3c99f087UL, 0x82552224UL, 0x00066238UL, 0x0da05571UL, 0x3c998d4dUL,
    0x667f3bccUL, 0x0006a09eUL, 0x86ce4786UL, 0x3ca52bb9UL, 0x3c651a2eUL,
    0x0006dfb2UL, 0x206f0dabUL, 0x3ca32092UL, 0xe8ec5f73UL, 0x00071f75UL,
    0x8e17a7a6UL, 0x3ca06122UL, 0x564267c8UL, 0x00075febUL, 0x461e9f86UL,
    0x3ca244acUL, 0x73eb0186UL, 0x0007a114UL, 0xabd66c55UL, 0x3c65ebe1UL,
    0x36cf4e62UL, 0x0007e2f3UL, 0xbbff67d0UL, 0x3c96fe9fUL, 0x994cce12UL,
    0x00082589UL, 0x14c801dfUL, 0x3c951f14UL, 0x9b4492ecUL, 0x000868d9UL,
    0xc1f0eab4UL, 0x3c8db72fUL, 0x422aa0dbUL, 0x0008ace5UL, 0x59f35f44UL,
    0x3c7bf683UL, 0x99157736UL, 0x0008f1aeUL, 0x9c06283cUL, 0x3ca360baUL,
    0xb0cdc5e4UL, 0x00093737UL, 0x20f962aaUL, 0x3c95e8d1UL, 0x9fde4e4fUL,
    0x00097d82UL, 0x2b91ce27UL, 0x3c71affcUL, 0x82a3f090UL, 0x0009c491UL,
    0x589a2ebdUL, 0x3c9b6d34UL, 0x7b5de564UL, 0x000a0c66UL, 0x9ab89880UL,
    0x3c95277cUL, 0xb23e255cUL, 0x000a5503UL, 0x6e735ab3UL, 0x3c846984UL,
    0x5579fdbfUL, 0x000a9e6bUL, 0x92cb3387UL, 0x3c8c1a77UL, 0x995ad3adUL,
    0x000ae89fUL, 0xdc2d1d96UL, 0x3ca22466UL, 0xb84f15faUL, 0x000b33a2UL,
    0xb19505aeUL, 0x3ca1112eUL, 0xf2fb5e46UL, 0x000b7f76UL, 0x0a5fddcdUL,
    0x3c74ffd7UL, 0x904bc1d2UL, 0x000bcc1eUL, 0x30af0cb3UL, 0x3c736eaeUL,
    0xdd85529cUL, 0x000c199bUL, 0xd10959acUL, 0x3c84e08fUL, 0x2e57d14bUL,
    0x000c67f1UL, 0x6c921968UL, 0x3c676b2cUL, 0xdcef9069UL, 0x000cb720UL,
    0x36df99b3UL, 0x3c937009UL, 0x4a07897bUL, 0x000d072dUL, 0xa63d07a7UL,
    0x3c74a385UL, 0xdcfba487UL, 0x000d5818UL, 0xd5c192acUL, 0x3c8e5a50UL,
    0x03db3285UL, 0x000da9e6UL, 0x1c4a9792UL, 0x3c98bb73UL, 0x337b9b5eUL,
    0x000dfc97UL, 0x603a88d3UL, 0x3c74b604UL, 0xe78b3ff6UL, 0x000e502eUL,
    0x92094926UL, 0x3c916f27UL, 0xa2a490d9UL, 0x000ea4afUL, 0x41aa2008UL,
    0x3c8ec3bcUL, 0xee615a27UL, 0x000efa1bUL, 0x31d185eeUL, 0x3c8a64a9UL,
    0x5b6e4540UL, 0x000f5076UL, 0x4d91cd9dUL, 0x3c77893bUL, 0x819e90d8UL,
    0x000fa7c1UL
};
static const juint _ALLONES[] = { 0xffffffffUL, 0xffffffffUL, 0xffffffffUL, 0xffffffffUL };
static const juint _ebias[] = { 0x00000000UL, 0x3ff00000UL, 0x00000000UL, 0x3ff00000UL };
static const juint _XMAX[] = { 0xffffffffUL, 0x7fefffffUL };
static const juint _XMIN[] = { 0x00000000UL, 0x00100000UL };
static const juint _INF[]  = { 0x00000000UL, 0x7ff00000UL };

static double ld(const juint* p) {  // little-endian u64 -> double
    std::uint64_t u = (std::uint64_t)p[1] << 32 | p[0];
    double d; std::memcpy(&d, &u, 8); return d;
}
static double ld_hi(const juint* p) {  // 8 bytes at p interpreted LE
    return ld(p);
}
static std::uint64_t lu(const juint* p) {
    return (std::uint64_t)p[1] << 32 | p[0];
}
static void hex(const char* tag, double v) {
    std::uint64_t u; std::memcpy(&u, &v, 8);
    std::printf("%s %016llx\n", tag, (unsigned long long)u);
}

// SHIFTER = 2^52 * 1.5
static const double SHIFTER = 6755399441055744.0;
static const double ONE = 1.0;

double libm_exp(double x) {
    // --- constants; NOTE: cv+16 etc are BYTE offsets -> juint index = offset/4
    double log2e_hi = ld(_cv);            // cv+0   = 64*log2e = 92.33248261688866
    double ln2_64_hi = ld(_cv + 4);       // cv+16  = ln2/64 hi = 0.010830424...
    double ln2_64_lo = ld(_cv + 8);       // cv+32  = ln2/64 lo (residual)
    double c2        = ld(_cv + 12);      // cv+48  = 0.49999999999999994
    double c3_0      = ld(_cv + 16);      // cv+64  = lane0 c3
    double c3_1      = ld(_cv + 18);      // cv+72  = lane1 c3
    double c4_0      = ld(_cv + 20);      // cv+80  = lane0 c4
    double c4_1      = ld(_cv + 22);      // cv+88  = lane1 c4

    // xmm1 = x*log2e_hi ; xmm1 += SHIFTER ; m = xmm1 - SHIFTER
    double t1 = x * log2e_hi + SHIFTER;
    double m  = t1 - SHIFTER;
    // r = x - ln2_64_hi*m ; then r -= ln2_64_lo*m
    double r  = x - ln2_64_hi * m;
    r = r - ln2_64_lo * m;

    // j = low 6 bits of t1 low32 ; n = t1 low32 >> 6 (arithmetic)
    std::uint64_t bits; std::memcpy(&bits, &t1, 8);
    int j = (int)(bits & 63);
    int n = (int)(std::uint32_t)bits >> 6;   // low 32 bits, arithmetic shift
    // 2^n bit pattern: ((t1_bits & 0xffffffc0) + 0xffc0) << 46
    std::uint64_t expbits = ((bits & 0xffffffc0ULL) + 0xffc0ULL) << 46;

    // T[j]: lane0 = T_hi, lane1 = T_lo
    const juint* T = _Tbl_addr + 4 * j;
    double T_hi = ld(T);
    double T_lo = ld(T + 2);

    // poly: xmm4=[c3_0*r, c3_1*r]; xmm5=[c4_0+c3_0*r, c4_1+c3_1*r]
    double p0 = c3_0 * r;              // xmm4 lane0
    double p1 = c3_1 * r;              // xmm4 lane1
    double q0 = c4_0 + p0;             // xmm5 lane0
    double q1 = c4_1 + p1;             // xmm5 lane1
    double rr2 = r * r;                // xmm6 = r^2
    double rr3 = r * rr2;              // xmm0 = r^3
    double rr5 = rr3 * rr2;            // mulsd xmm0, xmm6 -> r^5 (lane0); lane1 stays r^3
    double s1  = r + T_hi;             // addsd xmm1, xmm2
    double s2  = rr5 * q0;             // mulpd xmm0, xmm5 lane0
    double s3  = s1 + s2;              // addsd xmm1, xmm0
    double s4  = rr3 * q1;             // lane1 of xmm0
    double s5  = s3 + s4;              // addsd xmm0, xmm1
    double s6  = s5 + rr2 * c2;        // addsd xmm0, xmm6

    // 2^n * T_lo
    std::uint64_t tlo_bits = lu(T + 2);
    double tlo_scaled;
    std::uint64_t scaled = tlo_bits | expbits;
    std::memcpy(&tlo_scaled, &scaled, 8);

    double res = s6 * tlo_scaled + tlo_scaled;   // mulsd+addsd
    return res;
}

#ifdef LIBMMAIN_TEST
int main() {
    // sanity: table j=0
    double T0hi = ld(_Tbl_addr), T0lo = ld(_Tbl_addr + 2);
    std::printf("T[0] hi=%.17g lo=%.17g\n", T0hi, T0lo);
    std::uint64_t u; std::memcpy(&u, &T0lo, 8);
    std::printf("T[0] lo bits=%016llx\n", (unsigned long long)u);
    // targets
    // exact hex inputs from FTN probes
    double x19 = 0; { std::uint64_t u = 0x4016f7dc1cb0f397ULL; std::memcpy(&x19, &u, 8); }
    double x41 = 0; { std::uint64_t u = 0x402FE014FD2EED6FULL; std::memcpy(&x41, &u, 8); }
    double x11 = 0; { std::uint64_t u = 0x3fe4d23aa95c1f91ULL; std::memcpy(&x11, &u, 8); }
    hex("exp(-19) libm", libm_exp(-x19));
    std::printf("  target(-19) = 3f6a480c0ad2d7a1\n");
    // detailed intermediates for block41
    _putenv_s("LMDBG", "1");
    hex("exp(-41) libm", libm_exp(-x41));
    std::printf("  target(-41) = 3e8013664404ca67\n");
    _putenv_s("LMDBG", "");
    hex("exp(+025) libm", libm_exp(x11));
    std::printf("  target(+025)= 3ffeab41bff3c6d0\n");
    hex("exp(-41) std ", std::exp(-x41));
    hex("exp(-19) std ", std::exp(-x19));
    hex("exp(+025) std", std::exp(x11));
    return 0;
}
#endif  // LIBMMAIN_TEST
