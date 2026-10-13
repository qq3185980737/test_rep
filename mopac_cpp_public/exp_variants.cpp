// exp_variants.cpp — test small variants of the LIBM exp to match FTN on
// characteristic x: block19 (ok), block41 (need ...67), block11 (need ...D0)
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>

typedef std::uint32_t juint;
static const juint _cv[] = {
    0x652b82feUL, 0x40571547UL, 0x652b82feUL, 0x40571547UL, 0xfefa0000UL,
    0x3f862e42UL, 0xfefa0000UL, 0x3f862e42UL, 0xbc9e3b3aUL, 0x3d1cf79aUL,
    0xbc9e3b3aUL, 0x3d1cf79aUL, 0xfffffffeUL, 0x3fdfffffUL, 0xfffffffeUL,
    0x3fdfffffUL, 0xe3289860UL, 0x3f56c15cUL, 0x555b9e25UL, 0x3fa55555UL,
    0xc090cf0fUL, 0x3f811115UL, 0x55548ba1UL, 0x3fc55555UL
};
static const juint _mmask[] = { 0xffffffc0UL, 0x00000000UL };
static const juint _bias[]  = { 0x0000ffc0UL, 0x00000000UL };
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

static double ld(const juint* p) {
    std::uint64_t u = (std::uint64_t)p[1] << 32 | p[0];
    double d; std::memcpy(&d, &u, 8); return d;
}
static double nextup(double d) {  // +1 ULP
    std::uint64_t u; std::memcpy(&u, &d, 8);
    if (d > 0) u += 1; else u -= 1;
    double r; std::memcpy(&r, &u, 8); return r;
}
static double nextdn(double d) {
    std::uint64_t u; std::memcpy(&u, &d, 8);
    if (d > 0) u -= 1; else u += 1;
    double r; std::memcpy(&r, &u, 8); return r;
}
static void hex(const char* tag, double v) {
    std::uint64_t u; std::memcpy(&u, &v, 8);
    std::printf("%s %016llx\n", tag, (unsigned long long)u);
}

static const double SHIFTER = 6755399441055744.0;

// variant: ln2_hi ulp shift (+1/-1/0), m_round (0=nearest-even,1=floor), rmode (0=two-step,1=sum-first,2=fma,3=fma2), morder (0-7), tail (0=mul+add,1=(s6+1)*tlo)
double vexp(double x, int hi_ulp, int m_round, int rmode, int morder, int tail) {
    double log2e_hi = ld(_cv);
    double ln2_64_hi = ld(_cv + 4);
    double ln2_64_lo = ld(_cv + 8);
    double c2 = ld(_cv + 12);
    double c3_0 = ld(_cv + 16), c3_1 = ld(_cv + 18);
    double c4_0 = ld(_cv + 20), c4_1 = ld(_cv + 22);
    if (hi_ulp > 0) ln2_64_hi = nextup(ln2_64_hi);
    if (hi_ulp < 0) ln2_64_hi = nextdn(ln2_64_hi);

    double t1 = x * log2e_hi + SHIFTER;
    double m  = t1 - SHIFTER;
    if (m_round == 1) { m = std::floor(x * log2e_hi); }
    double r;
    if (rmode == 0) { r = x - ln2_64_hi * m; r = r - ln2_64_lo * m; }
    if (rmode == 1) { r = x - (ln2_64_hi * m + ln2_64_lo * m); }
    if (rmode == 2) { r = std::fma(-ln2_64_hi, m, x); r = r - ln2_64_lo * m; }
    if (rmode == 3) { r = std::fma(-ln2_64_hi, m, x); r = std::fma(-ln2_64_lo, m, r); }
    if (rmode == 4) { r = x - ln2_64_hi * m; }
    if (rmode == 5) { r = std::fma(-ln2_64_hi, m, x); }
    if (rmode == 6) { r = x - (ln2_64_hi * m); r = r - ln2_64_lo * m; double rr = r; (void)rr; }

    std::uint64_t bits; std::memcpy(&bits, &t1, 8);
    int j = (int)(bits & 63);
    int n = (int)(std::uint32_t)bits >> 6;
    std::uint64_t expbits = ((bits & 0xffffffc0ULL) + 0xffc0ULL) << 46;
    if (m_round != 0) {  // recompute expbits from m directly
        std::uint64_t mb; std::memcpy(&mb, &m, 8);
        // rebuild t1-equivalent low bits from integer m: t1_bits_low32 = (uint32_t)m
        std::uint32_t mlow = (std::uint32_t)(std::int64_t)m;
        expbits = ((std::uint64_t)(mlow & 0xffffffc0U) + 0xffc0ULL) << 46;
    }

    const juint* T = _Tbl_addr + 4 * j;
    double T_hi = ld(T);
    double T_lo = ld(T + 2);

    double p0 = c3_0 * r, p1 = c3_1 * r;
    double q0 = c4_0 + p0, q1 = c4_1 + p1;
    double rr2 = r * r, rr3 = r * rr2, rr5 = rr3 * rr2;
    double s1 = r + T_hi;
    double s2 = rr5 * q0;
    double s4 = rr3 * q1;
    double rr2c = rr2 * c2;
    double s6;
    switch (morder) {
        case 0: s6 = (s1 + s2) + s4 + rr2c; break;           // ((s1+s2)+s4)+rr2c (left assoc) == current
        case 1: s6 = s1 + (s2 + s4) + rr2c; break;
        case 2: s6 = (s1 + s4) + (s2 + rr2c); break;
        case 3: s6 = s1 + (s4 + (s2 + rr2c)); break;
        case 4: s6 = (s1 + (s4 + rr2c)) + s2; break;
        case 5: s6 = s1 + ((s2 + s4) + rr2c); break;
        case 6: s6 = (s1 + s2) + (s4 + rr2c); break;
        case 7: s6 = s1 + (s2 + (s4 + rr2c)); break;
        default: s6 = (s1 + s2) + s4 + rr2c;
    }

    std::uint64_t tlo_bits = ((std::uint64_t)T[3] << 32) | T[2];
    double tlo_scaled; std::uint64_t scaled = tlo_bits | expbits;
    std::memcpy(&tlo_scaled, &scaled, 8);
    if (tail == 0) return s6 * tlo_scaled + tlo_scaled;
    return (s6 + 1.0) * tlo_scaled;
}

int main() {
    double x19 = 0; { std::uint64_t u = 0x4016f7dc1cb0f397ULL; std::memcpy(&x19, &u, 8); }
    double x41 = 0; { std::uint64_t u = 0x402FE014FD2EED6FULL; std::memcpy(&x41, &u, 8); }
    double x11 = 0; { std::uint64_t u = 0x3fe4d23aa95c1f91ULL; std::memcpy(&x11, &u, 8); }
    const char* target19 = "3f6a480c0ad2d7a1";
    const char* target41 = "3e8013664404ca67";
    const char* target11 = "3ffeab41bff3c6d0";

    std::printf("targets: 19=%s 41=%s 11=%s\n", target19, target41, target11);
    // phase 1: rmode sweep (hi=0 mr=0 morder=0 tail=0)
    for (int rm = 0; rm <= 6; rm++) {
        double e19 = vexp(-x19, 0, 0, rm, 0, 0);
        double e41 = vexp(-x41, 0, 0, rm, 0, 0);
        double e11 = vexp(x11, 0, 0, rm, 0, 0);
        std::printf("rm=%d  19=", rm); hex("", e19);
        std::printf("    41="); hex("", e41);
        std::printf("    11="); hex("", e11);
    }
    // phase 2: morder sweep (rm=0)
    for (int mo = 0; mo <= 7; mo++) {
        double e19 = vexp(-x19, 0, 0, 0, mo, 0);
        double e41 = vexp(-x41, 0, 0, 0, mo, 0);
        double e11 = vexp(x11, 0, 0, 0, mo, 0);
        std::printf("mo=%d  19=", mo); hex("", e19);
        std::printf("    41="); hex("", e41);
        std::printf("    11="); hex("", e11);
    }
    // phase 3: tail sweep (rm=0 mo=0)
    for (int tl = 0; tl <= 1; tl++) {
        double e19 = vexp(-x19, 0, 0, 0, 0, tl);
        double e41 = vexp(-x41, 0, 0, 0, 0, tl);
        double e11 = vexp(x11, 0, 0, 0, 0, tl);
        std::printf("tl=%d  19=", tl); hex("", e19);
        std::printf("    41="); hex("", e41);
        std::printf("    11="); hex("", e11);
    }
    return 0;
}
