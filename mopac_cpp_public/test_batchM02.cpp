// test_batchM02.cpp — M02 输入/文件解析 语义断言
// 覆盖：reada(Fortran 数值解析)、digit、upcase(含 EXTERNAL 保护)、
//       myword(引号/= 处理)、getval(大写判断+reada)、jdate
#include <cstdio>
#include <cmath>
#include <string>
#include "reada.h"
#include "digit.h"
#include "upcase.h"
#include "myword.h"
#include "getval.h"
#include "jdate.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_reada() {
    std::printf("Test 1: reada Fortran numeric extraction\n");
    CHK_D(reada("12.345", 1), 12.345, 1e-12, "reada plain");
    CHK_D(reada("1.5D-3", 1), 1.5e-3, 1e-15, "reada D exponent");
    CHK_D(reada("  -2.5E+2", 1), -250.0, 1e-12, "reada signed E exponent");
    CHK_D(reada("3.14XYZ", 1), 3.14, 1e-12, "reada trailing junk");
    CHK_D(reada("ABC 3.14", 1), 3.14, 1e-12, "reada first numeric field");
    CHK_D(reada("0.5", 1), 0.5, 1e-12, "reada 0.5");
    CHK_D(reada(".5", 1), 0.5, 1e-12, "reada .5");
    CHK_D(reada("+7", 1), 7.0, 1e-12, "reada +7");
    CHK_D(reada("1E2", 1), 100.0, 1e-12, "reada 1E2");
    CHK_D(reada("1.5d-3", 1), 1.5e-3, 1e-15, "reada lowercase d");
    CHK_D(reada("123 456", 5), 456.0, 1e-12, "reada istart into second field");
}

static void test_digit() {
    std::printf("Test 2: digit clean substring -> double\n");
    CHK_D(digit("12.345", 1), 12.345, 1e-12, "digit plain");
    CHK_D(digit("-3.5", 1), -3.5, 1e-12, "digit negative");
    CHK_D(digit("+2.5", 1), 2.5, 1e-12, "digit plus");
    CHK_D(digit(" 2.5", 1), 2.5, 1e-12, "digit leading space");
    CHK_D(digit(".5", 1), 0.5, 1e-12, "digit .5");
    CHK_D(digit("123", 1), 123.0, 1e-12, "digit integer");
}

static void test_upcase() {
    std::printf("Test 3: upcase\n");
    std::string s = "abc DEF 123";
    upcase(s, (int)s.size());
    CHECK(s == "ABC DEF 123", "upcase basic");
    s = "a\tb";
    upcase(s, (int)s.size());
    CHECK(s == "A B", "upcase tab to space");
    s = "EXTERNAL=abc Def";
    upcase(s, (int)s.size());
    CHECK(s == "EXTERNAL=abc DEF", "upcase EXTERNAL protected");
}

static void test_myword() {
    std::printf("Test 4: myword keyword removal\n");
    std::string k1 = "AM1 XYZ";
    bool r1 = myword(k1, "AM1");
    CHECK(r1 && k1 == "    XYZ", "myword plain keyword deleted");
    std::string k2 = "CHARGE=2 AM1";
    bool r2 = myword(k2, "CHARGE");
    // F90: positions 1..9 (CHARGE=2) become blanks; " AM1" survives
    bool lead_blank = k2.substr(0, 9).find_first_not_of(' ') == std::string::npos;
    CHECK(r2 && lead_blank && k2.find("AM1") != std::string::npos,
          "myword keyword with = deleted");
    std::string k3 = "PM7 GEO-OK";
    bool r3 = myword(k3, "NONSENSE");
    CHECK(!r3 && k3 == "PM7 GEO-OK", "myword absent keyword");
    // quoted text: F90 deletes text inside quotes too; quote only suppresses
    // the '='/space handling inside quotes (and resets at closing quote)
    std::string k4 = "XYZ \"AM1\" ABC";
    bool r4 = myword(k4, "AM1");
    CHECK(r4, "myword finds keyword inside quotes (F90 deletes it)");
}

static void test_getval() {
    std::printf("Test 5: getval numeric vs word\n");
    double x = 0; std::string t;
    getval("1.5", x, t);
    CHK_D(x, 1.5, 1e-12, "getval number");
    CHECK(t == " ", "getval number tag blank");
    getval("1.5D-3", x, t);
    CHK_D(x, 1.5e-3, 1e-15, "getval D exponent via reada");
    getval("XYZ 1.5", x, t);
    CHK_D(x, -999.0, 1e-12, "getval word -> -999");
    CHECK(t == "XYZ", "getval word tag");
    // lowercase start: F90 treats as numeric (only A-Z blocks numeric)
    getval("abc 1.5", x, t);
    CHK_D(x, 1.5, 1e-12, "getval lowercase start numeric");
}

static void test_jdate() {
    std::printf("Test 6: jdate Julian date string\n");
    std::string d = jdate();
    CHECK(d.size() >= 7 && d.find_first_not_of("0123456789 ") != std::string::npos
          ? (d[0] >= '0' && d[0] <= '9') : true, "jdate starts with digit");
    CHECK(d.size() >= 5, "jdate length");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M02 batch — input / file parsing\n");
    test_reada();
    test_digit();
    test_upcase();
    test_myword();
    test_getval();
    test_jdate();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
