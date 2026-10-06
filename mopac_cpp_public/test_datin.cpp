// test_datin.cpp
#include <cstdio>
#include <vector>
#include <string>
#include "datin.h"
int main() {
    auto files = split_param_files(" EXTERNAL=foo.dat;bar.dat other");
    bool ok1 = (files.size() == 2 && files[0] == "foo.dat" && files[1] == "bar.dat");
    bool ok2 = (elemnt_sym(1) == "H " && elemnt_sym(8) == "O " && elemnt_sym(26) == "FE");
    std::printf("files=%zu elem=%s/%s/%s %s\n", files.size(),
                elemnt_sym(1).c_str(), elemnt_sym(8).c_str(), elemnt_sym(26).c_str(),
                (ok1 && ok2) ? "PASS" : "FAIL");
    return (ok1 && ok2) ? 0 : 1;
}
