// test_convert_storage.cpp
#include <cstdio>
#include <vector>

#include "convert_storage.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

int main() {
    using namespace common_arrays_C;
    using namespace molkst_C;
    using namespace MOZYME_C;
    numat = 1; norbs = 2; mpack = 3; nelecs = 2;
    nfirst = {0, 1, 2}; nlast = {0, 1, 2};
    std::vector<double> pk(4, 1.0), tri;
    convert_mat_packed_to_triangle(pk, tri);
    std::vector<std::vector<double>> csq(3, std::vector<double>(3, 9.0));
    isort.clear();
    ncocc.assign(3, 0); ncvir.assign(3, 0); ncf.assign(3, 0); nce.assign(3, 0);
    nncf.assign(3, 0); nnce.assign(3, 0);
    convert_lmo_packed_to_square(csq);
    bool ok = (tri.size() == 4 && csq[1][1] == 0.0 && isort[1] == 1 && isort[2] == 1);
    std::printf("tri.size=%zu c[1][1]=%g isort=(%d,%d) %s\n",
                tri.size(), csq[1][1], isort[1], isort[2], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
