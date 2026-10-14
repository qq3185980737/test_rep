// test_findn1.cpp
#include <cstdio>
#include <vector>
#include "findn1.h"
#include "common_arrays_C.h"
int main() {
    using namespace common_arrays_C;
    // 6-atom ring: 1-2-3-4-5-6-1
    nbonds.assign(7,0); ibonds.assign(7,std::vector<int>(7,0));
    auto add=[&](int a,int b){ ibonds[++nbonds[a]][a]=b; };
    add(1,2);add(1,6);
    add(2,1);add(2,3);
    add(3,2);add(3,4);
    add(4,3);add(4,5);
    add(5,4);add(5,6);
    add(6,5);add(6,1);
    bool ring = ring_6(1,2,6); // neighbors of 1: 2 and 6; ring closure via 3..5? 2->3,6->5,3->4,5->4 => kka=atom4=4 == atom3 path 3->4. true
    std::printf("ring_6(1,2,6)=%d %s\n",(int)ring, ring?"PASS":"FAIL");
    return ring?0:1;
}
