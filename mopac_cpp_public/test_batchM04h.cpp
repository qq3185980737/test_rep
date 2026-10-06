// test_batchM04h.cpp  M04 bonds_for_MOZYME (bond orders/valencies) batch.
// MOZYME-style scene with two 1-orbital atoms: valency formula
// 2p-p^2 per atom, bond sum = sum p^2 over the i-j block, threshold 0.01.
// Output text is captured to a file (stdout redirection) and asserted.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>

#include "bonds_for_MOZYME.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;
using namespace MOZYME_C;

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg)                                     \
  do {                                                    \
    if (c) {                                              \
      ++n_pass;                                           \
      std::printf("  [PASS] %s\n", msg);                  \
    } else {                                              \
      ++n_fail;                                           \
      std::printf("  [FAIL] %s\n", msg);                  \
    }                                                     \
  } while (0)

static std::string run_and_capture() {
  std::fflush(stdout);
  int saved = _dup(_fileno(stdout));
  int fd = _open("_m04h_out.txt", _O_CREAT | _O_WRONLY | _O_TRUNC,
                  _S_IREAD | _S_IWRITE);
  if (fd < 0) return "";
  _dup2(fd, _fileno(stdout));
  _close(fd);
  bonds_for_MOZYME();
  std::fflush(stdout);
  _dup2(saved, _fileno(stdout));
  _close(saved);
  std::ifstream in("_m04h_out.txt");
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

int main() {
  std::printf("M04h batch - bonds_for_MOZYME\n");
  numat = 2;
  keywrd = " ";
  iorbs.assign({0, 1, 1});
  lijbo = true;
  nijbo.assign(3, std::vector<int>(3, 0));
  nijbo[1][1] = 0; nijbo[1][2] = 1; nijbo[2][1] = 1; nijbo[2][2] = 2;
  p.assign(4, 0.0);
  p[1] = 0.8;   // atom 1 diagonal (kk=ijbo(1,1)+1 = 1)
  p[2] = 0.5;   // 1-2 block (kl=2, ku=2)
  p[3] = 0.7;   // atom 2 diagonal (kk=ijbo(2,2)+1 = 3)
  nat.assign({0, 6, 6});  // C, C

  std::string out = run_and_capture();
  // valenc(1) = 2*0.8 - 0.8^2 = 0.96; valenc(2) = 2*0.7 - 0.7^2 = 0.91.
  // bond 1-2 sum = 0.5^2 = 0.25 > 0.01 -> printed as C=0.250.
  CHECK(out.find("valenc=0.960") != std::string::npos,
        "bonds valenc(1)=0.96");
  CHECK(out.find("valenc=0.910") != std::string::npos,
        "bonds valenc(2)=0.91");
  CHECK(out.find("C =0.250") != std::string::npos, "bonds order 0.25 printed");
  // Sum of squares: 0.5^2 = 0.25 exactly.
  CHECK(out.find("0.250") != std::string::npos, "bonds square sum 0.25");
  // Threshold: a 0.01-sum bond is NOT reported. Set p[2]=0.0999.
  p[2] = 0.0999;
  out = run_and_capture();
  CHECK(out.find("C =0.010") == std::string::npos,
        "bonds below threshold filtered");
  std::remove("_m04h_out.txt");
  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
