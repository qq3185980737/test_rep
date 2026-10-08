#include <cstdio>
int main(){
  char buf[32];
  std::snprintf(buf, sizeof buf, " CHARGE=%2d", -1);
  std::printf("field=[%s] len=%zu\n", buf, __builtin_strlen(buf));
  return 0;
}
