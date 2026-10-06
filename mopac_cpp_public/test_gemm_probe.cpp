// Minimal dgemm_ correctness probe: 1-based column-major layout like our data
#include <cstdio>
#include <vector>
extern "C" void dgemm_(const char* ta, const char* tb, const int* m, const int* n,
  const int* k, const double* alpha, const double* a, const int* lda,
  const double* b, const int* ldb, const double* beta, double* c, const int* ldc);
int main() {
  int n = 3;
  // 1-based column-major: element(row=j, col=i) at linear j+(i-1)*n
  std::vector<double> a(10), b(10), c(10);
  // A(row=j,col=i) = (i-1)*3 + j  -> A = [1 2 3 / 4 5 6 / 7 8 9] row-major display
  // B(row=j,col=i) = (j-1)*3 + i  -> B = [1 4 7 / 2 5 8 / 3 6 9]
  for (int i = 1; i <= n; i++) for (int j = 1; j <= n; j++) {
    a[j + (i - 1) * n] = (i - 1) * 3 + j;
    b[j + (i - 1) * n] = (j - 1) * 3 + i;
  }
  std::printf("A:\n");
  for (int j = 1; j <= n; j++) { for (int i = 1; i <= n; i++) std::printf(" %g", a[j + (i - 1) * n]); std::printf("\n"); }
  std::printf("B:\n");
  for (int j = 1; j <= n; j++) { for (int i = 1; i <= n; i++) std::printf(" %g", b[j + (i - 1) * n]); std::printf("\n"); }
  char N = 'N'; double one = 1.0, zero = 0.0;
  dgemm_(&N, &N, &n, &n, &n, &one, &a[1], &n, &b[1], &n, &zero, &c[1], &n);
  std::printf("C = A*B (expect [30 36 42 / 66 81 96 / 102 126 150]):\n");
  for (int j = 1; j <= n; j++) { for (int i = 1; i <= n; i++) std::printf(" %g", c[j + (i - 1) * n]); std::printf("\n"); }
  return 0;
}
