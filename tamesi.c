#include <stdio.h>
#include <stdlib.h>

#define N 10

void func(int A[N+1][N+1]) {
  A[0][0] = 1;
}

int main(void) {

  int i, j, m, n;
  double ans;

  int A[N+1][N+1];

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i][j] = 0;
    }
  }

  // func(A);

  for (i = 0; i < N; i+=2) {
    printf("%d\n", i);
  }

  // printf("%d\n", 10/3);

  return 0;
}