/*
*   Gauss's elimination method (lu)
+
*   Student Number : 1w193077
*   Author         : T. Tabuchi
*   Last Update    : 2022/02/04
*   Since          : 2022/02/04
*
*   Compile        : $ gcc -O2 gauss.c read_matrix_mtx.c calc_residual_norm.c -o gauss
*   Execution      : $ ./gauss
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "gauss.h"

double A[N+1][N+1], b[N+1], c[N+1], x[N+1];
double A0[N+1][N+1], b0[N+1]; // Copy the initial values of A and b for Relative residual norm
int s[N+1]; // 交換情報

void left_func(double A[N+1][N+1], double b[N+1], int m, int n) {
  int i, j, k, p;
  double amax, tmp;
  int tmpi;

  // sの初期化
  for (i = 0; i < n; i++) {
    s[i] = i;
  }

  // 部分軸選択
  for (k = 0; k < n-1; k++) {
    amax = fabs(A[k][k]);
    p = k;
    for (i = k+1; i < n; i++) {
      if(fabs(A[i][k]) > amax) {
        amax = fabs(A[i][k]);
        p = i;
      }
    }
    // amaxが0なら終了する
    if (amax < 1.0e-12) {
      printf("Error. pivot is 0.\n");
      exit(0);
    }
    // 第k行と第p行の交換
    if (p != k) {
      for (i = k; i < n; i++) {
        tmp = A[k][i]; // 並列化可能?
        A[k][i] = A[p][i];
        A[p][i] = tmp;
      }
      s[k] = p;
    }

    // 前進消去
    for (i = k+1; i < n; i++) {
      A[i][k] = - A[i][k] / A[k][k]; // alphaの代わりにA[i][k]を使う
      for (j = k+1; j < n; j++) {
        A[i][j] = A[i][j] + A[i][k] * A[k][j]; // 外積形式?
      }
    }
  }
}

void right_func(double A[N+1][N+1], double b[N+1], int m, int n) {
  int i, j, k;
  int l;
  double tmp;

  // 前進消去
  for (k = 0; k < m-1; k++) {
    // b[k]とb[l]の交換
    l = s[k];
    tmp = b[k];
    b[k] = b[l];
    b[l] = tmp;
    for (i = k+1; i < m; i++) {
      b[i] = b[i] + A[i][k] * b[k];
    }
  }

  // 後退代入
  for (i = m-1; i >= 0; i--) {
    for (j = i+1; j < n; j++) {
      b[i] = b[i] - A[i][j] * b[j];
      // A[i][j] = 0.0;
    }
    b[i] = b[i] / A[i][i];
    // A[i][i] = 1.0;
  }
}

void gauss(double A[N+1][N+1], double b[N+1], int m, int n)
{
  left_func(A, b, m, n);
  right_func(A, b, m, n);
}


int main(void)
{
  int i, j, m, n, nz;
  double t, Mflops, tmp, rz;
  clock_t tic, toc;

  read_matrix_mtx(A, &m, &n, &nz, "bcsstk14.mtx");

  m -= 2; n -= 2; nz -= 0; // 端数処理を避けるため

  /* Store an integer of 0 or 1 in b */
  srand(100);
	for (i = 0; i < m; i++) {
    tmp = (double)(rand()%2);
    b[i] = tmp;
    b0[i] = tmp;
		for (j = 0; j < n; j++) {
      A0[i][j]= A[i][j];
		}
    c[i] = 0.0;
    x[i] = 0.0;
	}

  /******************************/
  /* Measure the execution time */
  /******************************/
  tic = clock();
  gauss(A, b, m, n);
  toc = clock();
  /******************************/

  t = (double)(toc - tic)/CLOCKS_PER_SEC;
  Mflops = (n*n*n/3+n*n)/t/1e+6;

  /* Calculate the relative residual norm */
  rz = calc_residual_norm(m, n, b, A0, b0);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("%.1f ", A[i][j]);
    }
    printf("\n");
  }

  printf("=========================================\n");
  printf("Time                   : %.2f [s]\n", t);
  printf("Performance            : %.2f [MFLOPS]\n", Mflops);
  printf("Relative residual norm : %.2lg\n", rz);
  printf("=========================================\n\n");

  return 0; 
}