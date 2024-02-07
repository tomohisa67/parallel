#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <mpi.h>
#include "main.h"

void left_func_parallel(int ib, int istart, int iend, int idiagPE, int myid, int nodes, MPI_Status istatus) 
{
  int i, j, k, p;
  double amax, tmp;
  int tmpi;

  for (k = 0; k < n-1; k++) {
    /************* 部分軸選択 **********************/
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
    /********************************************/

    /************* 前進消去 **********************/
    for (i = k+1; i < n; i++) {
      A[i][k] = - A[i][k] / A[k][k]; // alphaの代わりにA[i][k]を使う
      for (j = k+1; j < n; j++) {
        A[i][j] = A[i][j] +  A[k][j] * A[i][k]; // 外積形式?
      }
    }
    /********************************************/
  }
}

void right_func_parallel(int ib, int istart, int iend, int idiagPE, int myid, int nodes, MPI_Status istatus) 
{
  int i, j, k;
  int l;
  double tmp;

  /************* 前進消去 **********************/
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
  /********************************************/

  /************* 後退代入 **********************/
  for (i = m-1; i >= 0; i--) {
    for (j = i+1; j < n; j++) {
      b[i] = b[i] - A[i][j] * b[j];
      A[i][j] = 0.0;
    }
    b[i] = b[i] / A[i][i];
    A[i][i] = 1.0;
  }
  /********************************************/
}

void gauss_parallel(double A[N+1][N+1], double b[N+1], int m, int n, int myid, int nodes, MPI_Status istatus)
{
  int ib, istart, iend, idiagPE;
  double buf[N+1], c[N+1];

  ib = n/nodes; // 端数処理は?
  istart = myid * ib;
  iend = (myid+1) + ib;

  left_func_parallel(ib, istart, iend, idiagPE, myid, nodes, istatus);
  right_func_parallel(ib, istart, iend, idiagPE, myid, nodes, istatus);

}