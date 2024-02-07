/*
*   Gauss's elimination method (lu)
+
*   Student Number : 1w193077
*   Author         : T. Tabuchi
*   Last Update    : 2022/02/04
*   Since          : 2022/02/02
*
*   Compile        : $ gcc -O2 LUSolver.c read_matrix_mtx.c calc_residual_norm.c -o lu
*   Execution      : $ ./lu
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "LUSolver.h"

double A[N+1][N+1], b[N+1], c[N+1], x[N+1];
double A0[N+1][N+1], b0[N+1]; // Copy the initial values of A and b for Relative residual norm

void LUSolver(double A[N+1][N+1], double b[N+1], double c[N+1], double x[N+1], int m, int n) 
{
  int i, j, k;
  double dtmp, dakj;

  // A = LU
  /********************/
  /* LU Decomposition */
  /********************/
  for (k = 0; k < n; k++) {
    dtmp = 1.0 / A[k][k];
    for (i = k+1; i < n; i++) {
      A[i][k] = A[i][k]*dtmp;
    }
    for (j = k+1; j < n; j++) {
      dakj = A[k][j];
      for (i = k+1; i < n; i++) {
        A[i][j] = A[i][j] - A[i][k]*dakj;
      }
    }
  }

  // Lc = b
  /*************************/
  /* Forward substitution */
  /*************************/
  for (k = 0; k < n; k++) {
    c[k] += b[k];
    for (j = 0; j < k; j++) {
      c[k] -= A[k][j]*c[j];
    }
  }

  // Ux = c
  /*************************/
  /* Backward substitution */
  /*************************/
  for (k = n-1; k >= 0; k--) {
    x[k] += c[k];
    for (j = n-1; j > k; j--) {
      x[k] -= A[k][j]*x[j];
    }
    x[k] /= A[k][k];
  }
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
  LUSolver(A, b, c, x, m, n);
  toc = clock();
  /******************************/

  t = (double)(toc - tic)/CLOCKS_PER_SEC;
  Mflops = (n*n*n/3+n*n)/t/1e+6;

  /* Calculate the relative residual norm */
  rz = calc_residual_norm(m, n, x, A0, b0);

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