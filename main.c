/*
*   Gauss's elimination method (para)
+
*   Student Number : 1w193077
*   Author         : T. Tabuchi
*   Last Update    : 2022/02/04
*   Since          : 2022/02/02
*
*   Compile        : $ mpicc main.c read_matrix_mtx.c LUSolver_parallel.c calc_residual_norm.c -o main_parallel
*   Execution      : $ mpirun -np 4 ./main_parallel
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <mpi.h>
#include "main.h"

double A[N+1][N+1], b[N+1], c[N+1], x[N+1];
double A0[N+1][N+1], b0[N+1]; // Copy the initial values of A and b for Relative residual norm

int main(int argc, char* argv[])
{
  int i, j, m, n, nz;
  double t, Mflops, tmp, rz;
  clock_t tic, toc;

  int myid, nodes;
  MPI_Status istatus;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);

  read_matrix_mtx(A, &m, &n, &nz, "bcsstk14.mtx"); // Loading data

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
	}

  /******************************/
  /* Measure the execution time */
  /******************************/
  tic = clock();
  LUSolver_parallel(A, b, c, x, m, n, myid, nodes, istatus);
  // mpi_sample(myid, nodes, istatus); // test for parallization
  toc = clock();
  /******************************/

  if (myid == 0) {
    t = (double)(toc - tic)/CLOCKS_PER_SEC;
    Mflops = (n*n*n/3+n*n)/t/1e+6;

    /* Calculate the relative residual norm */
    rz = calc_residual_norm(m, n, x, A0, b0);

  // for (i = 0; i < m; i++) {
  //   for (j = 0; j < n; j++) {
  //     printf("%.1f ", A[i][j]);
  //   }
  //   printf("\n");
  // }
  for (i = 0; i < 100; i++) {
    printf("%.1f\n", c[i]);
  }

    printf("=========================================\n");
    printf("Time                   : %.2f [s]\n", t);
    printf("Performance            : %.2f [MFLOPS]\n", Mflops);
    printf("Relative residual norm : %.2lg\n", rz);
    printf("=========================================\n\n");
  }

  // if (myid == 0) {
  //   for (i = 0; i < 100; i++) {
  //     printf("%.1f\n", c[i]);
  //   }
  // }

  MPI_Finalize();
	
  return 0;
}