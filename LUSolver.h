#ifndef SUB_H
#define SUB_H

#define N 2000
#define MM_MAX_LINE_LENGTH 1025

void read_matrix_mtx(double A[N+1][N+1], int *m, int *n, int *nz, char *filename);
double calc_residual_norm(int m, int n, double x[N+1], double A0[N+1][N+1], double b0[N+1]);

#endif