#include <stdio.h>
#include <math.h>
#include "main.h"
// #include "LUSolver.h"
// #include "gauss.h"

double calc_residual_norm(int m, int n, double x[N+1], double A0[N+1][N+1], double b0[N+1]) 
{
  int i, j;
  double tmp, btmp = 0.0, z = 0.0, rz;

  for (i = 0; i < n; i++) {
    tmp = 0.0;
    for (j = 0; j < n; j++) {
      tmp += A0[i][j] * x[j];
    }
    tmp = tmp - b0[i];
    z += tmp*tmp;
  }
  z = sqrt(z);

  for (i = 0; i < n; i++) {
    btmp += x[i]*x[i];
  }
  btmp = sqrt(btmp);
  rz = z / btmp;

  return rz;
}