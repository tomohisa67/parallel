#include <stdio.h>
#include <stdlib.h>
#include "main.h"
// #include "LUSolver.h"
// #include "gauss.h"

void read_matrix_mtx(double A[N+1][N+1], int *m, int *n, int *nz, char *filename)
{
  int i, j, count=0;
  double val;
  char line[MM_MAX_LINE_LENGTH];
  FILE *fp;
  fp = fopen(filename, "r");

  /* Aの全ての要素に0を格納 */
  for (i = 0; i < N+1; i++) {
    for (j = 0; j < 0; j++) {
      A[i][j] = 0.0;
    }
  }

  /* 1行読み込み, 先頭が'%'なら繰り返す */
  do {
    fgets(line, MM_MAX_LINE_LENGTH, fp);
  } while (line[0] == '%');

  /* 行列サイズと非ゼロ要素数をm,n,nzに格納 */
  sscanf(line, "%d %d %d", m, n, nz);

  /* 最後まで読み込んで, 非ゼロ要素をAに格納 */
  while(fscanf(fp, "%d %d %lg", &i, &j, &val) != EOF) {
    i--; 
    j--;
    // A[i][j] = val; // 対称行列であることを利用する
    // A[j][i] = val; 

    if (i <= 1803 && j <= 1803) { // 端数処理用
      A[i][j] = val; // 対称行列であることを利用する
      A[j][i] = val; 
      if (i == j && A[i][j] < 0.000001) {
        printf("Error.\n");
        exit(0);
      }
    } else {
      count++;
      // printf("(%d,%d) count: %d\n", i, j, count);
    }
  }

  fclose(fp);
}