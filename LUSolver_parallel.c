#include <stdio.h>
#include <stdlib.h>
// #include <mpi.h>
#include "main.h"

void LUSolver_parallel(double A[N+1][N+1], double b[N+1], double c[N+1], double x[N+1], int m, int n, int myid, int nodes, MPI_Status istatus) 
{
  int i, j, k, kk, ib, istart, iend, idiagPE;
  double buf[N+1], dtemp;

  ib = n/nodes; 
  istart = myid * ib;
  iend = (myid+1) * ib;
  // if (myid == nodes-1) { // 端数処理 (n % nodes == 2)
  //   iend = n;
  // } else {
  //   iend = (myid+1) * ib;
  // }

  for (k = 0; k < m; k++) {
    buf[k] = 0.0;
  }

  /********************/
  /* LU Decomposition */
  /********************/
  for (k = 0; k < iend; k++) {
    idiagPE = k / ib; // kが0からibならば, idaigPEは0になる // 端数処理の時は注意!
    if (myid == idiagPE) { // 枢軸列を持つプロセス
      dtemp = 1.0 / A[k][k];

      /* 枢軸列の計算と, buf[]へ枢軸列をコピー */
      for (j = k+1; j < m; j++) { 
        buf[j] = A[j][k] * dtemp; // A[j][k] = A[j][k] * dtemp
      }

      for (i = myid+1; i < nodes; i++) { // 枢軸列の転送
        MPI_Send(&buf[k+1], m-k+1, MPI_DOUBLE, i, k, MPI_COMM_WORLD); // buf[k], ibでもいいけど, 冗長
      }
      istart = k+1; // 担当範囲の縮小 // istart += 1 でも意味は同じ あまり綺麗な書き方ではない?
    } else { // 枢軸列を持たないプロセス
      MPI_Recv(&buf[k+1], m-k+1, MPI_DOUBLE, idiagPE, k, MPI_COMM_WORLD, &istatus); // buf[k], ibでもいいけど, 冗長
    }
    /* 共通消去部分 */
    for (j = k+1; j < n; j++) {
      dtemp = buf[j];
      for (i = istart; i < iend; i++) {
        A[j][i] = A[j][i] - A[k][i]*dtemp;
      }
    }
  }

  // printf("break %d\n", myid); // ここまで到達していることを確認

  /* 前進消去にLが被らないように同期する */
  MPI_Barrier(MPI_COMM_WORLD);

  istart = myid * ib; 
  iend = (myid+1) * ib; // 担当範囲の初期化
  
  /*************************/
  /* Forward substitution */
  /*************************/
  for (k = 0; k < m; k++) {
    c[k] = 0.0;
  }

  for (k = 0; k < m; k += ib) { // 1ループ目 diagPE == 0, 2ループ目 1ループ目 diagPE == 1, ...
    if (k >= istart) { // プロセスの限定
      idiagPE = k / ib; // 端数処理の時は注意!
      /* ランク番号左隣りのプロセスからデータを受け取る */
      if (myid != 0) {
        MPI_Recv(&c[k], ib, MPI_DOUBLE, myid-1, k, MPI_COMM_WORLD, &istatus);
      } 
      if (myid == idiagPE) { // 対角ブロックをもつプロセス
        /* 対角ブロックだけ先に計算して値を確定させる */
        for (kk = 0; kk < ib; kk++) {
          c[k+kk] = b[k+kk] + c[k+kk]; // kはistartを表す. k+kkはistartからiendの間の値 送られてきた途中結果を加算する
          for (j = istart; j < istart+kk; j++) {
            c[k+kk] -= A[k+kk][j]*c[j];
          }
        }
      } else { // 対角ブロックを持たないプロセス
        /* 自分の所有範囲のデータのみ計算 (まだ最終結果ではない) */
        for (kk = 0; kk < ib; kk++) {
          for (j = istart; j < iend; j++) {
            c[k+kk] -= A[k+kk][j]*c[j];
          }
        }
        /* ランク番号右隣のプロセスに, 自分の担当範囲のデータを用いた演算結果を送る */
        if (myid != nodes-1) {
          MPI_Send(&c[k], ib, MPI_DOUBLE, myid+1, k, MPI_COMM_WORLD);
        }
      }
    }
  }

  /* 一旦同期する */
  MPI_Barrier(MPI_COMM_WORLD);

  /*************************/
  /* Backward substitution */
  /*************************/
  for (k = 0; k < m; k++) {
    x[k] = 0.0;
  }
  
  for (k = m-1; k >= 0; k -= ib) {
    if (k <= iend-1) {
      idiagPE = k / ib;
      if (myid != nodes-1) {
        MPI_Recv(&x[k-ib+1], ib, MPI_DOUBLE, myid+1, k, MPI_COMM_WORLD, &istatus); // 配列
      }
      if (idiagPE == myid) {
        for (kk = ib-1; kk >= 0; kk--) {
          x[k+kk] = c[k+kk] + x[k+kk];
          for (j = iend-1; j >= istart; j--) {
            x[k+kk] = x[k+kk] - A[k+kk][j]*c[j];
          }
          x[k+kk] = x[k+kk] / A[k+kk][k+kk];
        }
      } else {
        for (kk = ib-1; kk >= 0; kk--) {
          for (j = iend-1; j >= istart; j--) {
            x[k+kk] = x[k+kk] - A[k+kk][j]*c[j];
          }
          x[k+kk] = x[k+kk] / A[k+kk][k+kk];
        }
        if(myid != 0) {
          MPI_Send(&x[k-ib+1], ib, MPI_DOUBLE, myid-1, k, MPI_COMM_WORLD);
        }
      }
    }
  }

  /* 一旦同期する */
  MPI_Barrier(MPI_COMM_WORLD); // 必要ない

  // 結果をプロセス0のxに集約する
  if (myid == 0) {
    for (i = 1; i < nodes; i++) {
      MPI_Recv(&x[i*ib], ib, MPI_DOUBLE, i, i+9990, MPI_COMM_WORLD, &istatus);
    } 
  } else {
    MPI_Send(&x[myid*ib], ib, MPI_DOUBLE, 0, myid+9990, MPI_COMM_WORLD);
  }
}