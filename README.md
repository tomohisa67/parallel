
# LU分解法の並列化

# Requirement

* MPICH 4.0

# Usage

```bash
$ mpicc main.c read_matrix_mtx.c LUSolver_parallel.c calc_residual_norm.c -o main_parallel
$ mpirun -np 4 ./main_parallel
```

# Note
* gauss.c              : 以前作成した部分軸選択ありのガウスの消去法
* LUSolver.c           : 単純なLU分解
* LUSolver_parallel.c  : LUSolverの並列化バージョン(未完)
* read_matrix_mtx.c    : 行列の読み込み
* calc_residual_norm.c : 相対残差ノルムを求める