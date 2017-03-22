#!/bin/bash
#SBATCH -N 8
#SBATCH -p GPU
#SBATCH --gre=gpu:k80:4
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kgills@gmail.com
#SBATCH --ntasks-per-node=1

NODES=8
ITERS=$(seq 1 10)
SIZE=8192
FILE_NAME="floyd_acc_mpi_k80"

touch "${FILE_NAME}_${NODES}_${SIZE}.txt"
echo "matrix_dim, etime, flops, cores">>"${FILE_NAME}_${NODES}_${SIZE}.txt"

for ITER in ${ITERS}
do
    mpirun -n ${NODES} ./floyd.out>>"${FILE_NAME}_${NODES}_${SIZE}.txt"
done
