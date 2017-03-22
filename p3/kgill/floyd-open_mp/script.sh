#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kgills@gmail.com

ITERS=$(seq 1 1)
CORES=$(seq 19 28)
SIZE=8192
FILE_NAME="floyd_omp"

for CORE in ${CORES}
do

    export OMP_NUM_THREADS=${CORE}
    touch "${FILE_NAME}_${CORE}_${SIZE}.txt"
    echo "matrix_dim, etime, flops, cores">>"${FILE_NAME}_${CORE}_${SIZE}.txt"

    for ITER in ${ITERS}
    do
        ./floyd_omp.out>>"${FILE_NAME}_${CORE}_${SIZE}.txt"
    done
done
