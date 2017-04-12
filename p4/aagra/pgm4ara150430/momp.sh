#!/bin/bash
#SBATCH -N 28
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 5:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=ayush.agrawal464@gmail.com
# echo commands to stdout 
set -x


# run OpenMP program
export OMP_NUM_THREADS=2
./matrix_omp

export OMP_NUM_THREADS=3
./matrix_omp

export OMP_NUM_THREADS=4
./matrix_omp

export OMP_NUM_THREADS=5
./matrix_omp

export OMP_NUM_THREADS=6
./matrix_omp

export OMP_NUM_THREADS=8
./matrix_omp

export OMP_NUM_THREADS=12
./matrix_omp

export OMP_NUM_THREADS=16
./matrix_omp

export OMP_NUM_THREADS=20
./matrix_omp

export OMP_NUM_THREADS=24
./matrix_omp

export OMP_NUM_THREADS=28
./matrix_omp
