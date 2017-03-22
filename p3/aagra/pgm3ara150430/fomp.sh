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
./floyd_omp

export OMP_NUM_THREADS=3
./floyd_omp

export OMP_NUM_THREADS=4
./floyd_omp

export OMP_NUM_THREADS=5
./floyd_omp

export OMP_NUM_THREADS=6
./floyd_omp

export OMP_NUM_THREADS=8
./floyd_omp

export OMP_NUM_THREADS=12
./floyd_omp

export OMP_NUM_THREADS=16
./floyd_omp

export OMP_NUM_THREADS=20
./floyd_omp

export OMP_NUM_THREADS=24
./floyd_omp

export OMP_NUM_THREADS=28
./floyd_omp
