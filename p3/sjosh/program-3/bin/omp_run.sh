#! /bin/bash

#SBATCH -N 28
#SBATCH -p RM
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=shriroop.joshi@gmail.com

mkdir -p ./output
rm -f ./output/*

for((i = 2; i <= 28; i += 2)); do
    export OMP_NUM_THREADS=$i
    ./omp_floyd > ./output/output-$i
done
