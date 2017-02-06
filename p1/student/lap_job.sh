#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 00:05:00
export OMP_NUM_THREADS=28
./stripe
