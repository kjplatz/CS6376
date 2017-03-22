#!/bin/bash
#SBATCH -p RM
#SBATCH -t 0:10:00
#SBATCH -N 28
#SBATCH --ntasks-per-node 28

#echo commands to stdout
set -x

module unload icc

module load pgi
module load mpi/pgi_openmpi

mpicc -O floyd_mpi.c

mpirun -np 2 ./a.out

mpirun -np 4 ./a.out

mpirun -np 10 ./a.out

mpirun -np 20 ./a.out

