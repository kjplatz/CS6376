#!/bin/bash
##Author: Anay Sonawane
##File Name: job_accmpi.sh
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 28

#echo commands to stdout
set -x

module unload icc

module load pgi
module load mpi/pgi_openmpi

mpicc -acc -ta=tesla:cuda8.0 floyd_accmpi.c

 for var in 1 2 4 5 8 10 16 20 
 do
   mpirun -n $var ./a.out 2000
 done


