#!/bin/bash
##Author: Anay Sonawane
##File Name: job_openmpi.sh
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 28

#echo commands to stdout
set -x

MAX_NUM_THREADS=28;
NUM_THREADS=1;

module unload icc

module load pgi
module load mpi/pgi_openmpi

mpicc -mp floyd_openmpi.c

while [ $NUM_THREADS -lt $MAX_NUM_THREADS ]
do

 NUM_THREADS=`expr $NUM_THREADS + 1`
 export OMP_NUM_THREADS=$NUM_THREADS;
 for var in 1 2 4 5 8 10 16 20 
 do
   mpirun -n $var ./a.out 2000
 done
done


