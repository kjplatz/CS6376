#!/bin/bash
##Author:Anay Sonawane
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

mpicc -mp -O1 cannon_openmpi.c

while [ $NUM_THREADS -lt $MAX_NUM_THREADS ]
do

 NUM_THREADS=`expr $NUM_THREADS + 1`
 export OMP_NUM_THREADS=$NUM_THREADS;
 for var in 1 4 16 
 do
	for N in 256 512 1000 1024 2000 2048 3000 4096 8192 16384 
	do	
		mpirun -n $var ./a.out $N 
   	done
 done
done

