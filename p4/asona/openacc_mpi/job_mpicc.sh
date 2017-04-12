#!/bin/bash
##Author:Anay Sonawane
##File Name: job_mpicc.sh
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node 28

#echo commands to stdout
set -x

module unload icc

module load pgi
module load mpi/pgi_openmpi

mpicc -acc -O  -ta=telsa:cuda8.0  cannon_mpicc.c

for var in 1 4 16 
do
	N=100
	while [ $N -le 4000 ]
	do
		mpirun -n $var ./a.out $N
		N=`expr $N + 100` 
	done

	for p in 256 512 1024 2048 4096
	do
		mpirun -n $var ./a.out $p
	done

done


