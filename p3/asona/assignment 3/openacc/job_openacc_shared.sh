#!/bin/bash
##Author: Anay Sonawane
##File Name: job_openacc_shared.sh
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH --ntasks-per-node 7
#SBATCH --gres=gpu:1
#SBATCH -t 5:00:00

#echo commands to stdout
set -x

pgcc -acc -Minfo=accel -ta=tesla:cuda8.0 floyd_openacc.c

echo 2000 | ./a.out
