#!/bin/sh

#Author: ANAY Sonawane
#File Name: job_openmp.sh

echo Executing Program 1
pgcc -mp floyd_openmp.c

file="pgm1_out.out"
MAX_NUM_THREADS=28
NUM_THREADS=0
program_time=0

parse_output(){
var4=1.2
while read LINE
do
	var4="$(echo $LINE | cut -d' ' -f 4)"
	echo $1 " " $2 " " $var4 >> time.out
done < $file
}


while [ $NUM_THREADS -lt $MAX_NUM_THREADS ]
do

NUM_THREADS=`expr $NUM_THREADS + 1`
echo $NUM_THREADS

for var in 100 1000 2000 3000 4000 8000 
do

export OMP_NUM_THREADS=$NUM_THREADS
echo $var | ./a.out | tail -1 > $file
parse_output $NUM_THREADS $var $program_time

done

done
