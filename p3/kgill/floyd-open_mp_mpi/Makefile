CC=mpicc
C_FLAGS=-omp -O4 -Munroll -mcmodel=medium
OUTPUT=floyd_mpi.out

all: floyd.c
	$(CC) $(C_FLAGS) $^ -o $(OUTPUT)

clean:
	@rm -f $(OUTPUT)
	@rm -f slurm*
