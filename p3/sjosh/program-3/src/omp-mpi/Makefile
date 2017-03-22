CC = mpicc
DEPS = my-mpi.h
CFLAGS = -omp -O4 -Munroll
OBJ = my-mpi.o omp_mpi_floyd.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

../../bin/omp_mpi_floyd: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f *.o
	rm -f ../../bin/omp_mpi_floyd
