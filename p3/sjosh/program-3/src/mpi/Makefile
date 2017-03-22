CC = mpicc
DEPS = my-mpi.h
OBJ = my-mpi.o mpi_floyd.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

../../bin/mpi_floyd: $(OBJ)
	$(CC) -o $@ $^

clean:
	rm -f *.o
	rm -f ../../bin/mpi_floyd
