CC = mpicc
DEPS = my-mpi.h
CFLAGS = -acc -ta=tesla,cuda8.0
OBJ = my-mpi.o acc_mpi_floyd.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

../../bin/acc_mpi_floyd: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f *.o
	rm -f ../../bin/acc_mpi_floyd
