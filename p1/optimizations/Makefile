# CFLAGS = -O0 -fopenmp
# CFLAGS = -O1 -fopenmp
# CFLAGS = -O2 -fopenmp
# CFLAGS = -O3 -fopenmp
# CFLAGS = -O3 -funroll-loops -fopenmp
# CFLAGS = -O3 -funroll-loops -march=native -fopenmp
CFLAGS = -O0
LDFLAGS = -lm

clean: 
	$(RM) -f *.o lap

lap: driver.o jac1.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

lap2: driver.o jac2.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

lap2-2: driver2.o jac2.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<
