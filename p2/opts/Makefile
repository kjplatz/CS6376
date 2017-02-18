CFLAGS = -O3 -funroll-loops -fopenmp 
LDFLAGS = -lm

lap1: driver.o jac1.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o lap1
