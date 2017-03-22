CC=pgcc
C_FLAGS=-O4 -acc -ta=tesla,8.0 -Munroll -mcmodel=medium
# C_FLAGS+=-Minfo=accel
OUTPUT=floyd_acc_8192_.out

all: floyd.c
	$(CC) $(C_FLAGS) $^ -o $(OUTPUT)

clean:
	@rm -f $(OUTPUT)
	@rm -f slurm*
