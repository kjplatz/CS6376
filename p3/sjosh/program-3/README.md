# Execution of Floyd's algorithm

### Serial execution

To compile:
<code>
    gcc -o bin/serial_floyd src/serial_floyd.c
</code>

### Execution with OpenMP

To compile:
<code>
    pgcc -omp -o bin/omp_floyd src/omp_floyd.c
</code>

### Execution with OpenACC

To compile:
<code>
    pgcc -acc -o bin/acc_floyd src/acc_floyd.c -ta=tesla,cuda8.0
</code>

### Execution with MPI

To compile and build MPI binaries, we need to load the module <code>mpi/pgi_openmpi</code>.
To load the module:
<code>
    module load pgi
    module load mpi/pgi_openmpi
</code>

To compile:
<code>
    cd src/mpi
    make clean
    make
</code>

### Execution with MPI and OpenMP

To compile:
<code>
    cd src/omp_mpi
    make clean
    make
</code>

### Execution with MPI and OpenACC

To compile:
<code>
    cd src/acc_mpi
    make clean
    make
</code>
