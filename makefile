.PHONY: clean

executable : main.o tools.o my_mpi_aa.o serial_aa.o
	mpicc -fopenmp main.o tools.o my_mpi_aa.o serial_aa.o -o executable -lm

main.o : main.c my_aa.h
	mpicc -Wall -c main.c

tools.o : tools.c tools.h
	mpicc -Wall -c tools.c

my_mpi_aa.o : my_mpi_aa.c tools.h
	mpicc -fopenmp -Wall -c my_mpi_aa.c

serial_aa.o : serial_aa.c
	mpicc -fopenmp -Wall -c serial_aa.c

clean:
	rm *.o executable
