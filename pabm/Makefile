# © 2013, Eric Shook, CIGI, University of Illinois at Urbana-Champaign. All rights reserved.
PROG=abm-example 
FILES=main.c sim.c mpi.c env.c agent.c sugar.c lb.c migrate.c 
COMPFLAGS=-O2 -lm -o $(PROG) -g -w 
MPIRUN=mpirun
GCC=mpicc

# Uncomment to include IPM profiling if ipm is compiled in local directory
# Available at : http://ipm-hpc.sourceforge.net/userguide.html
#IPM=-Lipm/lib -lipm

COMP=$(GCC) $(COMPFLAGS) $(FILES) $(IPM)

RUN=./$(PROG) 
MPI=$(MPIRUN) -np 4 $(RUN) 

all: 
	rm -f $(PROG)
	date
	$(COMP)

run: all 
	time $(MPI)

clean:  
	rm -f $(PROG)



