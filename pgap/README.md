# PGAP: A scalable Parallel Genetic Algorithm (PGA) solver for the Generalized Assignment Problem (GAP)
Author: **Yan Y. Liu <yanliu@illinois.edu>**

## Introduction
This is a PGA application that uses MPI non-blocking functions for migrating solutions among PGA demes. The implementation of PGA is based on the island model. Although this solver is written for solving the Generalized Assignment Problem (GAP), it can be easily modified to solve other combinatorial optimization problems. The scalability of this application has been tested on BlueWaters supercomputer at the National Center for Supercomputing Applications (NCSA), Stampede supercomputer at the Texas Advanced Computing Center (TACC), and Trestles supercomputer at the San Diego Supercomputer Center (SDSC). This code can scale to 262K processor cores on BlueWaters with marginal communcation cost. For more details of the underlying techniques used to achieve desirable scalability of this code, please refer to the following publication:

"Yan Y. Liu, Shaowen Wang, A scalable parallel genetic algorithm for the Generalized Assignment Problem, Parallel Computing, http://dx.doi.org/10.1016/j.parco.2014.04.008, 2014."

## Dependency
  - MPI: OpenMPI/MPICH/MVAPICH2/Cray MPI
  - SPRNG: a parallel random number generator (http://www.sprng.org/). Environment variable ```$SPRNG_HOME``` must be set to the path to sprng. sprng must be compiled with MPI support. Need to use SPRNG2.0 (http://www.sprng.org/Version2.0/index.html)

## Compile
make clean; make

Please read Makefile for the changes you need to make first. You can make sequential version, parallel version, and Intel Phi MIC version of the code by selecting different compiling options. The default compiling option is parallel and non-blocking.

## Data
Input data must conform OR-LIB GAP format. GAP problem instance data can be obtained from OR-LIB:
  - Small instances: http://people.brunel.ac.uk/~mastjjb/jeb/orlib/gapinfo.html
  - Large instances: http://www.al.cm.is.nagoya-u.ac.jp/~yagiura/gap/

## Run
To run the program, please use the mpirun command on your platform, e.g., aprun on Blue Waters; ibrun on Stampede; ```mpirun_rsh``` on Trestles; or mpirun using openmpi/mpich/mvapichi2 .

A sample run that uses 4 processes and applies binary selection, unfitness replacement, constraint-based populication initialization, population size: 200, execution time: 300 seconds, export interval: every 100 iterations, migration rate: 2, import interval: every 50 iterations, sending parallelism: 2 can be launched with the following command:

```
mpirun -np 4 ga-async -s b -r u -i c -p 200 -t t 300 -e -l jobname -f /path/to/dataset/gap/y/e801600 -M 100 2 50 -P 2
```

Depending on the configuration of underlying cluster interconnect, non-blocking communcation may or may not be stable. If you run this code using mvapich2, you are advised to increase the MV2_HYBRID_ENABLE_THRESHOLD by inserting the following statement into the command line:

```
MV2_HYBRID_ENABLE_THRESHOLD=2048
```

## Results
Every process outputs its solutions and logging info to a file. Final best solution is output the end of process output file.

**Contact: Yan Y. Liu <yanliu@illinois.edu>**
