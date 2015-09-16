#!/bin/bash

#PBS -S /bin/bash
#PBS -N Generated Structured
#PBS -l walltime=20:00:00
#PBS -l nodes=22:ppn=1
#PBS -l pmem=8gb
#PBS -o $JOB_ID.out
#PBS -e $JOB_ID.err
#PBS -A Jay_Testing
use openmpi-1.6.4
cd $PBS_O_WORKDIR
/packages/openmpi-1.6.4/bin/mpiexec -n 22 --hostfile $PBS_NODEFILE python gen_structured.py 1>gen_output 2>generror

