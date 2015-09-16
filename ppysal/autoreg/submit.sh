#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=2:ppn=8
#PBS -N Speculative Regression
#PBS -o out
#PBS -e err
#PBS -A Jay
#PBS -l walltime=00:20:00
use openmpi-1.6.4
cd $PBS_O_WORKDIR
/packages/openmpi-1.6.4/bin/mpiexec -n 16 python autoreg.py ~/pysal_examples/columbus.dbf pbsout.json HOVAL INC CRIME
