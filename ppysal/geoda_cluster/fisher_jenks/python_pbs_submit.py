#!/usr/bin/python

from popen2 import popen2
import time
import sys
import os
import subprocess

import numpy as np
import pysal as ps
from fj_vect import fisher_jenks as vFisher

rho = [0.9, 0.7, 0.5, 0.3, 0.1, 0.0, -0.1, -0.3, -0.5, -0.7, -0.9]
ns = [75, 100, 125, 150]
ks = [5,7,9]
ss = [0.05, 0.10, 0.15, 0.20, 0.25]
realizations = int(sys.argv[1])

#Loop through the parameter combinatio(k, y, pivots):
def get_tss(k, y, pivots):
    """
    Total sum of squares around class means

    Returns sum of squares over all class means
    """
    tss = 0
    y = np.sort(y)
    classes = np.arange(k)
    starts = [0]
    for i in range(len(pivots)):
        starts.append(pivots[i])
        pivots.append(len(y))
        breaks = zip(starts, pivots)
    for b in breaks:
        ymean = y[b[0]:b[1]].mean()
        css = y[b[0]:b[1]] - ymean
        css *= css
        tss += np.sum(css)
    return tss
#Launch a job set for each combination
for i, n in enumerate(ns):
    for r, p in enumerate(rho):
        for j, k in enumerate(ks):
            ta = time.time()
            #Compute all the samples
            output, input = popen2('qsub')
            jobname = 'FJ_{}_{}_{}'.format(n, p, k)
            if n <= 50:
                processors = 'nodes=8:ppn=8'
                walltime = '06:00:00'
                memory = 'pmem=1gb'
            elif n <= 100:
                walltime = '10:00:00'
                processors = 'nodes=32:ppn=2'
                memory = 'pmem=4gb'
            elif n == 100:
                walltime = '20:00:00'
                processors = 'nodes=32:ppn=2'
                memory = 'mem=3gb'
            elif n == 125:
                walltime = '280:00:00'
                processors = 'nodes=32:ppn=2'
                memory = 'pmem=8gb'
            elif n == 150:
                walltime = '00:30:00'
                processors = 'nodes=64:ppn=1'
                memory = 'pmem=16gb'
            elif n == 175:
                walltime = '100:00:00'
                processors = 'nodes=64:ppn=1'
                memory = 'pmem=8gb'
            elif n == 200:
                walltime = '100:00:00'
                processors = 'nodes=64:ppn=1'
                memory = 'pmem=8gb'
            # --hostfile $PBS_NODEFILE
            #command = '/packages/openmpi-1.6.4/bin/mpiexec -n 64 --hostfile nodelist.lis python fj_mpi_integrated.py {} {} {} {} 1>>output 2>>error'.format(realizations, n, k, p )
            command = '/packages/openmpi-1.6.4/bin/mpiexec -n 64 python fj_mpi_integrated.py {} {} {} {} 1>>output 2>>error'.format(realizations, n, k, p )
            
            job_string = '''#!/bin/bash
            #PBS -S /bin/bash
            #PBS -N {}
            #PBS -l walltime={}
            #PBS -l {}
            #PBS -l {}
            #PBS -o $JOB_ID.out
            #PBS -e $JOB_ID.err
            #PBS -A Jay_Testing
            #PBS -m ae
            #PBS -M jlaura@asu.edu
            use openmpi-1.6.4
            cd $PBS_O_WORKDIR
            echo $PBS_NODEFILE; cat $PBS_NODEFILE
            {}
            '''.format(jobname, walltime, processors, memory, command)

            input.write(job_string)
            input.close()

            #Prints the PBS job submission string
            print job_string
            job_id = output.read()
            #Little sleep to flush
            time.sleep(0.1)
    #At this indent level we compute all rho and all classes for some number of samples in parallel.
    processing = True
    while processing:
        time.sleep(30)
        output = subprocess.Popen(['qselect', '-u', 'jlaura', '-s', 'R'], stdout=subprocess.PIPE).communicate()[0]
        running_count = output.count('\n')
        output = subprocess.Popen(['qselect', '-u', 'jlaura', '-s', 'Q'], stdout=subprocess.PIPE).communicate()[0]
        running_count += output.count('\n')
        if running_count == 0:
            processing = False

