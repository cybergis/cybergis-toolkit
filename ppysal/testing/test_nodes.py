import subprocess
from popen2 import popen2
import time

pipe = subprocess.Popen(['pbsnodes', '-a'],stdout=subprocess.PIPE)
raw_nodes =  pipe.communicate()[0].split('\n')
nodes = []
for i,l in enumerate(raw_nodes):
    if 'magic-' in l and not 'uname' in l and not 'down' in raw_nodes[i+1]:
        nodes.append(l)
    elif 'magic-' in l and 'down' in raw_nodes[i+1]:
        print "Node {} is down and will not be tested.".format(l)

with open('test1.nodelist', 'w') as f:
    for n in nodes:
        f.write(n + '\n')
nnodes = len(nodes)
print "pbsnodes -a reports a total of {} nodes as up.".format(nnodes)

command = '/packages/openmpi-1.6.4/bin/mpiexec -n 8 python singletest.py 1>>test1.out 2>>test1.err'

for n in nodes:
    print n
    if n == 'magic-4-4-4.magic':
       print "Skipping magic-4-4-4.magic"
       continue
    
    output, input = popen2('qsub')
    singlejob='''#!/bin/bash
    #PBS -S /bin/bash
    #PBS -N Testing Node {}
    #PBS -l nodes={}:ppn=8
    #PBS -l walltime=00:05:00
    #PBS -A Testing
    #PBS -o $JOB_ID.out
    #PBS -e $JOB_ID.err
    #PBS -M jlaura@asu.edu
    cd $PBS_O_WORKDIR
    use openmpi-1.6.4
    {}
    '''.format(n,n, command)
    input.write(singlejob)
    input.close()
    time.sleep(0.1)
    
