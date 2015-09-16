#! /usr/bin/env python

from pysal.spreg import ols, twosls_sp, error_sp_hom, error_sp_het

from mpi4py import MPI

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()

def wrapper(func, args, kwargs):
    return func(*args, **kwargs)

#Treat the search tree like a processing Queue
def processingqueue(d):
    for key, value in d.iteritems():
        if not isinstance(value, dict):
            yield key, value
        elif isinstance(value, dict):
            for k in processingqueue(value):
                yield k
        else:
            return

def chunklist(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

#Dummy Variables
y = 0
x = 0
w = 0
gwk = 0
name_y = 'str'
name_x = 'str'
name_w = 'str'
name_gwk = 'str'
name_ds = 'str'

if rank == 0:
    searchtree = {'r1':[ols.OLS, [y,x], {'w':w,'gwk':gwk,'spat_diag':True,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_gwk':name_gwk,'name_ds':name_ds}],
            'lme_lag_LT_opval':
            {'rmle_lag_LT_opval': {
                'notcombo':[twosls_sp.GM_Lag, [y,x],{'w':w,'gwk':gwk,'robust':'hac','name_y':name_y,'name_x':name_x,'name_w':name_w,'name_gwk':name_gwk,'name_ds':name_ds}],
                'het_true':[error_sp_het.GM_Combo_Het,[y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                'homo_true':[error_sp_hom.GM_Combo_Hom, [y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                'rmle_LT_opval':{'het_true':[error_sp_het.GM_Error_Het,[y,x,w,],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                    'homo_true':[error_sp_hom.GM_Error_Hom,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                'rlag_LT_opval':{'het_true':[twosls_sp.GM_Lag,[y,x],{'w':w,'robust':'white','name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                    'homo_true':[twosls_sp.GM_Lag,[y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                'else':"Robust Tests not Significant - Check Model"},
            'lme_LT_opval':
            {'het_true':[error_sp_het.GM_Error_Het,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                'het_false':[error_sp_hom.GM_Error_Hom,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
            'lag_LT_opval':
            {'het_true':[twosls_sp.GM_Lag,[y,x],{'w':w,'robust':'white','name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                'homo_true':[twosls_sp.GM_Lag, [y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
            'lme_lag_GTE_opval':
            {'het_true':[ols.OLS,[y,x],{'robust':'white','name_y':name_y,'name_x':name_x,'name_ds':name_ds}],
                'homo_true':[ols.OLS,[y,x],{'w':w,'gwk':gwk,'spat_diag':True,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_gwk':name_gwk,'name_ds':name_ds}]}
            }

#Assume that all models require the same amount of compute time.
jobs = [job for job in processingqueue(searchtree)]
njobs = len(jobs)
pjobs = njobs // (ncores - 1)  # Rank 0 is the manager.
segmentedjobs = [i for i in chunklist(jobs, pjobs)]

for i, tasks in enumerate(segmentedjobs):
    if rank == 0:
        comm.send(tasks, dest=i + 1)
    elif rank == i+1:
        tasks = comm.recv(source=0)

for r in range(ncores):
    if r != 0:
        print tasks
