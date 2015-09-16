#!/usr/bin/env python

"""
Automatic spatial regression using a speculative parallelization technique.
Original serial implementation by Luc Anselin.
"""

__author__ = "Jason Laura jason.laura@asu.edu, Luc Anselin luc.anselin@asu.edu"

#USAGE: mpirun -np # auto_mpi.py

from pysal.spreg import ols, twosls_sp, error_sp_hom, error_sp_het
from mpi4py import MPI
import numpy as np
import pysal as ps
import collections


#MPI tag key
#    READY = 1
#    DONE  = 2
#    EXIT = 3
#    START = 4


def wrapper(func, args, kwargs):
    return func(*args, **kwargs)


#Treat the search tree like a processing Queue
def processingqueue(d):
    parent = None
    for key, value in d.iteritems():
        if not isinstance(value, dict):
            yield key, value
        elif isinstance(value, dict):
            parent = key
            for k in processingqueue(value):
                yield key, k
        else:
            return

def chunklist(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable)\
                and not isinstance(el, basestring)\
                and not isinstance(el, np.ndarray)\
                and not isinstance(el, dict)\
                and not isinstance(el, list):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def getfromdict(datadict, maplist):
    return reduce(lambda d, k: d[k], maplist, datadict)


def populatesearchtree(searchtree, results):
    result = results[1]
    searchkey = results[0]
    getfromdict(searchtree, searchkey[:-1])[searchkey[-1]] = result
    return searchtree


def traversetree(searchtree, opvalue, combo):
    results = {'Final Model':None, 'heteroskedasticity':False,
              'spatial lag':False, 'spatial error': False,
              'regression1':None, 'regression2':None}

    r1 = searchtree['r1']
    results['regression1'] = r1
    het = r1.koenker_bassett['pvalue']
    if het < opvalue:
        hetflag = True
    else:
        hetflag = False
    results['heteroskedasticity'] = hetflag
    LMError1 = r1.lm_error[1]
    LMLag1 = r1.lm_lag[1]
    if LMError1 < opvalue and LMLag1 < opvalue:
        RLMError1 = r1.rlm_error[1]
        RLMLag1 = r1.rlm_lag[1]
        if RLMError1 < opvalue and RLMLag1 < opvalue:
            results['spatial lag']=True
            results['spatial error']=True
            if not combo:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rmle_lag_LT_opval']['notcombo']
                results['Final Model']="Spatial Lag with Spatial Error - HAC"
            elif hetflag:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rmle_lag_LT_opval']['het_true']
                results['Final Model']="Spatial Lag with Spatial Error - Heteroskedastic"
            else:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rmle_lag_LT_opval']['homo_true']
                results['Final Model']="Spatial Lag with Spatial Error - Homoskedastic"
        elif RLMError1 < opvalue:
            results['spatial error']=True
            if hetflag:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rmle_LT_opval']['het_true']
                results['Final Model']="Spatial Error - Heteroskedastic"
            else:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rmle_LT_opval']['homo_true']
                results['Final Model']="Spatial Error - Homoskedastic"
        elif RLMLag1 < opvalue:
            results['spatial lag']=True
            if hetflag:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rlag_LT_opval']['het_true']
                results['Final Model']="Spatial Lag - Heteroskedastic"
            else:
                results['regression2'] = searchtree['lme_lag_LT_opval']['rlag_LT_opval']['homo_true']
                results['Final Model']="Spatial Lag - Homoskedastic"
        else:
            results['regression2'] = None
            results['Final Model']="Robust Tests not Significant - Check Model"
    elif LMError1 < opvalue:
        results['spatial error']=True
        if hetflag:
            results['regression2'] = searchtree['lme_LT_opval']['het_true']
            results['Final Model']="Spatial Error - Heteroskedastic"
        else:
            results['regression2'] = searchtree['lme_LT_opval']['homo_true']
            results['Final Model']="Spatial Error - Homoskedastic"
    elif LMLag1 < opvalue:
        results['spatial lag']=True
        if hetflag:
            results['regression2'] = searchtree['lag_LT_opval']['het_true']
            results['Final Model']="Spatial Lag - Heteroskedastic"
        else:
            results['regression2'] = searchtree['lag_LT_opval']['homo_true']
            results['Final Model']="Spatial Lag - Homoskedastic"
    else:
        if hetflag:
            results['regression2'] = searchtree['lme_lag_GTE_opval']['het_true']
            results['Final Model'] = "No Space - Heteroskedastic"
        else:
            results['regression2'] = searchtree['lme_lag_GTE_opval']['homo_true']
            results['Final Model']="No Space - Homoskedastic"

    return results

def main(opvalue=0.01, combo=False, datafile='columbus.dbf'):

    #MPI Boilerplate
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncores = comm.Get_size()
    status = MPI.Status()

    if rank == 0:
        #Read the data
        db = ps.open(ps.examples.get_path(datafile), 'r')

        y = np.array(db.by_col('CRIME')).reshape((49,1))
        x = []
        x.append(db.by_col('INC'))
        x.append(db.by_col('HOVAL'))
        x = np.array(x).T

        #ROD Errors - the W is not pickable...generate it local to the child.
        #w = ps.weights.user.queen_from_shapefile(ps.examples.get_path('columbus.shp'))
        #w.transform = 'r'
        #gwk = ps.weights.insert_diagonal(w, diagonal=1.0)

        w = None
        gwk = None

        #Def. variables
        name_x = ['Income', 'Housing Value']
        name_y = 'Crime'
        name_w = 'Queen'
        name_ds = 'Columbus'
        name_gwk = 'Kernel from Docs'

        #Generate the search tree
        searchtree = {'r1':[ols.OLS, [y,x], {'w':w,'gwk':gwk,'spat_diag':True,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_gwk':name_gwk,'name_ds':name_ds}],
                    'lme_lag_LT_opval':
                    {'rmle_lag_LT_opval': {
                        'notcombo':[twosls_sp.GM_Lag, [y,x],{'w':w,'gwk':gwk,'robust':'hac','name_y':name_y,'name_x':name_x,'name_w':name_w,'name_gwk':name_gwk,'name_ds':name_ds}],
                        'het_true':[error_sp_het.GM_Combo_Het,[y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                        'homo_true':[error_sp_hom.GM_Combo_Hom, [y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                        'rmle_LT_opval':{'het_true':[error_sp_het.GM_Error_Het,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                                        'homo_true':[error_sp_hom.GM_Error_Hom,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                        'rlag_LT_opval':{'het_true':[twosls_sp.GM_Lag,[y,x],{'w':w,'robust':'white','name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                                        'homo_true':[twosls_sp.GM_Lag,[y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                    'else':"Robust Tests not Significant - Check Model"},
                    'lme_LT_opval':
                        {'het_true':[error_sp_het.GM_Error_Het,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                        'homo_true':[error_sp_hom.GM_Error_Hom,[y,x,w],{'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                    'lag_LT_opval':
                            {'het_true':[twosls_sp.GM_Lag,[y,x],{'w':w,'robust':'white','name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}],
                            'homo_true':[twosls_sp.GM_Lag, [y,x],{'w':w,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_ds':name_ds}]},
                    'lme_lag_GTE_opval':
                            {'het_true':[ols.OLS,[y,x],{'robust':'white','name_y':name_y,'name_x':name_x,'name_ds':name_ds}],
                            'homo_true':[ols.OLS,[y,x],{'w':w,'gwk':gwk,'spat_diag':True,'name_y':name_y,'name_x':name_x,'name_w':name_w,'name_gwk':name_gwk,'name_ds':name_ds}]}
                    }

        tasks = [list(flatten(j)) for j in processingqueue(searchtree)]
        taskidx = 0
        nworkers = ncores - 1
        closed_workers = 0

        while closed_workers < nworkers:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == 1:
                #Worker needs work
                if taskidx < len(tasks):
                    comm.send(tasks[taskidx], dest=source, tag=4)
                    taskidx += 1
                #All work is done.
                else:
                    comm.send(None, dest=source, tag=3)
            elif tag == 2:
                results = data
                searchtree = populatesearchtree(searchtree, results)
                #print results
            elif tag == 3:
                closed_workers += 1
        #All potential paths computed.  Now traverse
        results = traversetree(searchtree, opvalue=0.1, combo=False)
        #TODO: Why do the workers return as well when this is wrapped in an if rank == 0  block?
        return results

    else:
        while True:
            comm.send(None, dest=0, tag=1)
            task = comm.recv(source=0, tag=MPI.ANY_SOURCE, status=status)
            tag = status.Get_tag()
            if tag == 4:
                #Work time
                result = []
                t = task[-1]
                if isinstance(t, basestring):
                    result.append(task[:-1])
                    result.append(t)
                else:
                    result.append(task[:-1])
                    func = t[0]
                    args = t[1]
                    kwargs = t[2]
                    if 'w' in kwargs.keys():
                        w = ps.weights.user.queen_from_shapefile(ps.examples.get_path('columbus.shp'))
                        w.transform = 'r'
                        kwargs['w'] = w
                    elif len(args) == 3:
                        #Spatial error models require the w to be positional - why is this not standardized?
                        w = ps.weights.user.queen_from_shapefile(ps.examples.get_path('columbus.shp'))
                        w.transform = 'r'
                        args[2] = w
                    if 'gwk' in kwargs.keys():
                        gwk = ps.weights.user.kernelW_from_shapefile(ps.examples.get_path('columbus.shp'), 5, diagonal=True)
                        kwargs['gwk'] = gwk
                    res = wrapper(func, args, kwargs)
                    result.append(res)
                comm.send(result, dest=0, tag=2)
            elif tag == 3:
                break
        comm.send(None, dest=0, tag=3)

if __name__ == '__main__':
    import time
    t1 = time.time()
    results = main()
    if results != None:
        print results
        t2 = time.time()
        print t2 - t1

