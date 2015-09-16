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
import sys
import cPickle
import os
import json

import generateW
import argumentparser

#Hack to get the KDTree used in GWK to pickle.
# KDTree uses nested classes that the pickler can not find.
from scipy.spatial import kdtree
kdtree.node = kdtree.KDTree.node
kdtree.leafnode = kdtree.KDTree.leafnode
kdtree.innernode = kdtree.KDTree.innernode

#MPI tag key
#    READY = 1
#    DONE  = 2
#    EXIT = 3
#    START = 4

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()
status = MPI.Status()

class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.generic):
            return obj.item()
        return json.JSONEncoder.default(self, obj)

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


def loadW(datafile, adjacency):
    t1 = time.time()
    datafile = datafile.replace('.dbf', '.shp')
    basefile = datafile.split('.')[0]

    regenerateW = False

    if adjacency == 'rook':
        baseW = basefile + '_WR.pkl'
    else:
        baseW = basefile + '_WQ.pkl'

    if os.path.exists(baseW):
        with open(baseW, 'rb') as f:
            w = cPickle.load(f)
        if not w.adjacency == adjacency:
            print w.adjacency, adjacency
            regenerateW = True
    else:
        regenerateW = True

    baseGWK = basefile + '_GWK.pkl'
    if not os.path.exists(baseGWK):
        regenerateW = True
    else:
        with open(baseGWK, 'rb') as f:
            gwk = cPickle.load(f)

    if regenerateW:
        wf, gwkf = generateW.main(datafile, adjacency)
        with open(wf, 'rb') as f:
            w = cPickle.load(f)
        with open(gwkf, 'rb') as f:
            gwk = cPickle.load(f)
    t2 = time.time()
    print "W: ",t2 - t1
    return w, gwk

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

def main(dep, indep, opvalue=0.01, combo=False, datafile='columbus.dbf', adjacency='QUEEN'):

    if rank == 0:
        #Read the data
        db = ps.open(datafile, 'r')
    
        y = np.array(db.by_col(dep)).reshape(-1, 1)

        #y = np.array(db.by_col('CRIME')).reshape((49,1))
        x = []
        for i in indep:
            x.append(db.by_col(i))
        x = np.array(x).T

        #ROD Errors - the W is not pickable...generate it local to the child.
        #Setup the W
        datafile = datafile.replace('.dbf', '.shp')
        w, gwk = loadW(datafile, adjacency)

        #Def. variables
        name_x = indep
        name_y = dep
        name_w = adjacency
        name_ds = datafile
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
                    res = wrapper(func, args, kwargs)
                    result.append(res)
                comm.send(result, dest=0, tag=2)
            elif tag == 3:
                break
        comm.send(None, dest=0, tag=3)

if __name__ == '__main__':
    import time
    t1 = time.time()
    args = argumentparser.argparser()
    y = args['dependent_variable']
    x = args['independent_variable(s)']
    adjacency = args['adjacency']
    datafile = args['input']
    outjson = args['output']
    
    if datafile.split('.')[-1] == 'shp':
        datafile = datafile.replace('.shp', '.dbf')
    results = main(y, x, datafile=datafile, adjacency=adjacency)

    if results != None:
        #Setup output files
        try:
            baseout = outjson.split('.')[0]
            outsummary = baseout + '_summary.txt'
            outcsv = baseout + '_vectors.csv'
        except:
            outsummary = outjson + '_summary.txt'
            outcsv = outjson + '_vectors.csv'
        print outsummary
        print outcsv
        #Unpack the class attributes into a dict
        r1 = vars(results['regression1'])
        r1_summary = r1['summary']
        r2 = vars(results['regression2'])
        r2_summary = r2['summary']

        #Write the summary information
        with open(outsummary, 'w') as outtxt:
            outtxt.write(r1_summary)
            outtxt.write('\n')
            outtxt.write(r2_summary)
            
        results['regression1'] = r1
        results['regression2'] = r2
        reg1 = r1['u'][:,0]
        reg2 = r2['u'][:,0]
        reg1class = ps.esda.mapclassify.Quantiles(reg1, k=4).yb
        reg2class = ps.esda.mapclassify.Quantiles(reg2, k=4).yb
       
        print r1['u'].shape
        print r2['u'].shape
        print reg1class

        outvect = np.column_stack((reg1,reg1class,reg2, reg2class))
        header = 'r1_residuals r1_residuals_bins r2_residuals r2_residuals_bins'

        np.savetxt(outcsv, outvect, delimiter=',', fmt='%10.4f', header=header)

        jsonstring = json.dumps(results, cls=NumpyAwareJSONEncoder, indent=4)
        try:
            os.remove(outjson)
        except OSError:
            pass
        with open(outjson, 'w') as f:
            f.write(jsonstring)
        with open(outjson + '_done', 'w') as f:
	    f.write('done');
        t2 = time.time()
