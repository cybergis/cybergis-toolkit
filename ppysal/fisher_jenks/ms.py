"""
merge sort parallelization
"""

import math


import os, sys, time
import math, random
from multiprocessing import Process, Manager


def mergeSort(a):
    len_a = len(a)
    if len_a <= 1:
        return a
    m = int(math.floor(len_a) / 2)
    left = a[0:m]
    right = a[m:]
    left = mergeSort(left)
    right = mergeSort(right)
    return merge(left, right)

def merge(left, right):
    a = []
    while len(left) > 0 or len(right) > 0:
        if len(left) > 0 and len(right) > 0:
            if left[0] <= right[0]:
                a.append(left.pop(0))
            else:
                a.append(right.pop(0))
        elif len(left) > 0:
            a.append(left.pop(0))
        elif len(right)> 0:
            a.append(right.pop(0))
    return a


def mergeSortP(sub_list):
    responses.append(mergeSort(sub_list))


def mergeP(sub_list_left, sub_list_right):
    responses.append(merge(sub_list_left, sub_list_right)) 


if __name__ == '__main__':



    l = [ 1, 2, 7]
    r = [ 3, 5, 20]

    s = l[:]
    s.extend(r[:])

    print merge(l, r)
    print s
    print mergeSort(s)

    # mp solution
    manager = Manager()
    responses = manager.list()

    max_n = 5 * 10**5
    print max_n

    a = [random.randint(0, n*100) for n in range(0, max_n)]

    t0 = time.time()
    s = mergeSort(a)
    t1 = time.time()
    print 'sequential ms: ', t1-t0
    s_p = a[:]
    t2 = time.time()
    s_p.sort()
    t3 = time.time()
    print 'python: ', t3-t2

    cores = 4

    if cores > 1:
        t4 = time.time()
        step = int( math.floor(1/ cores))
        offset = 0
        p = []
        for n in range(0, cores):
            if n < cores - 1:
                proc = Process(target=mergeSortP,
                        args=(a[n*step:(n+1)*step],))
            else:
                proc = Process(target=mergeSortP, args=( a[n*step:],))
            p.append(proc)
        for proc in p:
            proc.start()
        for proc in p:
            proc.join()
        t5 = time.time()
        print 'Final merge'
        t6 = time.time()
        p = []
        if len(responses) > 2:
            while len(responses) > 0:
                proc = Process(target=mergeP, args=(responses.pop(0),
                    responses.pop(0)))
                p.append(proc)
            for proc in p:
                proc.start()
            for proc in p:
                proc.join()
        a = merge(responses[0], responses[1])
        t7 = time.time()
        print 'mp time: ', t7-t4
        print 'final merge time: ', t7-t6




