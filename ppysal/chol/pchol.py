import numpy as np
import multiprocessing as mp
import ctypes

def main(arr, lock):
    """
    arr - an n x n matrix
    lock - a semaphore to synchronize computation
    """
    print arr

if __name__ == '__main__':
    n = 10
    carr = mp.Array(ctypes.c_double, n*n, lock=True)
    lock = carr.get_lock()

    a = np.random.randint(0,10, size=(n,n)).astype(np.float64)
    a = np.dot(a, a.T)
    with lock:
        arr = np.frombuffer(carr.get_obj(), dtype=np.float).reshape(n,n)
        arr[:] = a
    print np.linalg.cholesky(arr)
    main(arr, lock)
