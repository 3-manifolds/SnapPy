import snappy
import multiprocessing
import time

n = 1000

def complex_volume(M):
    return M.high_precision().complex_volume()

def normal_surfaces(M):
    return len(M.normal_surfaces())

def test_one_core():
    return [complex_volume(M)
            for M in snappy.OrientableCuspedCensus(tet=9)[:n]]

def test_many_cores(cores=8):
    pool = multiprocessing.Pool(8)
    manifolds = list(snappy.OrientableCuspedCensus(tet=9)[:n])
    return pool.map(complex_volume, manifolds)
    
if __name__ == '__main__':
    start = time.time()
    one_ans = test_one_core()
    elapsed = time.time() - start
    print('One core took %.1f seconds' % elapsed)

    start = time.time()
    quick_ans = test_many_cores()
    elapsed = time.time() - start
    print('Eight cores took %.1f seconds' % elapsed)

    print('Max difference was %f' % max(abs(a - b) for a, b in zip(one_ans, quick_ans)))
    


