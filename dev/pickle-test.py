import cPickle as pickle
import snappy
import multiprocessing
import 

def complex_volume(M):
    return M.high_precision().complex_volume()

def normal_surfaces(M):
    return len(M.normal_surfaces())

def test_one_core():
    i = 0
    for M in snappy.OrientableCuspedCensus(tet=9)[:1000]:
        complex_volume(M)
	i += 1
    print i

def test_many_cores(cores=8):
    p = multiprocessing.Pool(8)
    manifolds = list(snappy.OrientableCuspedCensus(tet=9)[:1000])
    print len(manifolds)
    print len(p.map(complex_volume, manifolds))
    
if __name__ == '__main__':
    #test_many_cores()
    test_one_core()



