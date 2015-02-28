import cPickle as pickle
import snappy
import multiprocessing


#def reduce(self):
#    return (snappy.Manifold, (self._to_string(),))
#
#snappy.Manifold.__reduce__ = reduce


manifolds = [M.high_precision() for M in snappy.OrientableCuspedCensus[:100]]

def volume(M):
    return M.complex_volume()

if __name__ == '__main__':
    p = multiprocessing.Pool(4)
    print p.map(volume, manifolds)




