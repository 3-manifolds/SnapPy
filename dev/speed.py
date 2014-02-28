# Initial timings, quad_double is 100 times slower than doubles.
# New timings, quad_double is 37 times slower than doubles.
# 

import snappy

def test_double():
    for M in snappy.OrientableCuspedCensus[:5000]:
        snappy.Manifold(M.name())

def test_none():
    for M in snappy.OrientableCuspedCensus[:5000]:
        pass

def test_quad():
    """"
    50 / 1:40 
    """
    for M in snappy.OrientableCuspedCensus[:5000]:
        M.high_precision() #.volume()
        

def test_snap():
    for M in snappy.OrientableCuspedCensus[:500]:
        snappy.snap.polished_tetrahedra_shapes(M, bits_prec=230)

        
test_quad()
