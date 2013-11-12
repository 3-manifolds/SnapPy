# Initial timings, quad_double is 100 times slower than doubles.
# New timings, quad_double is 37 times slower than doubles.
# 

import snappy

def test_double():
    for M in snappy.OrientableCuspedCensus[:500]:
        snappy.Manifold(M.name())

def test_none():
    for M in snappy.OrientableCuspedCensus[:500]:
        pass

def test_quad():
    for M in snappy.OrientableCuspedCensus[:500]:
        M.high_precision()

def test_snap():
    for M in snappy.OrientableCuspedCensus[:500]:
        snappy.snap.polished_tetrahedra_shapes(M, bits_prec=230)

        
