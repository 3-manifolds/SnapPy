"""
With new code, can't find trace field of s075, s230 [...]

With old (2015/7/19) code, couldn't do:

m077
m125
m136
m203
m267
m306
m307
m319
m364
m403
s019
s045
s066
s074
s075
s079
[...]

"""

import snappy
from sage.all import vector

def find_trace_field(M, max_prec=500,  optimize_field=False):
    """
    Starts with 100 bits of precision and degree 10 and then doubles
    both successively until it succeeds or max_prec is bit.  The ratio
    of bits/degree is roughly the one recommended in [CGHN]
    """
    traces = M.trace_field_gens()
    prec, deg = 100, 10
    ans = None
    while ans is None and prec <= max_prec:
        ans = traces.find_field(prec, deg, optimize_field, verbosity=False)
        prec, deg = 2 * prec, 2 * deg
    #if ans is None:
    #    raise ValueError('Could not compute trace field')
    return ans

def test_trace_fields():
    for M in snappy.OrientableCuspedCensus:
        ans = find_trace_field(M)
        if ans is None:
            print(M.name())

def test_shapes():
    for M in snappy.OrientableCuspedCensus:
        shapes = M.tetrahedra_field_gens()
        ans = shapes.find_field(300, 20, verbosity=False)
        if ans is None:
            print(M.name())

test_trace_fields()
#M = snappy.Manifold('m320')
#print(find_trace_field(M))
