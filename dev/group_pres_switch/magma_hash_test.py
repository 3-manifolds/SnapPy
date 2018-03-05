import snappy, ntools
from sage.all import *
from sage.version import version


def hash_magma_group(G, index):
    def subgroup_hash(H):
        #ans = [G.Index(H), G.Index(G.Core(H)), H.AQInvariants(), G.Core(H).AQInvariants()]
        ans = [G.Index(H), G.Index(G.Core(H)), H.AQInvariants()]
        return [x.sage() for x in ans]

    sgs = G.LowIndexSubgroups("<1,%d>" % index)
    return sorted([subgroup_hash(H) for H in sgs])

def hash_fundamental_group_presentation(M, index):
    G = magma(M.fundamental_group())
    return hash_magma_group(G, index)
   
def test_closed():
    out = ntools.DataOutFile('/tmp/snappy-out-' + version)
    for M in snappy.OrientableClosedCensus():
        out.write( [M, hash_fundamental_group_presentation(M, 8) ] )

test_closed()
    
