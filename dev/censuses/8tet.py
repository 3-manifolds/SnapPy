import snappy, os, sys, re, gzip
from nsagetools import homological_longitude

def raw_data():
    def manifold(line):
        k_index, t_index = re.match('\$k8_{(\d+)}\$ \& (\d+) \&', line).groups()
        census_name = 't'+t_index.rjust(5, '0')
        M = snappy.Manifold(census_name)
        M.set_name('K8_'+k_index + '\t' + M.name())
        return M
    return [manifold(line) for line in gzip.open('table-8tet.gz')]

def meridian(M):
    T = M.without_hyperbolic_structure()
    for slope in [ (1,0), (0,1), (1,1), (1,-1)]:
        T.dehn_fill(slope)
        for i in range(15):
            if T.fundamental_group().num_generators() == 0:
                return slope
            T.randomize()

def one_line(M):
    m, l = meridian(M), homological_longitude(M)
    if m[0]*l[1] - m[1]*l[0] < 0:
        l = -l
    return M.name() + '\t%s\n' % repr([m, tuple(l)])
    
def main():
    file = gzip.open('knots_8_tet.gz', 'w')
    file.write(''.join(one_line(M) for M in raw_data()))

def test():
    added = [M for M in snappy.CensusKnots if M.name().startswith('K8')]
    assert len(added) == 299
    for M in snappy.CensusKnots:
        T = M.without_hyperbolic_structure()
        T.dehn_fill( (1,0) )
        found = False
        for i in range(100):
            if T.fundamental_group().num_generators() == 0:
                found = True
                break
            else:
                T.randomize()
        assert found
        T.dehn_fill( (0, 1) )
        assert T.homology().betti_number() == 1

test()
