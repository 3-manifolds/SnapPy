import os, sys, re, collections
import snappy
from sage.all import PolynomialRing, QQ
from snappy.snap.find_field import *

snap_data_dir = '/pkgs/snap-pari/snap_data/'
snap_regexp = 'M:([msv]\d+)[(), ]+\n'
for kind in ['IF', 'TF', 'SF']:
    snap_regexp += kind + ':\[(.*)\]\n'

SnapField = collections.namedtuple('SnapField', ['manifold', 'inv_trace',
                                                 'trace', 'shape'])

R = PolynomialRing(QQ, 'x')
def to_R(data_str):
    return [R(poly) for poly in data_str.split(', ')]

def snap_cusp_fields():
    """
    Only gets the 1707 orientable manifolds where the file has
    complete data.
    """
    data = open(snap_data_dir + 'cusped.fields').read()
    #return  re.findall('M:([msv]\d+)', data)
    ans = []
    for match in re.findall(snap_regexp, data):
        fields = [match[0]]
        for i in range(1, 4):
            fields.append(to_R(match[i])[0])
        #fields.append(to_R(match[4]))
        ans.append(SnapField(*fields))
    return ans


def test_snappy():
    for datum in snap_cusp_fields():
        M = snappy.Manifold(datum.manifold)
        ts = M.tetrahedra_field_gens()
        ans = ts.find_field(300, 16)
        if ans is None:
            print M

def test_list_of_names(names, prec):
    for name in names:
        M = snappy.Manifold(name)
        ts = M.tetrahedra_field_gens()
        ans = ts.find_field(prec, 16, True, True)
        if ans is None:
            print M

CC = ComplexField(100)
p_z = R('x^14 + 23*x^13 - 50*x^12 + 453*x^11 - 1381*x^10 + 880*x^9 + 3113*x^8 - 8829*x^7 + 12005*x^6 - 10603*x^5 + 6456*x^4 - 2688*x^3 + 728*x^2 - 115*x + 8')
z = ExactAlgebraicNumber(p_z,
                         CC('0.57748338906663923744414542407 + 0.049579506439292237668026944791*I'))

p_elt = R('x^14 + 81*x^13 + 212*x^12 - 5040*x^11 + 25965*x^10 - 76593*x^9 + 153867*x^8 - 225060*x^7 + 246488*x^6 - 203376*x^5 + 124992*x^4 - 55488*x^3 + 16768*x^2 - 3072*x + 256')
elt = ExactAlgebraicNumber(p_elt,
                           CC('0.54070751609019766447791162339 + 0.99690919449210912885134267604*I'))

                         


bad = 's141, s164, s176, s301, s315, s456, s509, s616, s805, s902, s909,v0009, v0010'.split(',')

        
        
