from snappy import Manifold
import tarfile, os
import snappy.SnapPy
# Make the paths point in to the current directory
snappy.SnapPy.manifold_path = manifold_path = './'
snappy.SnapPy.closed_census_directory = os.path.join(manifold_path,
                                                     'ClosedCensusData')
snappy.SnapPy.link_directory = os.path.join(manifold_path, 'ChristyLinks')
snappy.SnapPy.link_archive = link_archive = os.path.join(manifold_path, 'ChristyLinks.tgz')
snappy.SnapPy.census_knot_archive = census_knot_archive = os.path.join(manifold_path, 'CensusKnots.tgz')
snappy.SnapPy.Census_Morwen8 = tarfile.open(os.path.join(manifold_path, 'morwen8.tgz'), 'r:*')
snappy.SnapPy.Christy_links = tarfile.open(link_archive, 'r:*')
snappy.SnapPy.Census_Knots = tarfile.open(census_knot_archive, 'r:*')

from snappy.SnapPy import ObsOrientableCuspedCensus, ObsNonorientableCuspedCensus, ObsLinkExteriors, ObsCensusKnots, ObsOrientableClosedCensus, ObsNonorientableClosedCensus

id = [[1,0],[0,1]]
minus_id = [[-1,0],[0,-1]]
plusorminus_id = [id, minus_id]

def is_identity(iso):
    if not iso.extends_to_link():
        return False
    if not iso.cusp_images() == range(iso.num_cusps()):
        return False
    for matrix in iso.cusp_maps():
        if matrix.tolist() not in plusorminus_id:
            return False
    return True
    
def compare(M, N):
    isos = M.isomorphisms_to(N)
    if True in [is_identity(iso) for iso in isos]:
        return True
    print '%s != %s'%(M.name(), N.name())
    return False

def test_cusped():
    print('testing cusped census')
    for M in ObsOrientableCuspedCensus():
        N = Manifold(M.name())
        result = compare(M,N)
    print('done')

def test_links():
    print('testing links')
    for n in range(1, 6):
        for M in ObsLinkExteriors(n):
            N = Manifold(M.name())
            compare(M,N)
    print('done')

def test_census_knots():
    print('testing census_knots')
    for M in ObsCensusKnots():
        N = Manifold(M.name())
        result = compare(M,N)
    print('done')

if __name__ == '__main__':
    pass
    test_cusped()
    test_links()
    test_census_knots()
