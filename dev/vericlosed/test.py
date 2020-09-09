import branchedDoubleCover
import veriClosed
from veriClosed import *
from veriClosed.verifyHyperbolicStructureEngine import *
from veriClosed.testing.cocycleTester import *

from veriClosed.testing import __path__ as testPaths
testPath = testPaths[0]

from snappy import Triangulation
from snappy.sage_helper import doctest_modules

import os
import sys

def testCocycles(isoSig):
    """
    Testing.

        >>> testCocycles("kLLLvQQkccfgighjijjlnannwnashp")
        True
        >>> testCocycles("qLLvLAzMAQAkcdjifjhlnkmlnnopphsksskcimafjwovqw")
        True

    """    

        
    T = Triangulation(isoSig, remove_finite_vertices = False)
    hyperbolic_structure = compute_verified_hyperbolic_structure(T, source = 'new', bits_prec = 106)
    engine = VerifyHyperbolicStructureEngine(hyperbolic_structure)
    tester = CocycleTester(engine)
    tester.test()

    return True

def testVerifyFromVertexGramMatrixFile():
    """
    Testing.
      
        >>> f = os.path.join(testPath, "m004_1_2.tri")
        >>> T = Triangulation(f, remove_finite_vertices = False)
        >>> bool(compute_verified_hyperbolic_structure(T, source = f + '.vgm'))
        True
    """

def _doctest():
    doctest_modules(
        [ sys.modules[__name__],
          veriClosed.computeVerifiedHyperbolicStructure,
          branchedDoubleCover])

if __name__ == '__main__':
    _doctest()
    
