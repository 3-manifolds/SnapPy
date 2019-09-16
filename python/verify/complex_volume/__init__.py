from ...sage_helper import sage_method

from .cusped import *
from .closed import *

@sage_method
def complex_volume_torsion(manifold, bits_prec = None):
    
    completeness = [
        cusp_info['complete?'] for cusp_info in manifold.cusp_info() ]

    if False in completeness:
        return complex_volume_closed_torsion(manifold, bits_prec)
    else:
        return complex_volume_cusped_torsion(manifold, bits_prec)
