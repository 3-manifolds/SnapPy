from ...sage_helper import sage_method

from .cusped import *
from .closed import *

__all__ = ['verified_complex_volume_torsion']

@sage_method
def verified_complex_volume_torsion(manifold, bits_prec = None):
    
    completeness = [
        cusp_info['complete?'] for cusp_info in manifold.cusp_info() ]

    if False in completeness:
        return verified_complex_volume_closed_torsion(manifold, bits_prec)
    else:
        return verified_complex_volume_cusped_torsion(manifold, bits_prec)
