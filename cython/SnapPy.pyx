# cython: language_level=3str
# cython: embedsignature = False

ctypedef double Real

include "SnapPy.pxi"
include "numbers/double.pyx"
include "core/basic_conversions.pyx"
include "core/basic.pyx"
include "core/triangulation.pyx"
include "core/manifold.pyx"
include "core/abelian_group.pyx"
include "core/fundamental_group.pyx"
include "core/symmetry_group.pyx"
include "core/dirichlet.pyx"
include "core/cusp_neighborhoods.pyx"
include "core/pickle.pyx"
include "core/tail.pyx"
