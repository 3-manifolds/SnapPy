# distutils: language = c++
# cython: language_level=3str
# cython: embedsignature = False

ctypedef qd_real Real

include "SnapPy.pxi"
include "numbers/qd.pyx"
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
