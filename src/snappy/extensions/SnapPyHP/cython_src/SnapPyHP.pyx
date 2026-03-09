# distutils: language = c++
# cython: language_level=3str
# cython: embedsignature = False

ctypedef qd_real Real

include "../../SnapPy/cython_src/SnapPy.pxi"
include "../../SnapPy/cython_src/precision/qd/number.pyx"
include "../../SnapPy/cython_src/core/basic.pyx"
include "../../SnapPy/cython_src/core/triangulation.pyx"
include "../../SnapPy/cython_src/core/manifold.pyx"
include "../../SnapPy/cython_src/core/abelian_group.pyx"
include "../../SnapPy/cython_src/core/fundamental_group.pyx"
include "../../SnapPy/cython_src/core/symmetry_group.pyx"
include "../../SnapPy/cython_src/core/dirichlet.pyx"
include "../../SnapPy/cython_src/core/cusp_neighborhoods.pyx"
include "../../SnapPy/cython_src/core/pickle.pyx"
include "../../SnapPy/cython_src/core/tail.pyx"
