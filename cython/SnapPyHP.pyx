# distutils: language = c++
# distutils: sources = SnapPyHP.cpp
DEF REAL_TYPE = "qd_real"
DEF HIGH_PRECISION = True
include "SnapPy.pxi"
include "core/basic.pyx"
include "core/triangulation.pyx"
include "core/manifold.pyx"
include "core/fundamental_group.pyx"
include "core/symmetry_group.pyx"
include "core/dirichlet.pyx"
include "core/cusp_neighborhoods.pyx"
include "core/tail.pyx"
