# cython: language_level=3str
# cython: embedsignature = False

ctypedef double Real

include "Orb.pxi"
include "core/basic_conversions.pyx"
include "core/basic.pyx"
include "core/orbifold.pyx"
