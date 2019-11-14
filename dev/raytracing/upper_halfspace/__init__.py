"""
Classes and functions for the upper half space model **H**\ :sup:`3` of hyperbolic space
using interval arithmetics.

"""

try:
    from extendedMatrix import *
    from idealPoint import *
    from projectivePoint import *
    from finitePoint import *
    from reflections import *
except:
    from .extendedMatrix import *
    from .idealPoint import *
    from .projectivePoint import *
    from .finitePoint import *
    from .reflections import *
    
