# Python modules
import os
import sys
import re
import time

# Sage interaction
from ..sage_helper import _within_sage, SageNotAvailable
try:
    from sage.groups.free_group import FreeGroup
except ImportError:
    pass

from ..matrix import matrix
from .. import number

# SnapPy components
from .. import snap

SolutionType = [
    'not attempted',
    'all tetrahedra positively oriented (Note that list_interface.cpp still checks that contains_flat_tetrahedra is false)',
    'contains negatively oriented tetrahedra',
    'contains flat tetrahedra',
    'contains degenerate tetrahedra',
    'unrecognized solution type',
    'no solution found',
    'Step failed (See my_hyperbolic_structure.c for interpretation)',
    'Invalid solution (See my_hyperbolic_structure.c for interpretation)'
]
