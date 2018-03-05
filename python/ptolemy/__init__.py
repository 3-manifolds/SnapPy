### ptolemy module implements methods to find and solve the Ptolemy variety
### 
### 2012 - "Matthias Goerner" <enischte@gmail.com> 

from .coordinates import PtolemyCoordinates, Flattenings, CrossRatios
from .processMagmaFile import solutions_from_magma, solutions_from_magma_file
from .processFileDispatch import parse_solutions, parse_solutions_from_file
from .ptolemyGeneralizedObstructionClass import PtolemyGeneralizedObstructionClass

import os as _os

_env = _os.environ.get('PTOLEMY_DATA_URL')

if _env:
    DATA_URL = _env
else:
    DATA_URL = "http://ptolemy.unhyperbolic.org/"


