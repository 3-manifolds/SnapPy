"""

>>> o = Orbifold(os.path.join(test_files_paths[0], '6_5^2.7.orb'))
>>> o.volume() # doctest: +NUMERIC9
0.117838420347115

>>> o = Orbifold(os.path.join(test_files_paths[0], '1_1^4.84.orb'))
>>> o.fundamental_group(False)
Generators:
   a,b,c,d,e,f,g
Relators:
   FF
   DDD
   DfDfDf
   Gbee
   gg
   cbDB
   cge
   AA
   ACAC
   Aef
>>> o.fundamental_group(True)
Generators:
   a,b,c
Relators:
   cc
   bbb
   bCbCbC
   abab
   acac
   aaaabAAAcaaaabAAAc

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

import os
from ..extensions.Orb import Orbifold
from .files import __path__ as test_files_paths
