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
>>> o.cone_fill([5.0, 6.0])
>>> o.singular_orders()
[5.0, 6.0, 2.0, 3.0, 3.0, 2.0]
>>> o.fundamental_group()
Generators:
   a,b,c
Relators:
   bAbAbAbAbA
   ccc
   aaaaBAAAcaaaaBAAAcaaaaBAAAc
   acac
   bb
   bcbcbcbcbcbc
>>> o.cone_fill(2.0, 0)
>>> o.cone_fill(3.0, 1)
>>> o.cone_fill(4.0, 2)
>>> o.cone_fill(5.0, 3)
>>> o.cone_fill(6.0, 4)
>>> o.cone_fill(2.0, 5)
>>> o.singular_orders()
[2.0, 3.0, 4.0, 5.0, 6.0, 2.0]
>>> o.volume() # doctest: +NUMERIC9
5.43335845048923

>>> o.fundamental_group()
Generators:
   a,b,c
Relators:
   cc
   bbbbbb
   bCbCbCbCbC
   abab
   acacacac
   aaaabAAAcaaaabAAAcaaaabAAAc

>>> o.cone_fill(2.1, 0)
>>> o.singular_orders()
[2.1, 3.0, 4.0, 5.0, 6.0, 2.0]
>>> o.volume() # doctest: +NUMERIC9
5.67904978263216

This give the free group of three generators in Orb, but not for us:

>>> o.fundamental_group() # doctest: +SKIP



"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

import os
from ..extensions.Orb import Orbifold
from .files import __path__ as test_files_paths
