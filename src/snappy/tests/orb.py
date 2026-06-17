"""

>>> from snappy.extensions import SnapPy
>>> from snappy.extensions.SnapPy import Triangulation
>>> from snappy.extensions.SnapPy import Orbifold
>>> SnapPy._orb_set_use_orb_conventions(True)
>>> O = Orbifold(os.path.join(test_files_paths[0], '6_5^2.7.orb'), remove_finite_vertices = False)
>>> O.solution_type()
'partially flat tetrahedra'
>>> O.volume() # doctest: +NUMERIC9
0.117838420347115
>>> O
6_5^2.7(2)(2)(2)
>>> O.cusp_info()
[]
>>> O._orb_vertex_info()
[Vertex 0 : spherical link with cone points of order [2.0, 2.0, 2.0] and orbifold Euler char 0.5,
 Vertex 1 : spherical link with cone points of order [2.0, 2.0, 2.0] and orbifold Euler char 0.5]

>>> O = Triangulation(os.path.join(test_files_paths[0], '1_1^4.84.orb'), remove_finite_vertices = False)
>>> O.fundamental_group(False)
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
>>> O.fundamental_group(True)
Generators:
   a,b,c
Relators:
   cc
   bbb
   bCbCbC
   abab
   acac
   aaaabAAAcaaaabAAAc
>>> O._orb_cone_fill([5.0, 6.0])
>>> O
1_1^4.84(5)(6)(2)(3)(3)(2)
>>> O.cusp_info()
[]
>>> O._orb_vertex_info()
[Vertex 0 : spherical link with cone points of order [5.0, 6.0, 2.0] and orbifold Euler char -0.133333,
 Vertex 1 : spherical link with cone points of order [6.0, 2.0, 3.0] and orbifold Euler char -5.55112e-17,
 Vertex 2 : spherical link with cone points of order [2.0, 3.0, 2.0] and orbifold Euler char 0.333333,
 Vertex 3 : spherical link with cone points of order [5.0, 3.0, 3.0] and orbifold Euler char -0.133333]
>>> O._orb_singular_edge_info('singular_order')
[5.0, 6.0, 2.0, 3.0, 3.0, 2.0]
>>> O.fundamental_group()
Generators:
   a,b,c
Relators:
   bAbAbAbAbA
   ccc
   aaaaBAAAcaaaaBAAAcaaaaBAAAc
   acac
   bb
   bcbcbcbcbcbc
>>> O._orb_cone_fill(2.0, 0)
>>> O._orb_cone_fill(3.0, 1)
>>> O._orb_cone_fill(4.0, 2)
>>> O._orb_cone_fill(5.0, 3)
>>> O._orb_cone_fill(6.0, 4)
>>> O._orb_cone_fill(2.0, 5)
>>> O._orb_singular_edge_info('singular_order')
[2.0, 3.0, 4.0, 5.0, 6.0, 2.0]
>>> O.fundamental_group()
Generators:
   a,b,c
Relators:
   cc
   bbbbbb
   bCbCbCbCbC
   abab
   acacacac
   aaaabAAAcaaaabAAAcaaaabAAAc
>>> O = Orbifold(os.path.join(test_files_paths[0], '1_1^4.84.orb'), remove_finite_vertices = False)
>>> O._orb_cone_fill(2.0, 0)
>>> O._orb_cone_fill(3.0, 1)
>>> O._orb_cone_fill(4.0, 2)
>>> O._orb_cone_fill(5.0, 3)
>>> O._orb_cone_fill(6.0, 4)
>>> O._orb_cone_fill(2.0, 5)
>>> O.solution_type()
'partially flat tetrahedra'
>>> O.volume() # doctest: +NUMERIC9
5.43335845048923

>>> O._orb_cone_fill(2.1, 0)
>>> O._orb_singular_edge_info('singular_order')
[2.1, 3.0, 4.0, 5.0, 6.0, 2.0]
>>> O
1_1^4.84(2.1)(3)(4)(5)(6)(2)
>>> O.cusp_info()
[]
>>> O._orb_vertex_info()
[Vertex 0 : spherical link with cone points of order [2.1, 3.0, 2.0] and orbifold Euler char 0.309524,
 Vertex 1 : spherical link with cone points of order [3.0, 4.0, 6.0] and orbifold Euler char -0.25,
 Vertex 2 : spherical link with cone points of order [4.0, 5.0, 2.0] and orbifold Euler char -0.05,
 Vertex 3 : spherical link with cone points of order [2.1, 5.0, 6.0] and orbifold Euler char -0.157143]
>>> O._orb_singular_edge_info()
[Edge 0 : Singular of order = 2.1,
 Edge 1 : Singular of order = 3,
 Edge 2 : Singular of order = 4,
 Edge 3 : Singular of order = 5,
 Edge 4 : Singular of order = 6,
 Edge 5 : Singular of order = 2]
>>> O.solution_type()
'all tetrahedra positively oriented'
>>> O.volume() # doctest: +NUMERIC9
5.67904978263216

Non-integral cone fillings. This give the free group of three generators in Orb, but not for us:

Really skip this: >>> O.fundamental_group() # doctest: +SKIP

>>> T = SnapPy._orb_test_triangulating_diagram(os.path.join(test_files_paths[0], '6_5^2.7.orb'))
>>> T._orb_cone_fill([2,2,2])
>>> T.fundamental_group(False)
Generators:
   a,b,c,d,e,f
Relators:
   cdFcf
   Fbce
   dAB
   cc
   BAe
   DeDe
   AA

>>> T = SnapPy._orb_test_triangulating_diagram(os.path.join(test_files_paths[0], '1_1^2.1.orb'))
>>> T._orb_cone_fill([2,3,4])
>>> T.fundamental_group(False)
Generators:
   a,b,c
Relators:
   CAB
   bcBAbcBAbcBAbcBA
   CC
   AAA

>>> SnapPy._orb_set_use_orb_conventions(False)

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

import os
from ..extensions.Orb import Orbifold
from .files import __path__ as test_files_paths
