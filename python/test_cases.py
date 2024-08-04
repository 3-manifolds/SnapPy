"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

Doc tests that we do not want to show up in the documentation.

isometry_signature
------------------

>>> M = Manifold('m004(1,2)')
>>> M.isometry_signature()
'cPcbbbiht(3,2)'

sage: M.isometry_signature(verified=True)
'cPcbbbiht(3,2)'

>>> M = ManifoldHP('m004(1,2)')
>>> M.isometry_signature()
'cPcbbbiht(3,2)'

sage: M.isometry_signature(verified=True) # Test bug reported by Nathan
'cPcbbbiht(3,2)'

# Cases where the drilled manifold had a canonical cell decomposition
# with non-tetrahedral cells.

>>> M = Manifold('m137(3,2)')
>>> M.isometry_signature()
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'

sage: M.isometry_signature(verified=True)
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'

>>> M = ManifoldHP('m137(3,2)')
>>> M.isometry_signature()
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'

sage: M.isometry_signature(verified=True)
'sLLvwzvQPAQPQccghmiljkpmqnoorqrrqfafaoaqoofaoooqqaf(3,2)'

Class hierarchy
---------------

>>> isinstance(Manifold("m004"), Triangulation)
True

>>> isinstance(ManifoldHP("m004"), TriangulationHP)
True

Low precision and high precision comparisons
--------------------------------------------
>>> M   = Manifold("m004")
>>> Mhp = Manifold("m004")
>>> N   = Manifold("m003")
>>> Nhp = Manifold("m003")
>>> M.is_isometric_to(M)
True
>>> M.is_isometric_to(Mhp)
True
>>> Mhp.is_isometric_to(M)
True
>>> Mhp.is_isometric_to(Mhp)
True
>>> M.is_isometric_to(N)
False
>>> M.is_isometric_to(Nhp)
False
>>> Mhp.is_isometric_to(N)
False
>>> Mhp.is_isometric_to(Nhp)
False

>>> O = Triangulation("mvvLALQQQhfghjjlilkjklaaaaaffffffff",
... remove_finite_vertices = False)
>>> O.has_finite_vertices()
True
>>> Ohp = TriangulationHP("mvvLALQQQhfghjjlilkjklaaaaaffffffff",
... remove_finite_vertices = False)
>>> Ohp.has_finite_vertices()
True

>>> len(O.isomorphisms_to(O))
8
>>> len(O.isomorphisms_to(Ohp))
8
>>> len(Ohp.isomorphisms_to(O))
8
>>> len(Ohp.isomorphisms_to(Ohp))
8
>>> len(M.isomorphisms_to(O))
0
>>> len(M.isomorphisms_to(Ohp))
0
>>> len(Mhp.isomorphisms_to(O))
0
>>> len(Mhp.isomorphisms_to(Ohp))
0

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from . import Manifold, ManifoldHP, Triangulation, TriangulationHP