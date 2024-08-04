"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

Doc tests that we do not want to show up in the documentation.

decorated isomorphism signature
-------------------------------

An example where the isosig has changed when we were using the permutation
consistently. The old value was:
'gLMzQbcdefffaelaaai_acbBaabCbbabbBC'

>>> T = Triangulation("L6n1(0,0)(0,0)(0,0)")
>>> T.triangulation_isosig()
'gLMzQbcdefffaelaaai_bCabBBbacBDc'

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

sage: M.isometry_signature(verified=True, exact_bits_prec_and_degrees=[]) # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
RuntimeError: Could not compute or verify canonical retriangulation of drilled manifold. Geodesic was: abCDaDAd.

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

Canonical retriangulation
-------------------------

Some cases that should be rejected

>>> M = Manifold("m004(3,4)")
>>> M.canonical_retriangulation() # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
ValueError: Canonical retriangulation needs at least one unfilled cusp.

sage: M.canonical_retriangulation(verified=True) # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
ValueError: Canonical retriangulation needs at least one unfilled cusp.

>>> M = Manifold("m125(10,12)")
>>> M.canonical_retriangulation() # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
ValueError: Could not compute filled triangulation. Are the filling coefficients co-prime integers?

sage: M.canonical_retriangulation(verified=True) # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
ValueError: Could not compute filled triangulation. Are the filling coefficients co-prime integers?

>>> M = Manifold("m125(11.1,12)")
>>> M.canonical_retriangulation() # doctest: +ELLIPSIS +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
...
ValueError: Could not compute filled triangulation. Are the filling coefficients co-prime integers?

sage: M.canonical_retriangulation(verified=True)
Traceback (most recent call last):
...
ValueError: Could not compute filled triangulation. Are the filling coefficients co-prime integers?

>>> M = Manifold("m125(3,4)")
>>> M.canonical_retriangulation().triangulation_isosig()
'eLAkbbcdddhrhj_BaaB'

sage: M.canonical_retriangulation(verified=True).triangulation_isosig()
'eLAkbbcdddhrhj_BaaB'

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from . import Manifold, ManifoldHP, Triangulation, TriangulationHP
