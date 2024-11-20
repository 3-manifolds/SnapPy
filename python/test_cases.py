"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

Doc tests that we do not want to show up in the documentation.

Decorated isomorphism signature
-------------------------------

An example where the isosig doesn't pick up the identity perm
because the lexicographic order of the perm interacts with the
matrices.

>>> T = Triangulation("L6n1(0,0)(0,0)(0,0)")
>>> T.triangulation_isosig()
'gLMzQbcdefffaelaaai_acbBaabCbbabbBC'

Test that slopes are computed correctly.

>>> M=Manifold("L14n63023(-5,1)(5,1)(10,1)")
>>> M.triangulation_isosig(decorated=False, ignore_orientation=False)
'vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt'

The canonical orientation (used to compute the unoriented isosig)
is the reverse of the actual orientation:

>>> M.triangulation_isosig(decorated=False)
'vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth'
>>> Mop = M.copy()
>>> Mop.reverse_orientation()
>>> Mop.triangulation_isosig(decorated=False, ignore_orientation=False)
'vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth'

So we expect the oriented isometry signature to flip when neither the isomorphism
signature nor its decoration capture the orientation.

The following is a baked version of snappy.decorated_isosig.test_slope_transformations().
The 32 calls to ManifoldHP.isometry_signature are just too expensive for the test suite.

>>> for ignore_cusp_ordering in [False, True]:
...     for ignore_curves in [False, True]:
...         for ignore_curve_orientations in [False, True]:
...             for ignore_filling_orientations in [False, True]:
...                 for ignore_orientation in [False, True]:
...                     isosig = M.triangulation_isosig(
...                         ignore_cusp_ordering = ignore_cusp_ordering,
...                         ignore_curves = ignore_curves,
...                         ignore_curve_orientations = ignore_curve_orientations,
...                         ignore_filling_orientations = ignore_filling_orientations,
...                         ignore_orientation = ignore_orientation)
...                     print(isosig)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcbaabBBBbBaBbCbBbCb(-5,1)(5,1)(10,1)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabdBacbbBCbaBBBbabB(-5,1)(5,1)(10,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcbaabBBBbBaBbCbBbCb(-5,1)(5,1)(10,1)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabdBacbbBCbaBBBbabB(-5,1)(5,1)(10,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcbaabBBbBbabBcBbBcB(-5,1)(-5,-1)(-10,-1)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabdbacbbBcBabBBbaBb(5,1)(5,-1)(-10,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcbaabBBbBbabBcBbBcB(-5,1)(5,1)(10,1)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabdbacbbBcBabBBbaBb(5,1)(-5,1)(-10,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcba(-1,-6)(-6,5)(-12,11)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabd(7,1)(3,-4)(-1,-11)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcba(1,6)(-6,5)(-12,11)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabd(7,1)(-3,4)(1,11)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcba(-1,-6)(-6,5)(-12,11)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabd(7,1)(3,-4)(-1,-11)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_dcba(1,6)(-6,5)(-12,11)(0,0)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_cabd(7,1)(-3,4)(1,11)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_BbCbBbCbBbBaabBB(0,0)(10,1)(5,1)(-5,1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_aBBBBacbbBCbbabB(10,1)(-5,1)(5,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_BbCbBbCbBbBaabBB(0,0)(10,1)(5,1)(-5,1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_aBBBBacbbBCbbabB(10,1)(-5,1)(5,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_bBcBbBcBbBbaabBB(0,0)(-10,-1)(-5,-1)(-5,1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_abBBbacbbBcBbaBb(-10,1)(5,1)(5,-1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt_bBcBbBcBbBbaabBB(0,0)(10,1)(5,1)(-5,1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth_abBBbacbbBcBbaBb(-10,1)(5,1)(-5,1)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt(0,0)(-12,11)(-6,5)(6,-1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth(-1,-11)(7,1)(3,-4)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt(0,0)(-12,11)(-6,5)(-6,1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth(1,11)(7,1)(-3,4)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt(0,0)(-12,11)(-6,5)(6,-1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth(-1,-11)(7,1)(3,-4)(0,0)
vLLvvLLMALQQzQQceillmnppqrlmrqtruututiivimllaelaqxrvdoxqltt(0,0)(-12,11)(-6,5)(-6,1)
vLLvLLPwPQLAMPQcefikkmnplkopqrsttutuuiixvimqlippawidlabavth(1,11)(7,1)(-3,4)(0,0)

The special case where we do not encode the matrices and have the identity perm
and need it:

>>> T = Triangulation('ovLMvvPQQQccddlmnijklmnmnlgvfamtvfblhaumx_BBbacBbacBbaBCaB(312,241)(332,441)(3322,1441)(726,541)')
>>> T.homology()
Z/4684873938228
>>> isosig = T.triangulation_isosig(ignore_curves=True, ignore_cusp_ordering=False)
>>> isosig
'ovLMvvPQQQccddlmnijklmnmnlgvfamtvfblhaumx_abcd(71,312)(-1105,332)(-8085,3322)(726,1993)'
>>> Triangulation(isosig).homology()
Z/4684873938228


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

Test isometry_signature's ignore_orientation

>>> M = Manifold("m006")
>>> M.isometry_signature(ignore_orientation=False)
'eLMkaccddjgbaj'
>>> M.isometry_signature(ignore_orientation=True)
'eLAkaccddngbak'
>>> M.isometry_signature() # default value
'eLAkaccddngbak'
>>> M.dehn_fill((3,4))
>>> M.isometry_signature(ignore_orientation=False)
'eLMkaccddjgbaj(-1,4)'
>>> M.isometry_signature(ignore_orientation=True)
'eLAkaccddngbak(-3,4)'

Test isometry_signature's of_link and filling.

>>> Manifold("o9_44206(2,3)").isometry_signature(of_link=True)
'jLLvMQQacggfiihhijkkjkehhtb_abBabBbabaab'

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
>>> M.canonical_retriangulation() # doctest: +ELLIPSIS
Traceback (most recent call last):
...
ValueError: Canonical retriangulation needs all cusps to be complete.

sage: M.canonical_retriangulation(verified=True) # doctest: +ELLIPSIS
Traceback (most recent call last):
...
ValueError: Canonical retriangulation needs all cusps to be complete.

Cusp areas
----------

>>> M = Manifold('o9_44210')
>>> M.cusp_areas(policy='greedy') # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,1]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,1,2]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,2,1]) # doctest: +NUMERIC9
[7.053940530873898, 2.3513135103, 3.7690945490]
>>> M.cusp_areas(policy='greedy', first_cusps=[1,]) # doctest: +NUMERIC9
[2.30025338030798, 10.0315765558665, 0.883442685721903]

Cusp translations
-----------------

>>> M = Manifold("s776")
>>> M.cusp_translations(policy = 'greedy', first_cusps = [], bits_prec = 100) # doctest: +NUMERIC21
[(0.70710678118654752440084436210 + 1.8708286933869706927918743662*I, 2.8284271247461900976033774484), (0.35355339059327376220042218105 + 0.93541434669348534639593718308*I, 1.4142135623730950488016887242), (0.35355339059327376220042218105 + 0.93541434669348534639593718308*I, 1.4142135623730950488016887242)]


"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from . import Manifold, ManifoldHP, Triangulation, TriangulationHP
