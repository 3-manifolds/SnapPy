"""

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

Cusp orientability
------------------
>>> ' '.join(my_cusp_orientabilities(T) for T in NonorientableCuspedCensus[:400])
'K K KK K K K KK K K KK K K K KK K K K KK K KK KK K K KK K KK K KK K K KK KK KK KK KK tKK KK KK Kt t K K KK K K K K K K K K K KK KK K K K K K K K K K K K K K K K KK KK K K K K K K K K K K K K K K K K KK t K K K K KK KK KK KK KK KK KK KK K K KK KK K Kt KK K KK Kt K K K K K KK K K K K K K KK K K KK KK K K K KK KK K KK KK K KK K KK KK KK KK KK KK KK K K K K K K K K K K K K K K K K K K K K K K K K K K K K K K K K t K K K K K K K K K K KK KK K K K K K K K K K K K K K K K K K K K K K K K K K KK K K K K K KK KK K K K K KK KK KK KK K K K t KK KK KK K K K K K K K K KK KK t KK K K KK KK K K KK KK K K KKK KKK Kt Kt KK KK Kt tK K K KK K K K K K K K K K K KK KK K K K K K K K K K K K K K K K K KK K KK K K KK KK KK K KKt KKt KK Kt KK K K KK K K KK KK KK K K K K K K KK KK K K K K K K K K KK K K KK K KK K KK KK K tK KK tK K KK K K K K KK KK K K KK KK KK KKK tK KK t tt KKKK tKK KK KK K K KK K K K K K K KK K K K KK K K KK K KK K K K K K KK KK K'

CuspInfo
--------
>>> M = Triangulation('v3227(0,0)(1,2)(3,3)')
>>> c = M.cusp_info(2)
>>> sorted(c.keys())
['cone_point_orders', 'cone_point_singular_edge_indices', 'cusp_index', 'euler_characteristic', 'filling', 'is_complete', 'is_cusp', 'is_finite', 'orbifold_euler_characteristic', 'orientable', 'singular_order', 'topology']
>>> c.cone_point_orders
[]
>>> c.is_cusp
True
>>> c['singular_order'] # The new key
3
>>> c.singular_order # The new key
3
>>> 'singular_order' in dir(c)
True
>>> 'singular_order' in c.keys()
True
>>> c['singularity_index'] # The old key should still work (via CuspInfo._obsolete)
3
>>> c.singularity_index # The old key should still work (via CuspInfo._obsolete)
3
>>> 'singularity_index' in dir(c)
False
>>> 'singularity_index' in c.keys()
False

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from .. import Manifold, ManifoldHP, Triangulation, TriangulationHP, NonorientableCuspedCensus

def my_cusp_orientabilities(T):
    T._testing_compute_cusp_orientabilities()
    return ''.join(topology[0] for topology in T.cusp_info('topology'))
