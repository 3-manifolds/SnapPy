.. SnapPy documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======
SnapPy
======

What is SnapPy?
===============

..  image:: images/SnapPy-196.png
    :align: right
    :alt: SnapPy logo


SnapPy is a program for studying the topology and geometry of
3-manifolds, with a focus on hyperbolic structures.  It runs on Mac OS
X, Linux, and Windows, and combines a link editor and 3D-graphics for
Dirichlet domains and cusp neighborhoods with a powerful command-line
interface based on the Python_ programming language. You can see it
:doc:`in action<screenshots>`, learn how to :doc:`install<installing>`
it, and watch the :doc:`tutorial<tutorial>`.

.. _Python: http://python.org

News
====

* Version 3.3 (January 2026) and 3.3.1 (February 2026):

  - :class:`Link <spherogram.Link>` now supports band moves and can
    search for ribbon disks and ribbon concordances. See
    :meth:`ribbon_concordant_links
    <spherogram.Link.ribbon_concordant_links>` and :meth:`add_band
    <spherogram.Link.add_band>`.  From `[DG]
    <https://arXiv.org/abs/2512.21825>`_.

  - New census :class:`RibbonLinks <snappy.RibbonLinks>`.  From `[DG]
    <https://arXiv.org/abs/2512.21825>`_.

  - Additional slice obstructions added to
    :meth:`slice_obstruction_HKL
    <snappy.Triangulation.slice_obstruction_HKL>`.
    From `[DG] <https://arXiv.org/abs/2512.21825>`_.

  - The Fox-Milnor slice obstruction is now available as
    :meth:`fox_milnor_test <snappy.Triangulation.fox_milnor_test>`.
    From `[DG] <https://arXiv.org/abs/2512.21825>`_.


  - The census :class:`OrientableCuspedCensus
    <snappy.OrientableCuspedCensus>` has been extended to 10 ideal
    tetrahedra by `[Li] <https://arXiv.org/abs/2512.02142>`_,
    adding 150,000 new manifolds.
	  
  - Margulis number can now be computed with
    :meth:`margulis <snappy.Manifold.margulis>`.

  - The upper bounds on the bridge number of a link from `[BKVV2020]
    <https://dx.doi.org/10.4310/CAG.2020.v28.n2.a2>`_ and `[BKP2025]
    <https://arxiv.org/abs/2504.10517>`_ are available as
    :meth:`bridge_upper_bound <spherogram.Link.bridge_upper_bound>`.

  - :meth:`triangulation_isosig <snappy.Manifold.triangulation_isosig>` now
    preserves the orientation of the core curve rather than the Dehn-filling
    curve - thus, preserving the (unoriented) spun-triangulation structure.
    This only changes the result when ``ignore_curves=True`` (and
    ``ignore_orientation`` and ``ignore_filling_orientations`` have their
    default values).

  - Fixing a bug in
    :meth:`length_spectrum_alt <snappy.Manifold.length_spectrum_alt>`
    and
    :meth:`length_spectrum_alt_gen <snappy.Manifold.length_spectrum_alt_gen>`.
    The estimate for the upper bound of the global spine radius was not
    correct in all cases. Note that we expect this bug to be very hard to hit
    and do not have a single example where the result of these methods was
    wrong.

  - :meth:`length_spectrum_alt <snappy.Manifold.length_spectrum_alt>` can be
    given ``count`` and ``max_len`` at the same time. In that case, it stops
    when enough geodesics have been listed to ensure that they include the
    ``count`` shortest geodesics or that they include all geodesics shorter
    than ``max_len``.

  - Acceleration to
    :meth:`length_spectrum_alt <snappy.Manifold.length_spectrum_alt>` when
    the next geodesic has large gap to the last geodesic in the result.
    It uses the new ``include_intermediates=True`` flag of
    :meth:`length_spectrum_alt_gen <snappy.Manifold.length_spectrum_alt_gen>`.

  - Computing maximal cusp areas in an unbiased way uses a simpler algorithm
    which also returns tighter intervals for verified computations.
    This affects :meth:`cusp_areas <snappy.Manifold.cusp_areas>` and related
    methods.

  - Support for Python 3.14 and SageMath 10.8.
    


* Version 3.2 (January 2025):

  - An alternative implementation of the length spectrum. It supports verified
    computations and is significantly faster in some cases. See
    :meth:`length_spectrum_alt <snappy.Manifold.length_spectrum_alt>`
    and
    :meth:`length_spectrum_alt_gen <snappy.Manifold.length_spectrum_alt_gen>`.

  - :meth:`isometry_signature <snappy.Manifold.isometry_signature>` now also
    works for closed manifolds.

  - A paper plane or eye ball in the
    :meth:`inside_view <snappy.Manifold.inside_view>` showing the position
    and bearing of the user in the hyperbolic manifold. By default, we now show
    the paper plane instead of the edges of the triangulation so that the view
    is intrinsic to the manifold. The geodesics window now also features a
    button to put the camera onto a geodesic. Here are some examples:

    .. list-table::
       :width: 100%
       :class: borderless

       * - .. image:: images/o9_00000_systole_paper_plane.jpg
              :width: 90%
              :align: center
              :alt: Paper plane near shortest geodesic of o9_00000
         - .. image:: images/o9_00000_systole_paper_plane_closer.jpg
              :width: 90%
              :align: center
              :alt: Paper plane very near shortest geodesic of o9_00000
         - .. image:: images/m004_paper_plane_on_systole.jpg
              :width: 90%
              :align: center
              :alt: Paper plane on shortest geodesic of m004
         - .. image:: images/m125_paper_plane.jpg
              :width: 90%
              :align: center
              :alt: Paper plane coming out of a cusp of m125

  - A faster and more robust algorithm to the compute maximal cusp area matrix.
    The new algorithm is now the default for
    :meth:`~snappy.Manifold.cusp_area_matrix`,
    :meth:`~snappy.Manifold.cusp_areas`,
    :meth:`~snappy.Manifold.short_slopes` and
    :meth:`~snappy.Manifold.cusp_translations`.

  - New options ``ignore_curves`` and ``ignore_filling_orientations``
    for :meth:`~snappy.Triangulation.triangulation_isosig`. Also
    fixing a subtle bug where the filling coefficients returned by
    :meth:`triangulation_isosig <snappy.Triangulation.triangulation_isosig>` were
    not canonical when ``ignore_curve_orientations = True``.

  - :meth:`~snappy.Manifold.canonical_retriangulation`
    and
    :meth:`~snappy.Manifold.isometry_signature` fail with
    exceptions rather than silently returning ``None``. In particular, it now
    safe to compare isometry signatures (without further checks) to determine
    whether M and N are isometric hyperbolic manifolds::

      >>> M.isometry_signature(verified=True) == N.isometry_signature(verified=True)

  - Bug fix to :meth:`slice_obstruction_HKL
    <snappy.Manifold.slice_obstruction_HKL>`: earlier versions
    incorrectly allowed looking at mod 2 homology of the branched
    cover, where the theory does not directly apply.

  - New self-contained SnapPy application for Linux.

  - Support for Python 3.13 and SageMath 10.5.

* :doc:`Complete version history <news>`.

Documentation
=============

.. toctree::
   :maxdepth: 1

   installing
   screenshots
   tutorial
   snappy
   plink
   spherogram
   snap
   verify
   other
   news
   credits
   bugs
   todo
   development

Credits
=======

Written by `Marc Culler <http://www.math.uic.edu/~culler>`_, `Nathan
Dunfield <http://dunfield.info>`_, and `Matthias Goerner
<http://www.unhyperbolic.org/>`_ using the SnapPea kernel written by
`Jeff Weeks <http://www.geometrygames.org>`_, with contributions from
:doc:`many others <credits>`.  If you use SnapPy in your work, please
:ref:`cite it as described here <credits:Citing SnapPy>`.  If you
encounter problems with SnapPy, :doc:`please report them <bugs>`.

Released under the terms of the `GNU General Public License
<http://www.gnu.org/licenses/gpl-2.0.txt>`_, version 2 or later.

The development of SnapPy was partially supported by grants from the
National Science Foundation, including DMS-0707136, DMS-0906155,
DMS-1105476, DMS-1510204, DMS-1811156, and DMS-2303572, and the
Institute for Computational and Experimental Research in
Mathematics. Any opinions, findings, and conclusions or
recommendations expressed on this site are those of the authors and do
not necessarily reflect the views of the National Science Foundation.
