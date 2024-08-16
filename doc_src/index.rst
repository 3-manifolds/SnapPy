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

* Version 3.2 (??? 2024):

  - :meth:`inside_view <snappy.Manifold.inside_view>` shows the user as a paper
    plane or eye ball. Also adding button to geodesics window to put camera
    onto a geodesic. TODO: PICTURE!

  - An alternative implementation of the length spectrum supporting verified
    computations and being significantly faster in some cases, see
    :meth:`length_spectrum_alt <snappy.Manifold.length_spectrum_alt>`
    and
    :meth:`length_spectrum_alt_gen <snappy.Manifold.length_spectrum_alt_gen>`.

  - :meth:`isometry_signature <snappy.Manifold.isometry_signature>` now also
    works for closed manifolds.

  - It is now safe to call::

        >>> M.isometry_signature(verified=True) == N.isometry_signature(verified=True)

    to determine whether M and N are isometric hyperbolic manifolds since
    :meth:`isometry_signature <snappy.Manifold.isometry_signature>` no longer
    fails silently returning ``None`` but raises an exception if the signature
    could not be computed.

  - Similarly, :meth:`canonical_retriangulation <snappy.Manifold.canonical_retriangulation>` can not longer fail silently returning ``None``.

  - The decoration :meth:`triangulation_isosig <snappy.Triangulation.triangulation_isosig>`
    clearly entangles the cusp indexing from the peripheral curves. Note that this is
    a very subtle change only occurring in very few cases. The result is fully backwards
    compatible in that it can be given to older SnapPy versions to recover the
    triangulation and its decoration.

  - New algorithm to compute maximal cusp area matrix which is faster and more robust
    (Example:
    :meth:`Manifold("otet10_00027").cusp_area_matrix(method='maximal') <snappy.Manifold.cusp_area_matrix>`).

* Versions 3.1 (May 2023) and 3.1.1 (June 2023):

  - A method :meth:`exterior_to_link <snappy.Manifold.exterior_to_link>`
    for going from a link exterior to a link diagram taken from
    `Dunfield-Obeidin-Rudd <https://arxiv.org/abs/2112.03251>`_.

  - Covers now computed by the stand-alone `low_index
    <https://pypi.org/project/low-index/>`_ module, which uses
    multiple processor cores and is typically much faster than the old
    code.  In some cases, it is dramatically faster than even GAP or
    Magma.

  - Added geodesics to the :meth:`inside_view
    <snappy.Manifold.inside_view>`.  Here are some intersecting tubes
    about closed geodesics in the manifold ``v3539(5,1)``:

    .. image:: images/geodesics.jpg
       :width: 50%
       :align: center
       :alt: Geodesic tubes for v3539(5,1)

  - Added drilling any simple geodesic with :meth:`drill_word
    <snappy.Manifold.drill_word>` and :meth:`drill_words
    <snappy.Manifold.drill_words>`, not just those that are
    :meth:`combinatorially simple <snappy.Manifold.dual_curves>`.

  - Added `ignore_orientation` flag to :meth:`triangulation_isosig
    <snappy.Triangulation.triangulation_isosig>`.

  - Added `include_words` flag to :meth:`length_spectrum
    <snappy.Manifold.length_spectrum>` for getting the word
    corresponding to a geodesic which can be given to
    :meth:`drill_word <snappy.Manifold.drill_word>`.

  - Support for Python 3.11 and SageMath 10.0.

  - Modernized styling of the documentation.

* Version 3.0.3 (December 2021):

  - Runs natively on Macs with Apple Silicon processors (M1, M2, and variants).

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
DMS-1105476, DMS-1510204, DMS-1811156, and the Institute for
Computational and Experimental Research in Mathematics. Any opinions,
findings, and conclusions or recommendations expressed on this site
are those of the authors and do not necessarily reflect the views of
the National Science Foundation.
