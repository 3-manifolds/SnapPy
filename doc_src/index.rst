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

* Version 3.0 (April 2021): New features include:

  - Incorporates Zoltán Szabó's `program
    <https://web.math.princeton.edu/~szabo/HFKcalc.html>`_ for
    computing Knot Floer homology, see :meth:`knot_floer_homology
    <spherogram.Link.knot_floer_homology>`.  This can compute the
    Seifert genus of a 25 crossing knot in mere seconds!

  - Topological slice obstructions of Herald-Kirk-Livingston, see
    :meth:`slice_obstruction_HKL <snappy.Manifold.slice_obstruction_HKL>`.

  - Faster "local" algorithm for :meth:`jones_polynomial
    <spherogram.Link.jones_polynomial>`.

  - `Cohomology fractals <https://arxiv.org/abs/2010.05840>`_ added to
    :meth:`inside_view <snappy.Manifold.inside_view>`.

  - Convention changes: Sign of knot signature (now positive knots have
    negative signatures), choice of braid generators (now positive
    generators give positive crossings).

  - Updates to methods :meth:`cusp_translations
    <snappy.Manifold.cusp_translations>`, :meth:`cusp_areas
    <snappy.Manifold.cusp_areas>`, :meth:`short_slopes
    <snappy.Manifold.short_slopes>`. Also :meth:`Link <spherogram.Link>`
    now accepts DT codes.

  - Support for SageMath 9.3, Python 3.9, and macOS Big Sur.

  - macOS app now code-signed and notarized.

  - SnapPy now requires Python 3.6 or newer.


* Version 2.8 (June 2020): New features include:

  - Raytraced interior views of a hyperbolic 3-manifold via the 
    :meth:`inside_view <snappy.Manifold.inside_view>` method, see also
    `images <https://im.icerm.brown.edu/portfolio/snappy-views/>`_ and
    `demo video <https://youtu.be/CAERhmUCkRs>`_.

  - :doc:`verify`: Several new features:

    * Complex volume (and thus the Chern-Simons invariant) for both
      cusped and closed manifolds, see
      :meth:`complex_volume <snappy.Manifold.complex_volume>`.
      
    * Disjoint cusp neighborhoods by the method :meth:`cusp_areas
      <snappy.Manifold.cusp_areas>` which uses
      :meth:`cusp_area_matrix <snappy.Manifold.cusp_area_matrix>`.

    * Cusp shapes (see :meth:`cusp_info <snappy.Manifold.cusp_info>`).
      
    * Finding all :meth:`short_slopes <snappy.Manifold.short_slopes>`
      in disjoint embedded cusp neighborhoods.

  - The census :class:`HTLinkExteriors <snappy.HTLinkExteriors>` has
    been extended to 15 crossing knots (contributed by Malik
    Obeidin).

  - The census :class:`CensusKnots <snappy.CensusKnots>` has been
    extended to triangulations with 9 ideal tetrahedra.

  - Support for SageMath 9.0 and macOS Catalina.

  - Development moved to `GitHub <https://github.com/3-manifolds>`_.

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
