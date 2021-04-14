.. SnapPy news

====
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

* Version 2.7 (July 2019): New features include:

  - Python 3 is now recommended over Python 2 on all platforms; the
    default Mac and Windows apps use Python 3 rather than
    Python 2. The only difference most users will notice is that one
    must type ``print(blah)`` instead of ``print blah``.

  - `Verified computations <verify.html>`_: performance improvements
    by switching to the Krawczyk test.

  - Support for SageMath 8.8.

  - Installation instructions extensively updated.

  - GUI improvements, especially on macOS. These include improved
    support for dark mode and tabs on macOS Mojave and preliminary
    support for macOS Catalina.
    
* Version 2.6.1 (August 2018): New features include:

  - Support for SageMath 8.3, Python 3.7, and macOS Mojave.

  - Computing `ideals defining SL(2, C) character varieties.
    <additional_classes.html#snappy.HolonomyGroup.character_variety_vars_and_polys>`_
    Contributed by Jean-Philippe Burelle, based on `this paper
    <https://arxiv.org/abs/1703.08241>`_.

  - Many bug fixes. 

* Version 2.6 (Nov 2017): New features include:

  - Support for macOS High Sierra, SageMath 8.1, and Windows systems
    using non-Latin alphabets.

  - Many bug fixes, including improved Python 3 support.

* Version 2.5 (Feb 2017): New features include:

  - Rigorous computation of `hyperbolic volume
    <manifold.html#snappy.Manifold.volume>`_.

  - STL export of Dirichlet domains for 3D printing, contributed by
    Jose Sanchez.

  - Support for Python 3, SageMath 7.5, 7.6, and 8.0, and many more 
    versions of Python on Windows.

  - Much improved infrastructure for testing and distributing SnapPy.

* Version 2.4 (May 2016): New features include:

  - Added `census of Platonic manifolds <platonic_census.html>`_. 

  - Rigorous computation of `cusp translations <manifold.html#snappy.Manifold.cusp_translations>`_.  
  
  - Added decorations to `triangulation isomorphism signatures
    <manifold.html#snappy.Manifold.triangulation_isosig>`_ for
    encoding peripheral curves.
    
  - Faster verification of non-tetrahedral canonical cell decompositions.
  
  - Improvements to the `link and planar diagram component
    <spherogram.html>`_, mostly contributed by Malik Obeidin, include:

    * Bar-Natan's super-fast `tangle-based algorithm
      <http://www.math.toronto.edu/drorbn/Talks/Aarhus-1507/>`_ for
      computing the Alexander polynomial.

    * Can now compute the `Seifert matrix
      <spherogram.html#spherogram.Link.seifert_matrix>`_ and express a
      link as a `braid closure <spherogram.html#spherogram.Link.braid_word>`_.

    * Conversion to/from `SageMath links and braids
      <spherogram.html#spherogram.Link.sage_link>`_.

    * Many under-the-hood improvements.  
    
  - New Windows installer. 

* Version 2.3 (March 2015):  New features include:

  - Major improvements to the `link and planar diagram component
    <spherogram.html>`_, including link simplification, random links,
    and better documentation.

  - Basic support for `spun normal surfaces
    <manifold.html#snappy.Manifold.normal_boundary_slopes>`_.

  - New extra features when used inside of Sage:

    * HIKMOT-style `rigorous verification of hyperbolic structures
      <verify.html>`_, 
      contributed by Matthias Goerner.  
      
    * Many `basic knot/link invariants
      <spherogram.html#the-link-class>`_, contributed by Robert
      Lipschitz and Jennet Dickinson.

    * Sage-specific functions are now more easily accessible as
      methods of Manifold and better documented.

    * Improved number field recognition, thanks to Matthias.  
      
  - Better compatibility with OS X Yosemite and Windows 8.1.

  - Development changes:

    * Major source code reorganization/cleanup.  

    * Source code repository moved to `Bitbucket
      <https://bitbucket.org/t3m>`_.

    * Python modules now hosted on `PyPI
      <https://pypi.python.org/pypi>`_, simplifying `installation <installing.html>`_.  

* Version 2.2 (June 2014): Includes Ben Burton's `census of
  orientable cusped manifolds with 9 tetrahedra. <http://arxiv.org/abs/1405.2695>`_

* Version 2.1 (February 2014): New `high-precision manifolds
  (ManifoldHP) <manifoldhp.html>`_ which compute hyperbolic structures
  (and everything related) in `quad-double (212 bit) <http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf>`_
  precision.

* Version 2.0 (September 2013): Many new features, including:

  - A `manifold browser <manifold.html#snappy.Manifold.browse>`_
    window for easily examining a particular manifold.  

  - Many improvements to the `link editor <plink.html#using-snappy-s-link-editor>`_, including

    * A smoothed view mode with image export to EPS/PDF/SVG/TikZ.

    * Producing a fully editable link from combinatorial data like a DT
      code. 
 
  - `Splitting manifolds <manifold.html#snappy.Manifold.split>`_ along surfaces of non-negative euler
    characteristic. 

  - Generalizing the ptolemy obstruction class to allow computation of
    PGL(3,C)-representations and improving usability of the `ptolemy module
    <http://www.unhyperbolic.org/ptolemy.html>`_.	     

  - `CensusKnots <censuses.html#snappy.CensusKnots>`_ now includes
    knot exteriors with 8 tetrahedra.  

* Version 1.8 (May 2013) improves handling of DT codes and adds the
  `HTLinkExteriors <censuses.html#snappy.HTLinkExteriors>`_ census,
  which provides identification for knots and links up to 14 crossings.

* Version 1.7 (November 2012) incorporates the `ptolemy module
  <http://www.unhyperbolic.org/ptolemy.html>`_ for studying
  representations of 3-manifold groups into pSL(*N*, **C**).  

* Version 1.6 (August 2012) includes a `new way to make links
  <spherogram.html>`_ and some support for `arbitrary precision calculation <snap.html>`_.  

* Version 1.5 (February 2012) includes `much improved manifold
  censuses <censuses.html>`_.  

* Version 1.4 (December 2011) uses the current release of IPython, which has been completely rewritten.

*  Version 1.3.10 (July 2011) incorporates `Twister
   <https://github.com/MarkCBell/twister/>`_.

* Version 1.3 (February 2011) has a completely redesigned cusp horoball viewer and many bug fixes!

* Version 1.2 (December 2010).

* Version 1.1 (February 2010).

* Version 1.0 (August 2009) Initial version. 
