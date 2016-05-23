.. SnapPy news

============
News
============

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
  (ManifoldHP) <manifoldhp.html>`_ which compute hyperbolic stuctures
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
   <http://bitbucket.org//Mark_Bell/twister/>`_; a useful `users guide <http://bitbucket.org//Mark_Bell/twister/raw/tip/docs/Twister.pdf>`_.

* Version 1.3 (February 2011) has a completely redesigned cusp horoball viewer and many bug fixes!

* Version 1.2 (December 2010).

* Version 1.1 (February 2010).

* Version 1.0 (August 2009) Initial version. 
