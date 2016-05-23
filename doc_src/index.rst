.. SnapPy documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


==================
SnapPy
==================

What is SnapPy?
==================

..  image:: images/SnapPy-196.png 
    :align: right


SnapPy is a program for studying the topology and geometry of
3-manifolds, with a focus on hyperbolic structures.  
It runs on Mac OS X, Linux, and Windows, and combines a link editor and 3D-graphics
for Dirichlet domains and cusp neighborhoods with a powerful
command-line interface based on the `Python <http://python.org>`_
programming language. You can see it `in action <screenshots.html>`_,
learn how to `install <installing.html>`_ it, and watch the `tutorial
<tutorial.html>`_.   

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

* `Complete version history <news.html>`_.

Documention
============

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
   todo
   development

Credits
=============

Written by `Marc Culler <http://www.math.uic.edu/~culler>`_,
`Nathan Dunfield <http://dunfield.info>`_, and `Matthias Goerner <http://www.unhyperbolic.org/>`_ using the SnapPea kernel
written by `Jeff Weeks <http://www.geometrygames.org>`_, with
contributions from `many others <credits.html>`_.   If you use SnapPy in your work, please `cite it as described here <credits.html#citing-snappy>`_.

Released under the terms of the `GNU General Public License
<http://www.gnu.org/licenses/gpl-2.0.txt>`_, version 2 or later.

The development of SnapPy was partially supported by grants from the
National Science Foundation, including, DMS-0707136, DMS-0906155,
DMS-1105476, and DMS-1510204.  Any opinions, findings, and conclusions
or recommendations expressed on this site are those of the author(s)
and do not necessarily reflect the views of the National Science
Foundation.
