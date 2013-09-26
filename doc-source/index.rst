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

  - A video tutorial of the new features `is available <http://youtu.be/bCYe_a48viA>`_. 

* Version 1.6 (August 2012) includes a `new way to make links
  <spherogram.html>`_ and some support for `arbitrary precision calculation <snap.html>`_.  

* Version 1.5 (February 2012) includes `much improved manifold
  censuses <censuses.html>`_.  

*  Version 1.3.10 (July 2011) incorporates `Twister
   <http://surfacebundles.wordpress.com/>`_.

* `Complete version history <news.html>`_.

Documention
============

.. toctree::
   :maxdepth: 1

   screenshots   
   news
   installing
   tutorial
   snappy
   plink
   spherogram
   snap
   todo

Credits
=============

Written by `Marc Culler <http://www.math.uic.edu/~culler>`_ and
`Nathan Dunfield <http://dunfield.info>`_ with additions from Matthias
Goerner and data from Morwen Thistlethwaite.  Uses the SnapPea kernel
written by `Jeff Weeks <http://www.geometrygames.org>`_, and includes
`Twister <http://surfacebundles.wordpress.com/>`_ by Mark Bell, Tracy
Hall and Saul Schleimer.

Released under the terms of the `GNU General Public License
<http://www.gnu.org/licenses/gpl-2.0.txt>`_, version 2 or later.
Please cite as:

M. Culler, N. M. Dunfield, and J. R. Weeks. SnapPy, a computer program
for studying the geometry and topology of 3-manifolds, http://snappy.computop.org 

The development of SnapPy was partially supported by grants from the
National Science Foundation, including DMS-0906155 and DMS-0707136.
Any opinions, findings, and conclusions or recommendations expressed
on this site are those of the author(s) and do not necessarily
reflect the views of the National Science Foundation.
