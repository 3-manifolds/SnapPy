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

  - A video tutorial of the new features `is available <http://youtu.be/bCYe_a48viA>`_. 

* Version 1.6 (August 2012) includes a `new way to make links
  <spherogram.html>`_ and some support for `arbitrary precision calculation <snap.html>`_.  

* Version 1.5 (February 2012) includes `much improved manifold
  censuses <censuses.html>`_.  

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
   other
   news
   credits
   todo
   development

Credits
=============

Written by `Marc Culler <http://www.math.uic.edu/~culler>`_ and
`Nathan Dunfield <http://dunfield.info>`_ using the SnapPea kernel
written by `Jeff Weeks <http://www.geometrygames.org>`_, with
contributions from `many others <credits.html>`_.   If you use SnapPy in your work, please `cite it as described here <credits.html#citing-snappy>`_.

Released under the terms of the `GNU General Public License
<http://www.gnu.org/licenses/gpl-2.0.txt>`_, version 2 or later.

The development of SnapPy was partially supported by grants from the
National Science Foundation, including DMS-0906155 and DMS-0707136.
Any opinions, findings, and conclusions or recommendations expressed
on this site are those of the author(s) and do not necessarily
reflect the views of the National Science Foundation.
