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
