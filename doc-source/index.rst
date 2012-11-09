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

* Version 1.7 (November 2012) incorporates the `ptolemy module
  <http://www.unhyperbolic.org/ptolemy.html>`_ for studying
  representations of 3-manifold groups into pSL(*N*, **C**).  

* Version 1.6 (August 2012) includes a `new way to make links
  <spherogram.html>`_ and some support for `arbitrary precision calculation <snap.html>`_.  

* Version 1.5 (February 2012) includes `much improved manifold
  censuses <censuses.html>`_.  

* Version 1.4 (December 2011) uses the current release of IPython, which has been completely rewritten.

*  Version 1.3.10 (July 2011) incorporates `Twister <http://surfacebundles.wordpress.com/>`_.

* Version 1.3 (February 2011) has a completely redesigned cusp horoball viewer and many bug fixes!

Documention
============

.. toctree::
   :maxdepth: 1

   screenshots   
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

Released under the terms of the GNU General Public License.  Please
cite as:

M. Culler, N. M. Dunfield, and J. R. Weeks. SnapPy, a computer program
for studying the geometry and topology of 3-manifolds, http://snappy.computop.org 

The development of SnapPy was partially supported by grants from the
National Science Foundation, including DMS-0906155 and DMS-0707136.
Any opinions, findings, and conclusions or recommendations expressed
on this site are those of the author(s) and do not necessarily
reflect the views of the National Science Foundation.
