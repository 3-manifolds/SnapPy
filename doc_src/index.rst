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

* `Complete version history <news.html>`_.

Documentation
==============

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
=============

Written by `Marc Culler <http://www.math.uic.edu/~culler>`_, `Nathan
Dunfield <http://dunfield.info>`_, and `Matthias Goerner
<http://www.unhyperbolic.org/>`_ using the SnapPea kernel written by
`Jeff Weeks <http://www.geometrygames.org>`_, with contributions from
`many others <credits.html>`_.  If you use SnapPy in your work, please
`cite it as described here <credits.html#citing-snappy>`_.  If you
encounter problems with SnapPy, `please report them <bugs.html>`_. 

Released under the terms of the `GNU General Public License
<http://www.gnu.org/licenses/gpl-2.0.txt>`_, version 2 or later.

The development of SnapPy was partially supported by grants from the
National Science Foundation, including, DMS-0707136, DMS-0906155,
DMS-1105476, and DMS-1510204.  Any opinions, findings, and conclusions
or recommendations expressed on this site are those of the author(s)
and do not necessarily reflect the views of the National Science
Foundation.
