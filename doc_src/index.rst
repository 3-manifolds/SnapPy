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
