.. SnapPy documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


==================
SnapPy
==================

What is SnapPy?
==================
SnapPy is a user interface to the SnapPea kernel which runs on Mac OS
X, Linux, and Windows.  SnapPy combines a link editor and 3D-graphics for
Dirichlet domains and cusp neighborhoods with a powerful command-line
interface based on the `Python <http://python.org>`_ programming language.


Installation
==================

- OS X: Just download SnapPy.app and place in applications 
- Linux: 
- Windows: 

Tutorial
==================

The easiest way to learn to use SnapPy is to watch the screencasts
available here:

- Intro and quickstart.
- More advanced features.

The **key** thing to remember when using the SnapPy command shell window is
that you can explore objects using introspection and tab-completion::

     In [1]: Manifold? <hit return-key>
     ...instructions for creating a manifold...

So now we create a manifold::

   In [2]: M = Manifold("m004")

But what can we do with it?  ::

    In [3]: M.<hit tab-key>
    ...list of methods...

What does the "cover" method do? ::
     
     In [7]: M.cover?
     ...description of cover method..


Contents
============

.. toctree::
   :maxdepth: 2
   
   snappy
   plink
   todo

Credits
=============

Written by `Marc Culler <http://math.uic.edu/~culler>`_ and `Nathan
Dunfield <http://dunfield.info>`_.  Uses the SnapPea kernel written by
`Jeff Weeks <http://geometrygames.org>`_.  


