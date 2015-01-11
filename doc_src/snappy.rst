.. Documentation of the Python part of SnapPy

The snappy module and its classes
=======================================

SnapPy is centered around a Python interface for SnapPea called
"snappy", and this is what you're interacting with in the main "SnapPy
command shell" window.  The main class is Manifold, which is an ideal
triangulation of the interior of a compact 3-manifold with torus
boundary, where each tetrahedron has has been assigned the geometry of
an ideal tetrahedron in hyperbolic 3-space.  A Dehn-filling can be
specified for each boundary component, allowing the description of
closed 3-manifolds and some orbifolds.  The class Manifold is derived
from the simpler Triangulation class which lacks any geometric
structure.  There are also some additional classes for things like
fundamental groups, Dirichlet domains, etc.  Snappy comes with a large
library of 3-manifolds, some of which are grouped together in censuses.  


.. toctree::
   :maxdepth: 2
   
   manifold
   manifoldhp
   triangulation
   additional_classes
   censuses 

