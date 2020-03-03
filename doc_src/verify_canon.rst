Canonical retriangulation and isometry signature
--------------------------------------------------------

The canonical retriangulation is a close relative to the canonical cell
decomposition defined by `Epstein and Penner 
<https://projecteuclid.org/euclid.jdg/1214441650>`_.
Like the canonical cell decomposition, it is intrinsic to
a hyperbolic manifold M and is (up to combinatorial isomorphism
relabeling the tetrahedra and vertices) completely determined by the
isometry type of a hyperbolic manifold. Unlike the canonical cell decomposition,
the canonical retriangulation always conists entirely of tetrahedra which makes
it more amenable for many computations by SnapPy.

If the canonical cell decompositon of manifold M has only tetrahedral cells,
we define the canonical retriangulation to be the canonical cell decomposition.
In this case, the canonical retriangulation consists of ideal hyperbolic
tetrahedra and the ``canonical_retriangulation`` method returns a
SnapPy manifold. Example::

   sage: M = Manifold("m015")
   sage: K = M.canonical_retriangulation(verified = True)
   sage: K.has_finite_vertices() # False iff all canonical cells tetrahedral
   False

If the canonical cell decomposition has non-tetrahedral cells, we turn it into
a topological triangulation as follows: pick a point (called center) in each
3-cell. "Suspend" each 2-cell (which is an ideal n-gon) between
the centers of the two neighboring 3-cells. These suspensions form a
decomposition of M into topological "diamonds". Each diamond can be split along
its central axis into n tetrahedra. This introduces finite vertices, thus
the ``verified_canonical_retriangulation`` method returns only a SnapPy
triangulation. Example (canonical cell is a cube)::
 
   sage: M = Manifold("m412")
   sage: K = M.canonical_retriangulation(verified = True)
   sage: K.has_finite_vertices()
   True

The canonical retriangulation can be used to certifiably find all isometries
of a manifold::
  
   sage: K.isomorphisms_to(K)
   [0 -> 1  1 -> 0
    [1 0]   [1 0] 
    [0 1]   [0 1] 
    Extends to link,
    ...
    Extends to link]
   sage: len(K.isomorphisms_to(K))
   8

Recall that the *isomorphism
signature* is a complete invariant of the combinatorial
isomorphism type of a triangulation that was defined by `Burton
<http://arxiv.org/abs/1110.6080>`_. We can compute the isomorphism signature
of the canonical retriangulation::

   sage: Manifold("m003").canonical_retriangulation(verified = True).triangulation_isosig()
   'cPcbbbdxm'

The resulting invariant was called *isometry signature* by
`Goerner <http://arxiv.org/abs/1502.00383>`_ and, for convenience, can be
accessed by::
   
   sage: Manifold("m003").isometry_signature(verified = True)
   'cPcbbbdxm'

It is a complete invariant of the isometry type of a hyperbolic manifold.
Thus it can be used to easily identify isometric manifolds
(here, the last two manifolds have the same isometry signature and thus
have to be isomorphic)::

   sage: Manifold("m003").isometry_signature(verified = True)
   'cPcbbbdxm'
   sage: Manifold("m004").isometry_signature(verified = True)
   'cPcbbbiht'
   sage: Manifold("4_1").isometry_signature(verified = True)
   'cPcbbbiht'
   sage: Manifold("m004").isometry_signature(verified = True) == Manifold("4_1").isometry_signature(verified = True)
   True


Other applications of the canonical retriangulation include the detection of
2-bridge knots.

========================================
Verifiying the canonical retriangulation
========================================

..   autofunction:: snappy.verify.verified_canonical_retriangulation
