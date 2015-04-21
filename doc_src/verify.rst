Verified computations
========================================

When used inside `Sage <http://sagemath.org>`_, SnapPy can verify the
following computations:

* Complex intervals for the shapes that are guaranteed to contain a true
  solution to the rectangular gluing equations::

      sage: M = Manifold("m015(3,1)")
      sage: M.tetrahedra_shapes('rect', intervals=True)
      [0.625222762246? + 3.177940133813?*I,
       -0.0075523593782? + 0.5131157955971?*I,
       0.6515818912107? - 0.1955023488930?*I]

  (Specify :py:attr:`bits_prec` or :py:attr:`dec_prec` for higher precision intervals.)
  
* Verify the hyperbolicity of an orientable 3-manifold using intervals::

   sage: M = Manifold("m015")
   sage: M.verify_hyperbolicity()
   (True,
    [0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I])

* Give the canonical retriangulation (a close relative to the canonical cell
  decomposition) of a cusped hyperbolic manifold using
  intervals or exact arithmetic if necessary::

   sage: M = Manifold("m412")
   sage: K = M.canonical_retriangulation(M, verified = True)
   sage: len(K.isomorphisms_to(K)) # Certified size of isometry group
   8
 
  **Warning:** For the case of non-tetrahedral canonical cell, this code can sometimes
  require setting high values for :py:attr:`exact_bits_prec_and_degrees` for the
  `LLL-algorithm 
  <http://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm>`_
  and be very slow. In the ``OrientableCuspedCensus``, for example,
  ``t11669`` and 9 manifolds with 9 tetrahedra required setting a very high
  :py:attr:`exact_bits_prec_and_degrees` and several minutes of computation time. Future
  work pending on `Sage bug 14164 <http://trac.sagemath.org/ticket/14164>`_ and an
  `untracked Sage bug <https://groups.google.com/forum/#!topic/sage-support/N-O8FAHBQTM>`_.
  will hopefully improve performance a lot.

* The isometry signature which is a complete invariant of the isometry type
  of a cusped hyperbolic manifold (i.e., two manifolds are isometric if and only
  if they have the same isometry signature)::

   sage: M = Manifold("m412")
   sage: M.isometry_signature(verified = True)
   'mvvLALQQQhfghjjlilkjklaaaaaffffffff'

  **Warning:** The isometry signature is based on the canonical
  retriangulation so the same warning applies.

This is all based on a reimplementation of `HIKMOT
<http://www.oishi.info.waseda.ac.jp/~takayasu/hikmot/>`_ which
pioneered the use of interval methods for hyperbolic manifolds. It
can be used in a way very similar to HIKMOT, but uses Sage's complex
interval types and the Newton interval method (instead of the Krawczyk
test) for certification. See
`Zgliczynski's notes <http://ww2.ii.uj.edu.pl/~zgliczyn/cap07/krawczyk.pdf>`_ for a quick
overview of these two tests. It furthermore makes use of code by 
`Dunfield, Hoffman, Licata <http://arxiv.org/abs/1407.7827/>`_. The code to
compute the isomorphism signature was ported over from
`Regina <http://regina.sf.net/>`_.

This verification code was contributed by Matthias Goerner.  

The canonical retriangulation and the isometry signature
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

   sage: Manifold("m003").canonical_retriangulation(verified = True).isomorphism_signature()
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

Methods for verified computaions
--------------------------------

..   autofunction:: snappy.verify.verify_hyperbolicity

..   autofunction:: snappy.verify.verified_canonical_retriangulation

Internals
---------

.. toctree::
   :maxdepth: 1
   
   verify_internals