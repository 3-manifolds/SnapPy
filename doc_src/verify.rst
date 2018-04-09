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
   sage: K = M.canonical_retriangulation(verified = True)
   sage: len(K.isomorphisms_to(K)) # Certified size of isometry group
   8
 
  **Remark:** For the case of non-tetrahedral canonical cell, exact values
  are used which are found
  using the   `LLL-algorithm 
  <http://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm>`_
  and then verified using exact computations. These computations can be slow. A massive speed-up was achieved by
  recent improvements so that the computation of the isometry signature of any manifold in ``OrientableCuspedCensus``
  takes at most a couple of seconds, typically, far less. Manifolds with more simplices might require setting
  a higher value for 
  :py:attr:`exact_bits_prec_and_degrees`.

* The isometry signature which is a complete invariant of the isometry type
  of a cusped hyperbolic manifold (i.e., two manifolds are isometric if and only
  if they have the same isometry signature)::

   sage: M = Manifold("m412")
   sage: M.isometry_signature(verified = True)
   'mvvLALQQQhfghjjlilkjklaaaaaffffffff'

  The isometry signature can be strengthened to include the peripheral curves
  such that it is a complete invariant of a hyperbolic link::

   sage: M = Manifold("L5a1")
   sage: M.isometry_signature(of_link = True, verified = True)
   'eLPkbdcddhgggb_baCbbaCb'

  **Remark:** The isometry signature is based on the canonical
  retriangulation so the same warning applies.

* Complex intervals for the translations of meridian and longitude with respect
  to disjoint cusp neighborhoods::

   sage: M = Manifold("s441")
   sage: M.cusp_translations()
   [(0.30456698? + 1.38179990?*I, 1.84652839?),
    (0.30456698? + 1.38179990?*I, 1.84652839?)]   

  These can be used to find all potential exceptional slopes which by the
  `Agol's <http://arxiv.org/abs/math/9906183>`_ and 
  `Lackenby's <http://arxiv.org/abs/math/9808120>`_ 6-Theorem must have
  a translation less or equal to 6.

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

   sage: Manifold("m003").canonical_retriangulation(verified = True).isometry_signature()
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
