Verified computations
========================================

Overview
--------

When used inside `Sage <http://sagemath.org>`_, SnapPy can verify the
following computations:

* Complex intervals for the shapes that are guaranteed to contain a true
  but not necessarily geometric solution to the rectangular gluing equations::

      sage: M = Manifold("m015(3,1)")
      sage: M.tetrahedra_shapes('rect', intervals=True)
      [0.625222762246? + 3.177940133813?*I,
       -0.0075523593782? + 0.5131157955971?*I,
       0.6515818912107? - 0.1955023488930?*I]

  (Specify :py:attr:`bits_prec` or :py:attr:`dec_prec` for higher precision intervals.)
  
* Verify the hyperbolicity
  of an orientable 3-manifold giving complex intervals for the
  shapes corresponding to a hyperbolic structure or holonomy representation with
  :py:meth:`~snappy.Manifold.verify_hyperbolicity`::

   sage: M = Manifold("m015")
   sage: M.verify_hyperbolicity()
   (True,
    [0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I])
   sage: M.verify_hyperbolicity(holonomy=True)[1].SL2C('a')
   [-0.324717957? - 1.124559024?*I -0.704807293? + 0.398888830?*I]
   [ 1.409614585? - 0.797777659?*I       -1.000000000? + 0.?e-9*I]

* Intervals for the volume and complex volume of a hyperbolic orientable 3-manifold::

   sage: M = Manifold("m003(-3,1)")
   sage: M.volume(verified=True, bits_prec = 100)
   0.942707362776927720921299603?
   sage: M = Manifold("m015")
   sage: M.complex_volume(verified_modulo_2_torsion=True)
   2.8281220883? + 1.9106738240?*I

  (Note that when using verified computation, the Chern-Simons invariant is only computed
  modulo pi^2/2 even though it is defined modulo pi^2.)

* Give the :doc:`canonical retriangulation <verify_canon>` (a close
  relative to the canonical cell decomposition) of a cusped hyperbolic
  manifold using intervals or exact arithmetic if necessary with
  :py:meth:`~snappy.Manifold.canonical_retriangulation`::

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

* The :doc:`isometry signature <verify_canon>` which is a complete invariant of the isometry type
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

  See :py:meth:`~snappy.Manifold.isometry_signature` for details.

  **Remark:** The isometry signature is based on the canonical
  retriangulation so the same warning applies.

* The maximal cusp area matrix which characterizes the configuration
  space of disjoint cusp neighborhoods with
  :py:meth:`~snappy.Manifold.cusp_area_matrix`::

   sage: M=Manifold("m203")
   sage: M.cusp_area_matrix(method='maximal', verified=True)
   [   27.000000? 9.0000000000?]
   [9.0000000000?   27.0000000?]

  In this example, the cusp neighborhood about cusp 0 or 1 is only
  embedded if and only if its area is less than sqrt(27). The cusp
  neighborhood about cusp 0 is only disjoint from the one about cusp 1
  if and only if the product of their areas is less than 9.

* Compute areas for disjoint cusp neighborhoods with
  :py:meth:`~snappy.Manifold.cusp_areas`::

   sage: M=Manifold("m203")
   sage: M.cusp_areas(policy = 'unbiased', method='maximal', verified = True)
   [3.00000000000?, 3.00000000000?]

  With the above parameters, the result is intrinsic to the hyperbolic
  manifold with labeled cusped.

* Find all slopes of length less or equal to 6 when measured on the boundary
  of disjoint cusp neighborhoods::

   sage: M=Manifold("m203")
   sage: M.short_slopes(policy = 'unbiased', method='maximal', verified = True)
   [[(1, 0), ...,  (1, 2)], [(1, 0), ...,  (1, 2)]]

  First block has all short slopes for first cusp, ..., see
  :py:meth:`~snappy.Manifold.short_slopes` for details.

  By `Agol's <http://arxiv.org/abs/math/9906183>`_ and `Lackenby's
  <http://arxiv.org/abs/math/9808120>`_ 6-Theorem any Dehn-filling
  resulting in a non-hyperbolic manifold must contain one of the above
  slopes.  Thus, :py:meth:`~snappy.Manifold.short_slopes` can be used
  to implement the techniques to find exceptional Dehn surgeries
  (`arXiv:1109.0903 <https://arxiv.org/abs/1109.0903>`_ and
  `arXiv:1310.3472 <https://arxiv.org/abs/1310.3472>`_).

This is all based on a reimplementation of `HIKMOT
<http://www.oishi.info.waseda.ac.jp/~takayasu/hikmot/>`_ which
pioneered the use of interval methods for hyperbolic manifolds (also see
`Zgliczynski's notes <http://ww2.ii.uj.edu.pl/~zgliczyn/cap07/krawczyk.pdf>`_). It
can be used in a way very similar to HIKMOT, but uses Sage's complex
interval types for certification. It furthermore makes use of code by 
`Dunfield, Hoffman, Licata <http://arxiv.org/abs/1407.7827/>`_. The code to
compute the isomorphism signature was ported over from
`Regina <https://regina-normal.github.io/>`_.

This verification code was contributed by Matthias Goerner.  


Verified computation topics
---------------------------

.. toctree::
   :maxdepth: 1

   verify_canon   
   verify_internals
