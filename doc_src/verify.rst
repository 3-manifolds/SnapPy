Verified computations
========================================

When used inside `Sage <http://sagemath.org>`_, SnapPy can use
interval arithmetic to rigorously verify the hyperbolicity of an
orientable 3-manifold::

   sage: M = Manifold("m015")
   sage: M.verify_hyperbolicity()
   (True,
    [0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I])   

It does so by producing intervals that are guaranteed to contain a
true solution to the rectangular gluing equations (give ``bits_prec`` or
``dec_prec`` for higher precision intervals).  This method works even when
not all tetrahedra are positively oriented::

   sage: M = Manifold("m015(3, 1)")
   sage: M.tetrahedra_shapes('rect', intervals=True)
   [0.625222762246? + 3.177940133813?*I,
    -0.0075523593782? + 0.5131157955971?*I,
    0.6515818912107? - 0.1955023488930?*I]

This is all achieved through a reimplementation of `HIKMOT
<http://www.oishi.info.waseda.ac.jp/~takayasu/hikmot/>`_ which
pioneered the use of interval methods for hyperbolic manifolds.  It
can be used in a very similar way than HIKMOT, but uses Sage's complex
interval types and the Newton interval method (instead of the Krawczyk
test) for certification. See `Zgliczynski's notes`_ for a quick
overview of these two tests.

This verification code was contributed by Matthias Goerner.  


Methods of verifying hyperbolicity
----------------------------------

..   autofunction:: snappy.verify.verify_hyperbolicity

Generating certified shape intervals
------------------------------------

The recommeded way to obtain certified intervals for the shapes is via
``manifold.tetrahedra_shapes(intervals=True)`` as `described earlier
<verify.html>`_. Here we document the ``CertifiedShapesEngine`` used
internally to generate these intervals. It is of interest for those
users who want to understand the underlying interval math and
experiment with the Newton interval method.


..   automodule:: snappy.verify
..   autoclass:: CertifiedShapesEngine
     :members:
     :inherited-members:

