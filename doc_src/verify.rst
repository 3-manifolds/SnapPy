Verified computations
========================================

When used inside `Sage <http://sagemath.org>`_, SnapPy can produce intervals for orientable manifolds that are guarenteed to contain a true solution to the rectangular gluing equations (give bits_prec or dec_prec for higher precision intervals):

   >>> M = Manifold("m015")
   >>> M.tetrahedra_shapes('rect', intervals=True)
   [0.6623589786224? + 0.5622795120623?*I,
    0.6623589786224? + 0.5622795120623?*I,
    0.6623589786224? + 0.5622795120623?*I]

Using interval arithmetics, it can also rigorously verify the hyperbolicity of an orientable manifold:

   >>> from snappy.verify import verify_hyperbolicty
   >>> M = Manifold("m015")
   >>> verify_hyperbolicty(M)
   (True,
    [0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I])   

This is achieved through a reimplementation of `HIKMOT <http://www.oishi.info.waseda.ac.jp/~takayasu/hikmot/>`_ which pioneered the use of interval methods for hyperbolic manifolds. 
It can be used in a very similar way than HIKMOT, but uses `Sage
<http://sagemath.org>`_'s complex interval types and the Newton
interval method (instead of the Krawczyk test) for certification. See
`this note <http://ww2.ii.uj.edu.pl/~zgliczyn/cap07/krawczyk.pdf>`_
for a quick overview of these two tests.

This verification code was contributed by Matthias Goerner.  


Methods of verifying hyperbolicity
==================================

..   autofunction:: snappy.verify.verify_hyperbolicity
..   autofunction:: snappy.verify.verify_logarithmic_gluing_equations_and_positively_oriented_tets

Generating certified shape intervals
=================================================

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

