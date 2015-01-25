hikmot2: verified computations of shapes
========================================

When used inside `Sage <http://sagemath.org>`_, SnapPy can produce intervals for orientable manifolds that are guarenteed to contain a true solution to the rectangular gluing equations (give bits_prec or dec_prec for higher precision intervals):

   >>> M = Manifold("m015")
   >>> M.tetrahedra_shapes('rect', intervals=True)
   [0.6623589786224? + 0.5622795120623?*I,
    0.6623589786224? + 0.5622795120623?*I,
    0.6623589786224? + 0.5622795120623?*I]

Using interval arithmetics, it can also rigorously verify the hyperbolicity of an orientable manifold:

   >>> from snappy.hikmot2 import verify_hyperbolicty
   >>> M = Manifold("m015")
   >>> verify_hyperbolicty(M)
   (True,
    [0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I,
     0.6623589786224? + 0.5622795120623?*I])   

This is achieved through hikmot2, a reimplementation of `hikmot <http://www.oishi.info.waseda.ac.jp/~takayasu/hikmot/>`_ which pioneered the use of interval methods for hyperbolic manifolds. 
It can be used in a very similar way than hikmot, but uses `Sage <http://sagemath.org>`_'s complex interval types and the Newton interval method (instead of the Krawczyk test) for certification.

.. toctree::
   
   hikmot2Methods
   hikmot2CertifiedShapesEngine
