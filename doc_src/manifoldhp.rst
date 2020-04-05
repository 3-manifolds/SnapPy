ManifoldHP: High-precision variant
==================================================

A ManifoldHP is a variant of :class:`Manifold <snappy.Manifold>` which
does all floating-point calculations in `quad-double precision
<http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf>`_,
which has four times as many significant digits as the ordinary
`double precision numbers
<http://en.wikipedia.org/wiki/Double_precision_floating-point_format>`_
used by Manifold.  More precisely, numbers used in ManifoldHP have 212
bits for the mantissa/significand (roughly 63 decimal digits) versus
53 bits with Manifold.

To the user, the only difference between Manifold and ManifoldHP is the extra precision::

   >>> L = Manifold('m004')
   >>> L.volume()
   2.02988321282
   >>> H = ManifoldHP('m004')
   >>> H.volume()
   2.029883212819307250042405108549040571883378615060599584034978214

and it is easy to go back and forth between the two types::

    >>> D = H.low_precision()
    >>> D.volume(), type(D)
    (2.02988321282, <class 'snappy.Manifold'>)
    >>> U = L.high_precision()
    >>> type(U)
    <class 'snappy.ManifoldHP'>

FAQ
---

Q. How does this differ from the program `Snap <http://snap-pari.sourceforge.net/>`_ or the :doc:`corresponding features <snap>` of SnapPy? 

A. Snap computes hyperbolic structures to whatever precision you specify, not just 212 bits.  However, only some aspects of that structure can be accessed at the higher precision.  In contrast, with ManifoldHP every part of the SnapPea kernel uses the high-precision structure.  Eventually, we hope to add a ManifoldAP which allows for arbitrary precision throughout the kernel.  

Q. Are there any negatives to using ManifoldHP over Manifold?

A. Yes, ManifoldHP is generally slower by a factor of 10 to 100.  Multiplying two quad-double numbers requires at least 10 ordinary double multiplications, so some of this is inevitable.  

Q. What is one place where the extra precision really helps?  

A. Computing Dirichlet domains and subsidiary things like the length spectrum. A ManifoldHP can find the Dirichlet domain of a typically 15 crossing knot exterior but Manifold can't.  

