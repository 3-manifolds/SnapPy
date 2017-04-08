.. Documentation of the Spherogram part of SnapPy

..   automodule:: spherogram
		     
Links: planar diagrams and invariants
=======================================



Tutorial
--------

SnapPy includes the `Spherogram
<https://bitbucket.org/t3m/spherogram>`_ module which allows one to
create links programmatically.  The graphical conventions used are
`summarized here
<https://bitbucket.org/t3m/spherogram/raw/tip/spherogram_src/links/doc.pdf>`_.

First, here is the figure-8 knot assembled manually from four crossings, with conventions similar to those used by `KnotTheory <http://katlas.org/wiki/Planar_Diagrams>`_::
       
       >>> a, b, c, d = [Crossing(x) for x in 'abcd']
       >>> a[0], a[1], a[2], a[3] = c[1], d[0], b[1], b[0]
       >>> b[2], b[3] = d[3], c[2]
       >>> c[3], c[0] = d[2], d[1]
       >>> L = Link([a,b,c,d])
       >>> E = L.exterior()
       >>> E.volume()
       2.029883212819
       >>> Manifold('4_1').is_isometric_to(E)
       True

We can also give the same knot as a rational tangle::

       >>> L = RationalTangle(3,5).denominator_closure()
       >>> L.PD_code()
       [[6, 3, 7, 4], [4, 2, 5, 1], [0, 6, 1, 5], [2, 7, 3, 0]]
       >>> L.DT_code(True)
       'DT[dadCDAB]'

The natural algebra of tangles `shown here
<https://bitbucket.org/t3m/spherogram/raw/tip/spherogram_src/links/doc.pdf>`_
all works.  For instance, we can build the (-2, 3, 7) pretzel knot by
adding together three rational tangles::
      
      >>> T = RationalTangle(-1, 2) + RationalTangle(1, 3) + RationalTangle(1, 7)
      >>> L = T.numerator_closure()
      >>> Manifold('m016').is_isometric_to(L.exterior())
      True
      
To create the figure-8 knot as a closed braid, we first mash tangles
together horizontally using "|" to make the standard braid generators;
then multiplication in the braid group is just tangle multiplication::

   >>> C, Id = RationalTangle(1), IdentityBraid(1)
   >>> x = sigma_1 = C | Id
   >>> y = sigma_2_inverse = Id | -C
   >>> L = (x*y*x*y).denominator_closure()
   >>> E = L.exterior()
   >>> Manifold('4_1').is_isometric_to(E)
   True

Here's the minimally-twisted five chain from Figure 2 of `this paper
<http://arxiv.org/abs/math.GT/0209214>`_::

      def twisted_chain(n, k):
           T = RationalTangle(1, 2)
	   m = (n+1)//2
	   base = (m*[T, -T])[:n]
           tangles = base + [RationalTangle(k)]
           return sum(tangles, RationalTangle(0) ).bridge_closure()

      >>> L = twisted_chain(5, -1)
      >>> L.exterior().volume()
      10.14941606410

Spherogram includes ways to create very large random links, see below.
When used inside `Sage <http://sagemath.org>`_, one can compute many
basic link invariants, including the Jones polynomial.  See the
complete list of Link methods below.  

  
Random Links
------------

.. autofunction:: spherogram.random_link


The Link class
--------------

.. autoclass:: spherogram.Link
   :members:
   :inherited-members:
   :undoc-members:

The ClosedBraid class
---------------------

The ClosedBraid class provides an alternative way to construct links
as closed braids.  It is a subclass of Link, and currently defines
the same methods and attributes.

.. autoclass:: spherogram.ClosedBraid
