.. Documentation of the Snap part of SnapPy

Number theory of hyperbolic 3-manifolds
=============================================

SnapPy has rudimentary support for arbitrary-precision computation and
for identifying number fields associated to hyperbolic
3-manifolds.  While this functionality is far less than that of `Snap
<http://snap-pari.sf.net/>`_, it will eventually improve and even now
may occasionally be useful.  Except for the first example, one
currently needs to use SnapPy inside of `Sage <http://sagemath.org>`_
to have access to these features.  Here's how to find the tetrahedra
shapes to high-precision::

       sage: import snappy
       sage: M = snappy.Manifold('m004')
       sage: M.tetrahedra_shapes('rect', bits_prec=100)
       [0.50000000000000000000000000000 + 0.86602540378443864676372317075*I, 0.50000000000000000000000000000 + 0.86602540378443864676372317075*I]
       sage: M.tetrahedra_shapes('rect', dec_prec=100)[0]
       0.500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +
       0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093378624288*I

One can also compute the holonomy representation to any precision::

    sage: import snappy.snap as snap
    sage: G = snap.polished_holonomy(M, bits_prec=100)
    sage: G.SL2C('a')
    [0.50000000000000000000000000000 + 0.86602540378443864676372317075*I                                   -1.0000000000000000000000000000*I]
    [                                   1.0000000000000000000000000000*I   1.0000000000000000000000000000 - 1.7320508075688772935274463415*I]

You can also try to guess the shapes exactly using an LLL-based
algorithm of the type pioneered by Snap::

          sage: T = snappy.snap.tetrahedra_field_gens(M)
	  sage: T.find_field(prec=100, degree=10, optimize=True)
	  (Number Field in z with defining polynomial x^2 - x + 1, <ApproxAN: 0.5 + 0.866025403784*I>, [x, x])

In more complicated examples, one needs to use higher precision (or
degree) to actually find the exact values::

	  sage: N = snappy.Manifold('m004(1,2)')
	  sage: K = snap.trace_field_gens(N)
	  sage: K
	  <SetOfAAN: [1.83341116859 - 0.676655426231*I, -2.13471864115 + 0.0209637751151*I, -2.13471864115 + 0.0209637751151*I]>
	  sage: K.find_field(prec=100, degree=10, optimize=True)    # Fails, so no output 
	  sage: K.find_field(prec=200, degree=10, optimize=True)
	  (Number Field in z with defining polynomial x^7 - 3*x^6 - x^5 + 5*x^4 + 3*x^3 - 3*x^2 - 2*x + 1, 
	   <ApproxAN: -0.723173028398 - 0.587151903177*I>, 
           [x^6 - 3*x^5 - x^4 + 5*x^3 + 3*x^2 - 3*x - 1, x^6 - 3*x^5 + 2*x^3 + 3*x^2 - 1, x^6 - 3*x^5 + 2*x^3 + 3*x^2 - 1])


We can also compute various hyperbolicly-twisted Alexander
polynomials, as described `here <http://dunfield.info/torsion>`_::

	sage: M = snappy.Manifold('5_2')
	sage: snap.alexander_polynomial(M)
	2*a^2 - 3*a + 2
	sage: snap.hyperbolic_torsion(M, bits_prec=100)
	(2.3376410213776269870195455729 - 0.56227951206230124389918214504*I)*a^2 
	- 4.0000000000000000000000000003*a 
	+ 2.3376410213776269870195455731 - 0.56227951206230124389918214477*I
	sage: snap.hyperbolic_SLN_torsion(M, 3, 100)   # Dubois-Yamagachi adjoint torsion
	(0.40431358073618481197132660504 +
	0.75939451500971650241038772223*I)*a^3 
	+ (2.9032849613891083021420278850 -
	4.1185388389935516999882632998*I)*a^2 
	+ (-2.9032849613891083021420278809 +
	4.1185388389935516999882633007*I)*a 
	- 0.40431358073618481197132661847 - 0.75939451500971650241038771418*I
	sage: snap.hyperbolic_SLN_torsion(M, 4, 100)   # Why not?
	(2.5890988184099251088892745185 + 3.5126610817613336586374292713*I)*a^4
	+ (10.357403823939297224437742077 - 13.378446302375451727042633120*I)*a^3
	+ (-26.821802363180149782221451472 + 7.0253221635226673172748587283*I)*a^2
	+ (10.357403823939297224437738856 - 13.378446302375451727042631346*I)*a 
	+ 2.5890988184099251088892549440 + 3.5126610817613336586374448040*I

