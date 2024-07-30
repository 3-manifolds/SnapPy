"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

>>> M = Manifold("m125(3,4)(0,0)")
>>> spec = M.length_spectrum_alt_gen()
>>> next(spec) # doctest: +NUMERIC9
Length                                       Word          Core curve
0.24208261435543 - 1.73621300277325*I        aBDcDcb       Cusp 0
>>> next(spec).length # doctest: +NUMERIC9
0.90986906036840 + 3.03574280072295*I
>>> next(spec).length # doctest: +NUMERIC9
0.98258854739348 - 2.33353878259198*I
>>> next(spec).length # doctest: +NUMERIC9
1.00610499287709 - 3.02617893116978*I
>>> next(spec).length # doctest: +NUMERIC9
1.13551437663552 - 2.12918861416187*I
>>> next(spec).length # doctest: +NUMERIC9
1.63203771292969 + 2.30009520293758*I

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")
