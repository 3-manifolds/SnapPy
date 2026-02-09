"""
IMPORTANT: Python only recognises this as a doc string if there is
nothing before it. In particular, add any includes after the doc string.

>>> M = Manifold("m125(3,4)(0,0)")
>>> spec = M.length_spectrum_alt_gen()
>>> next(spec) # doctest: +NUMERIC9
Length                                      Core curve  Word
0.24208261435543 - 1.73621300277325*I       Cusp 0      aBDcDcb
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

Examples with +SKIP from length_spectrum_alt

>>> M = Manifold("m202(3,4)(3,4)")
>>> spec = M.length_spectrum_alt(count = 3)
>>> len(spec)
4
>>> spec[0].length # doctest: +NUMERIC9
0.14820741547094 - 1.76955170166922*I
>>> spec[1].length # doctest: +NUMERIC9
0.14820741547097 - 1.76955170166923*I
>>> spec[2].length # doctest: +NUMERIC9
0.79356651781096 + 2.65902431489655*I
>>> spec[3].length # doctest: +NUMERIC9
0.79356651781096 + 2.65902431489655*I

Verified:

sage: spec = M.length_spectrum_alt(count = 3, verified = True, bits_prec = 110)
sage: len(spec)
4
sage: spec[0].length # doctest: +NUMERIC9
0.14820741547094772? - 1.76955170166923543?*I
sage: spec[1].length # doctest: +NUMERIC9
0.14820741547094772? - 1.76955170166923543?*I
sage: spec[2].length # doctest: +NUMERIC9
0.79356651781095741? + 2.65902431489655135?*I
sage: spec[3].length # doctest: +NUMERIC9
0.79356651781095741? + 2.65902431489655135?*I

>>> M = Manifold("m202(3,4)(0,0)")
>>> M.length_spectrum_alt(max_len = 1.1) # doctest: +NUMERIC9
[Length                                      Core curve  Word
 0.14742465268512 - 1.78287093565202*I       Cusp 0      aabcDabcB,
 0.81161414965958 + 2.72911699294426*I       -           b,
 0.84163270359334 + 2.61245944742151*I       -           aB,
 0.93461379591349 + 2.70060614107722*I       -           a]

sage: M.length_spectrum_alt(max_len = 1.1, verified=True, bits_prec=130) # doctest: +NORMALIZE_WHITESPACE
[Length                                      Core curve  Word
 0.14742465268515... - 1.78287093565201...*I Cusp 0      aabcDabcB,
 0.81161414965958... + 2.72911699294425...*I -           b,
 0.84163270359334... + 2.61245944742151...*I -           aB,
 0.93461379591349... + 2.70060614107721...*I -           a]


"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")
