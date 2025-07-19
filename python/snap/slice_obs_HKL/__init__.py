"""
Using metabelian representations to obstruct slicing
====================================================

Based on::

  Herald, Kirk, and Livingston, Math Zeit., 2010
  https://dx.doi.org/10.1007/s00209-009-0548-1
  https://arXiv.org/abs/0804.1355

and::

  Dunfield and Gong, "Ribbon concordances and slice obstructions:
  experiments and examples", https://arXiv.org/FILLIN

In in the basic case, the implementation follows the [HKL] paper closely, with a
few minor changes:

1. We use the default simplified presentation for pi_1(knot exterior)
   that SnapPy provides, rather than a Wirtinger presentation, as the
   former has many fewer generators (but of course much longer
   relators).

2. The construction of the metabelian representations in Section 7 of
   [HKL] is done the language of twisted cohomology and specifically
   twisted cocycles, whereas this viewpoint is not quite explicit
   (though clearly implied) in [HKL].

3. To match the conventions of the other twisted Alexander polynomials
   computed by SnapPy, we insist on all actions are on the left,
   rather than the mix of left and right actions used in [HKL].

4. For the final check of whether the twisted polynomial is a norm, we
   simply use that Sage can factor univariate polynomials over
   cyclotomic fields into irreducibles, rather than the various
   generalizations of Gauss's Lemma in [HKL].

Additionally, we implement the refinements and generalizations
introduced in [DG].
"""

from ... import SnapPy
from ...sage_helper import _within_sage, sage_method
if _within_sage:
    from ...sage_helper import ZZ, prime_range, prime_powers
from . import basics, rep_theory, direct


def expand_prime_power_spec(spec):
    if spec in ZZ:
        a, b = 0, spec
    else:
        if len(spec) != 2:
            raise ValueError(f'Spec {spec} does not specify a range')
        a, b = spec
    return prime_powers(a, b + 1)


def expand_prime_spec(spec):
    if spec in ZZ:
        a, b = 0, spec
    else:
        if len(spec) != 2:
            raise ValueError(f'Spec {spec} does not specify a range')
        a, b = spec
    return prime_range(a, b + 1)


@sage_method
def slice_obstruction_HKL(self,
                          primes_spec,
                          verbose=0,
                          check_in_S3=True,
                          method='advanced',
                          ribbon_mode=False):
    """

    For the exterior of a knot in the 3-sphere, search for a
    Herald-Kirk-Livingston (HKL) topological slice obstruction as
    described in:

    * [DG] Dunfield and Gong, Section 3 of https://arXiv.org/abs/FILLIN

    Specifically, it looks at the cyclic branched covers of the knot
    of prime-power order p and the F_q homology thereof where q is
    prime. The range of such (p, q) pairs searched is given by
    primes_spec as a list of (p_max, [q_min, q_max]).  It returns the
    pair (p, q) of the first nonzero obstruction found (in which case
    K is not slice), and otherwise returns ``None``::

       sage: M = Manifold('K12n813')
       sage: spec = [(10, [0, 20]), (20, [0, 10])]
       sage: M.slice_obstruction_HKL(spec, method='basic', verbose=1)
           Looking at (2, 3) ...
           Looking at (3, 2) ...
           Looking at (3, 7) ...
       (3, 7)

    You can also specify the p to examine by a range [p_min, p_max] or
    the q by just q_max::

       sage: spec = [([5, 10], 10)]
       sage: M.slice_obstruction_HKL(spec, method='advanced', verbose=1)
           Looking at (8, 3) ...
       (8, 3)

    If primes_spec is just a pair (p, q) then only that obstruction is
    checked::

       sage: M.slice_obstruction_HKL((3, 7))
       (3, 7)

    The ``method`` argument determines which HKL tests are employed:

    * ``method='basic'`` employs Lemma 3.5 of [DG] and corresponds
      roughly to the test used in SnapPy 3.1 and 3.2, though earlier
      versions required p to be prime and forbid q = 2.

    * ``method='advanced'`` is the default and employs Theorem 3.8
      from [DG]. This subsumes the `basic` method, but can be quite
      slow when the multiplicity of an irreducible V_i is large.
      Example::

        sage: M.slice_obstruction_HKL((2, 3), method='basic')  # returns None
        sage: M.slice_obstruction_HKL((2, 3), method='advanced')
        (2, 3)

    * ``method='direct'`` employs the method of Section 3.20 of [DG] to
      consider not just the F_q homology of the cover, but
      epimorphisms of its group to Z/q^e Z for e > 1.  This is the
      most computationally expensive method.  Example::

        sage: M = Manifold('K14n16945')
        sage: M.slice_obstruction_HKL((2, 3), method='advanced') # returns None
        sage: M.slice_obstruction_HKL((2, 3), method='direct')
        (2, 3)

    For any method, q = 2 is handled using Section 3.11 of [DG].

    If ``verbose`` is ``1`` or ``True``, it prints each pair (p, q) being considered;
    when ``verbose==2`` more is printed about each step.

    ADDD RIBBON MODE
    """
    if method not in ['basic', 'advanced', 'direct']:
        raise ValueError("Argument method is not 'basic', 'advanced', or 'direct'")

    M = self

    if M.cusp_info('is_complete') != [True]:
        raise ValueError('Need exactly one cusp which should be unfilled')
    if M.homology().elementary_divisors() != [0]:
        raise ValueError('Not the exterior of knot in S^3 as H_1 != Z')
    if check_in_S3:
        T = SnapPy.Triangulation(M)
        T.dehn_fill((1, 0))
        if T.fundamental_group().num_generators() != 0:
            raise ValueError('The (1, 0) filling is not obviously S^3')

    # Special case of only one (p, q) to check
    if len(primes_spec) == 2 and primes_spec[0] in ZZ and primes_spec[1] in ZZ:
        primes_spec = [([primes_spec[0]], [primes_spec[1]])]
    else:
        primes_spec = [(expand_prime_power_spec(a), expand_prime_spec(b))
                       for a, b in primes_spec]

    for ps, qs in primes_spec:
        for p in ps:
            d = basics.nonzero_divisor_product(M, p)
            for q in qs:
                if d % q == 0:
                    if verbose:
                        print('    Looking at', (p, q), '...')
                    if method=='direct':
                        success = direct.slicing_obstructed_by_larger_quotient(M, p, q, verbose)
                    else:
                        success = rep_theory.slicing_is_obstructed(M, p, q,
                                                                   skip_higher_mult=(method=='basic'),
                                                                   ribbon_mode=ribbon_mode,
                                                                   verbose=(verbose > 1))
                    if success:
                        return (p, q)
