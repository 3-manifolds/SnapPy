from ..sage_helper import _within_sage
import flint
from flint import arb, acb, fmpz, fmpz_mat, acb_mat

#from ..pari import pari, Gen
from snappy.number import (
    Number,
    bit_precision,
    prec_dec_to_bits,
    prec_bits_to_dec)

if _within_sage:
    from ..sage_helper import ComplexField

# A gluing equation has the form:
# z_0^a_0 * (1-z_0)^b_0 * ... z_N^a_N*(1-z_N)^b_N = +/-1
# A gluing euation is represented as a triple a, b, c where a and b are lists
# of the exponents of z_i and c = +/-1 is the right hand side.

# The following 3 functions take acb arguments.  You must set the precision
# before calling them!

def eval_gluing_equation(eqn, acb_shapes):
    """
    Evaluate the product of cross ratios in an edge equation.
    The result will be 1 if the equation is satisfied.
    """
    a, b, c = eqn
    ans = arb(c)
    for i, z in enumerate(acb_shapes):
        ans = ans * (z**int(a[i]) * (1 - z)**int(b[i]))
    return ans


def gluing_equation_errors(eqns, acb_shapes):
    """ A list containing the difference from 1 for each equation."""
    return [eval_gluing_equation(eqn, acb_shapes) - 1 for eqn in eqns]


def infinity_norm(L):
    result = max(map(abs, L))
    return result


def gluing_equation_error(eqns, shapes):
    """ The max of the equation errors."""
    return infinity_norm(gluing_equation_errors(eqns, shapes))


def enough_gluing_equations(manifold):
    """
    Row reduce the gluing equations and select a full-rank subsystem.
    """
    n_tet = manifold.num_tetrahedra()
    n_cusps = manifold.num_cusps()
    eqns = manifold.gluing_equations("rect")
    # The first n_tet equations are the edge equations.
    edge_eqns = fmpz_mat([a + b for a, b, _ in eqns[:n_tet]])
    edge_eqns_with_RHS = fmpz_mat([a + b + [(1 - c) // 2]
                                         for a, b, c in eqns[:n_tet]])
    # Row reduce the equation matrix to get: U * edge_eqns = H
    H, U = edge_eqns.hnf(transform=True)
    # Check that we have the expected rank:
    non_zero_rows = len([row for row in H.tolist() if any(row)])
    assert non_zero_rows == n_tet - n_cusps
    # Perform the same row ops on the auugmented matrix.
    edge_eqns_with_RHS = U * edge_eqns_with_RHS
    # Rewrite the equations as losts of triples.
    edge_eqns_with_RHS = [(e[:n_tet], e[n_tet: 2 * n_tet], (-1)**int(e[-1]))
                          for e in edge_eqns_with_RHS.tolist()[:non_zero_rows]]
    # Add the cusp equations.
    cusp_eqns = []
    j = n_tet
    for i in range(n_cusps):
        cusp_eqns.append(eqns[j])
        j += 2 if manifold.cusp_info(i)['complete?'] else 1
    ans_eqns = edge_eqns_with_RHS + cusp_eqns
    # Do sanity checks.
    assert len(ans_eqns) == n_tet
    ans_matrix = fmpz_mat([a + b for a, b, _ in ans_eqns])
    assert ans_matrix.rank() == ans_matrix.nrows()
    # Return the equations as lists of python ints.
    return [(list(map(int, A)), list(map(int, B)), int(c)) for A, B, c in ans_eqns]


def polished_tetrahedra_shapes(manifold, dec_prec=None, bits_prec=200,
                               certify_prec=None, ignore_solution_type=False):
    """
    Refine a solution to the gluing equations to the specified
    precision.
    """

    if dec_prec is None:
        dec_prec = prec_bits_to_dec(bits_prec)
    else:
        bits_prec = prec_dec_to_bits(dec_prec)

    max_working_prec = bits_prec + 30

    # This is a potentially long calculation, so we cache the result
    # and use the cached values as a start, if possible.
    if ('polished_shapes', bits_prec) in manifold._cache:
        return manifold._cache['polished_shapes', bits_prec]

    # Check to make sure the initial solution is reasonable
    if not ignore_solution_type and manifold.solution_type() not in [
            'all tetrahedra positively oriented',
            'contains negatively oriented tetrahedra']:
        raise ValueError('Initial solution to gluing equations has '
                         'flat or degenerate tetrahedra')
    initial_shapes = manifold.tetrahedra_shapes('rect')
    initial_shape_prec = initial_shapes[0].precision()

    original_equations = manifold.gluing_equations('rect')
    if gluing_equation_error(original_equations, initial_shapes) > 0.000001:
        raise ValueError('Initial solution not very good')

    # Now begin the actual computation.
    # This computation does not work with the algebraic shape
    # equations, but rather with the nearly linear system obtained by
    # setting the log of the left hand side to 0.
    eqns = enough_gluing_equations(manifold)
    n = len(eqns)

    original_precision = flint.ctx.prec
    flint.ctx.prec = min(2*initial_shape_prec, max_working_prec)

    target_epsilon = arb(2.0) ** -bits_prec
    initial_shapes = acb_mat(n, 1, [z.flint_obj for z in initial_shapes])
    initial_error = infinity_norm(gluing_equation_errors(eqns, initial_shapes))
    shapes = acb_mat(n, 1, initial_shapes.entries())
    for i in range(100):
        errors = gluing_equation_errors(eqns, shapes)
        error = infinity_norm(errors)
        if error < target_epsilon:
            break
        if error > 100 * initial_error:
            print('Diverging at step', i, 'with error', error)
            break

        # Note that these are logarithmic derivatives!
        derivative = acb_mat([[eqn[0][i] / z - eqn[1][i] / (1 - z)
                               for i, z in enumerate(shapes)]
                               for eqn in eqns])
        rhs = acb_mat(n, 1, errors)
        correction = derivative.solve(rhs, algorithm="approx")
        corrected_shapes = shapes - correction

        # When refining shapes, we don't want our intervals to expand
        # with each Newton iteration.  So we reset each shape to an exact
        # value before starting the next iteration.
        shapes = corrected_shapes.mid()

        # If shapes didn't move much, increase the working precision.
        if flint.ctx.prec < max_working_prec:
            if infinity_norm(correction) < arb(2.0) ** -int(0.8*flint.ctx.prec):
                flint.ctx.prec = min(2*flint.ctx.prec, max_working_prec)

    # Check that things worked out ok.
    error = gluing_equation_error(original_equations, shapes)
    total_change = infinity_norm(shapes - initial_shapes)
    flint.ctx.prec = original_precision
    if error > 1000 * target_epsilon or total_change > 0.0000001:
        raise ValueError('Did not find a good solution to the gluing equations')
    # Prepare the final result
    # We want to preserve the intervals, so use the corrected shapes here
    result = [Number(s, precision=bits_prec) for s in corrected_shapes]
    manifold._cache['polished_shapes', bits_prec] = result
    return result
