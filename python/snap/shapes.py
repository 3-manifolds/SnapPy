from ..sage_helper import _within_sage
from flint import arb, acb, fmpz, fmpz_mat, acb_mat
#from ..pari import pari, Gen
from snappy.number import (
    Number,
    bit_precision,
    prec_dec_to_bits,
    prec_bits_to_dec)

if _within_sage:
    from ..sage_helper import ComplexField

#def is_pari(x) -> bool:
#    return isinstance(x, Gen)

# A gluing equation has the form:
# z_0^a_0 * (1-z_0)^b_0 * ... z_N^a_N*(1-z_N)^b_N = +/-1
# A gluing euation is represented as a triple a, b, c where a and b are lists
# of the exponents of z_i and c = +/-1 is the right hand side.

def eval_gluing_equation(eqn, shapes):
    """
    Evaluate the product of cross ratios in an edge equation.
    The result will be 1 if the equation is satisfied.
    """
#    if is_pari(eqn):
#        raise RuntimeError('eval_gluing_equation: got pari eqn')
    a, b, c = eqn
    ans = int(c)
    for i, z in enumerate(shapes):
        ans = ans * (z**int(a[i]) * (1 - z)**int(b[i]))
    return ans


def gluing_equation_errors(eqns, shapes):
    """ A list containing the difference from 1 for each equation."""
    return [eval_gluing_equation(eqn, shapes) - 1 for eqn in eqns]


def infinity_norm(L):
#    if is_pari(L):
#        raise RuntimeError('infinity_norm: got pari vector')
    result = max(map(abs, L))
    return result


def gluing_equation_error(eqns, shapes):
    """ The max of the equation errors."""
    return infinity_norm(gluing_equation_errors(eqns, shapes))


def XXXenough_gluing_equations(manifold):
    """
    Row reduce the gluing equations and select a full-rank subsystem.
    """
    n_tet = manifold.num_tetrahedra()
    n_cusps = manifold.num_cusps()
    eqns = manifold.gluing_equations("rect")
    # The first n_tet equations are the edge equations.
    edge_eqns = fmpz_mat([a + b for a, b, _ in eqns[:n_tet]])
    RHS = fmpz_mat([[(1 - c) // 2] for _, _, c in eqns[:n_tet]])
    # Row reduce the equation matrix: U * edge_eqns = H
    H, U = edge_eqns.hnf(transform=True)
    # Remove the zero rows
    reduced_eqns = [row for row in H.tolist() if any(row)]
    # Check that we have the expected rank:
    assert  len(reduced_eqns) == n_tet - n_cusps
    # Apply the same row ops to the right hand side
    RHS = U * RHS
    # Convert RHS to a list of ints
    RHS = [( -1 if RHS[n,0] // 2 else 1) for n in range(len(reduced_eqns))]
    # Express the row reduced equations as triples. 
    edge_eqns_with_RHS = [(e[:n_tet], e[n_tet: 2 * n_tet], c)
                          for e, c in zip(reduced_eqns, RHS)]
    # Add triples for the row equations
    cusp_eqns = []
    j = n_tet
    for i in range(n_cusps):
        cusp_eqns.append(eqns[j])
        j += 2 if manifold.cusp_info(i)['complete?'] else 1
    ans_eqns = edge_eqns_with_RHS + cusp_eqns
    assert len(ans_eqns) == n_tet
    ans_matrix = fmpz_mat([a + b for a, b, _ in ans_eqns])
    assert ans_matrix.rank() == ans_matrix.nrows()
    # Return the equations as lists of python ints
    return [(list(map(int, A)), list(map(int, B)), int(c)) for A, B, c in ans_eqns]

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
    # Perform the same row ops on the auugmented matrix
    edge_eqns_with_RHS = U * edge_eqns_with_RHS
    edge_eqns_with_RHS = [(e[:n_tet], e[n_tet: 2 * n_tet], (-1)**int(e[-1]))
                          for e in edge_eqns_with_RHS.tolist()[:non_zero_rows]]
    # Add triples for the cusp equations
    cusp_eqns = []
    j = n_tet
    for i in range(n_cusps):
        cusp_eqns.append(eqns[j])
        j += 2 if manifold.cusp_info(i)['complete?'] else 1
    ans_eqns = edge_eqns_with_RHS + cusp_eqns
    assert len(ans_eqns) == n_tet
    ans_matrix = fmpz_mat([a + b for a, b, _ in ans_eqns])
    assert ans_matrix.rank() == ans_matrix.nrows()
    # Return the equations as lists of python ints
    return [(list(map(int, A)), list(map(int, B)), int(c)) for A, B, c in ans_eqns]


def polished_tetrahedra_shapes(manifold, dec_prec=None, bits_prec=200, ignore_solution_type=False):
    """
    Refines the current solution to the gluing equations to one with
    the specified accuracy.
    """
    if dec_prec is None:
        dec_prec = prec_bits_to_dec(bits_prec)
    else:
        bits_prec = prec_dec_to_bits(dec_prec)
    working_prec = 2 * bits_prec + 30

    # This is a potentially long calculation, so we cache the result
    # and use the cached values as a start, if possible.
    if "polished_shapes" in manifold._cache:
        curr_bits_prec, curr_sol = manifold._cache["polished_shapes"]
        if bits_prec <= curr_bits_prec:
            return curr_sol

    # Check to make sure the initial solution is reasonable
    if not ignore_solution_type and manifold.solution_type() not in [
            'all tetrahedra positively oriented',
            'contains negatively oriented tetrahedra']:
        raise ValueError('Initial solution to gluing equations has '
                         'flat or degenerate tetrahedra')
    initial_shapes = [z.flint_obj for z in manifold.tetrahedra_shapes('rect')]
    original_equations = manifold.gluing_equations('rect')
    if gluing_equation_error(original_equations, initial_shapes) > 0.000001:
        raise ValueError('Initial solution not very good')

    # Now begin the actual computation
    eqns = enough_gluing_equations(manifold)
    n_rows = len(eqns)
    n_cols = len(initial_shapes)
    with bit_precision(working_prec):
        target_epsilon = arb(2.0) ** -bits_prec
        initial_shapes = acb_mat(n_cols, 1, initial_shapes)
        initial_error = infinity_norm(gluing_equation_errors(eqns, initial_shapes))
        shapes = acb_mat(n_cols, 1, initial_shapes)
        for i in range(100):
            errors = gluing_equation_errors(eqns, shapes)
            error = infinity_norm(errors)
            if error < target_epsilon:
                break
            if error > 100 * initial_error:
                print('Diverging at step', i, 'with error', error)
                break
            derivative = acb_mat(
                [[eqn[0][i] / z - eqn[1][i] / (1 - z) for i, z in enumerate(shapes)]
                     for eqn in eqns])
            rhs = acb_mat(n_rows, 1, errors)
            correction = derivative.solve(rhs)
            # Replace each acb shape by its midpoint for the next iteration.
            shapes = acb_mat(n_rows, 1, [z.mid() for z in shapes - correction])

        # Check that things worked out ok.
        error = gluing_equation_error(original_equations, shapes)
        total_change = infinity_norm(shapes - initial_shapes)
        if error > 1000 * target_epsilon or total_change > 0.0000001:
            raise ValueError('Did not find a good solution to the gluing equations')
        result = [Number(s, precision=bits_prec) for s in shapes]
        if _within_sage:
            CC = ComplexField(bits_prec)
            result = [CC(repr(z)) for z in result]
        manifold._cache["polished_shapes"] = (bits_prec, result)
        return result
