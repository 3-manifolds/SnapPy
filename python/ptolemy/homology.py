from . import matrix

def _gcd(s, t):
    if t == 0:
        return s
    return _gcd(t, s % t)

def homology_basis_representatives_with_orders(d1, d2, N):
    """
    Given two matrices d1 and d2 such that d1 * d2 = 0, computes the
    homology H_1 = ker(d1) / im(d2) when using Z (when N = 0)
    or Z/N coefficients. The result is a list of pairs (vector, order)
    forming a basis of H_1 as Abelian group.
    A pair (vector, order) consists of a vector c_1 in ker(d1) representing a
    homology class [c_1] in H_1 with the integer order being the order of [c_1]
    if it is finite otherwise 0.
    """

    # two consecutive maps in a chain complex should give zero    
    assert matrix.is_matrix_zero(
        matrix.matrix_mult(d1, d2))

    # Change the basis of the chain groups C_0, C_1, C_2 by the matrices
    # basechangeN
    # transformed_d1 and transformed_d2 are the boundary maps d1 and d2 but
    # with respect to the new basis. In this new basis, the boundary maps will
    # have at most one non-zero entry per row, respectively, per column so that
    # we can directly read of whether a basis vector is in ker(d1) and im(d2).
    basechange0, basechange1, basechange2, transformed_d1, transformed_d2 = (
        matrix.simultaneous_smith_normal_form(d1, d2))

    # Perform consistency checks
    matrix.test_simultaneous_smith_normal_form(
        d1, d2, 
        basechange0, basechange1, basechange2,
        transformed_d1, transformed_d2)

    # Perform consistency check
    assert matrix.is_matrix_zero(
        matrix.matrix_mult(transformed_d1, transformed_d2))

    # Will hold the result
    homology_basis = []

    # A list of all the basis vectors
    basis_vectors = matrix.matrix_transpose(basechange1)

    # Iterate through the basis vectors
    for i, basis_vector in enumerate(basis_vectors):

        # Get the absolute value of the one non-zero entry in
        # the i-th column of d1 with respect to the new basis
        d1_entry = matrix.max_abs_of_col(transformed_d1, i)

        # i-th row of d2
        d2_entry = matrix.max_abs_of_row(transformed_d2, i)

        # Note that d1 * d2 = 0, so at most one of d1_entry and d2_entry
        # can be non-zero.

        # If using Z/N coefficients, apply gcd and modulo N to entries.
        if not N == 0:
            # We do this because we want to compute the image and the kernel
            # of maps m: Z/N -> Z/N, x |-> e*x.
            # For this consider the quantity gcd(e,N) % N.
            # The image of m consists of multiplies of gcd(e,N) % N.
            # If gcd(e, N) % N is 0, the kernel is Z/N. If 1, the kernel is 0.
            # Otherwise, N / (gcd(e, N) % N) is generating the kernel and the
            # kernel has size gcd(e, N) % N).
            
            # Example: Z/6 -> Z/6, x |-> 5 * x has full image Z/6 because
            # gcd(5, 6) = 1. The kernel is 0.

            # Z/6 -> Z/6, x |-> 4 * x has image all even elements in Z/6
            # because gcd(4, 6) = 2. Its kernel also consists of 2 elements,
            # namely 0 and 3 = 6 / 2 = 6 / gcd(4, 6).

            d1_entry = _gcd(d1_entry, N) % N
            d2_entry = _gcd(d2_entry, N) % N

        if d1_entry == 0:
            # The image of the basis vector is 0, hence the basis vector
            # is in the kernel of d1.

            # However, it might be entirely killed by the image of d2.
            # It is not entirely killed if and only if d2 entry is not 1.
            if not d2_entry == 1:
                # The basis vector times the d2 entry is in the image of d2
                # and thus zero in the homology.
                # Thus the d2 entry determines the order of the homology class
                # this basis vector represents.
                if N == 0:
                    # In the Z-coefficient case, the order is just given by the
                    # d2 entry. This is even when d2 entry is zero as infinite
                    # order is encoded by 0.
                    homology_basis.append( (basis_vector, d2_entry) )
                else:
                    if d2_entry == 0:
                        # If the d2 entry is 0, no multiple of the basis vector
                        # is killed by the image. But we have Z/N-coefficients,
                        # so the order is N.
                        homology_basis.append( (basis_vector, N) )
                    else:
                        # Order is just d2 entry.
                        homology_basis.append( (basis_vector, d2_entry) )

        elif not d1_entry == 1:
            # d1 entry is neither 0 or 1. In Z-coefficients this means
            # that no multiple of basis vector can be in the kernel.
            if not N == 0:
                # But in Z/N-coefficients, multiplying by N/d1_entry
                # gives a generator of the kernel that has order d1_entry
                # in the homology.
                kernel_vector = [ ( b * N / d1_entry ) % N
                                  for b in basis_vector ]
                homology_basis.append( (kernel_vector, d1_entry) )

        # d1 entry is 1, so the basis vector or any multiple of it
        # is in the kernel. Do not add it
            
    return homology_basis

def homology_representatives(d1, d2, N):
    """
    Given two matrices d1 and d2 such that d1 * d2 = 0, computes the
    homology H_1 = ker(d1) / im(d2) when using Z/N coefficients.

    The result is a list of vectors c_1 for each homology class [c_1].
    """

    assert N > 1

    # Compute a basis of the homology group
    homology_basis = homology_basis_representatives_with_orders(d1, d2, N)

    # Enumerate all linear combinations of the basis elements
    return _enumerate_from_basis(homology_basis, matrix.num_cols(d1), N)

def _enumerate_from_basis(basis, l, N):
    """
    Given a list of pairs (v_i, o_i) where each v1 is a vector of length l and 
    an integer N > 2, iterate through all linear combinations 
        k_0 v_0 + k_1 v_1 + ... + k_m v_m (modulo N) where k_i = 0, ..., o_i.
    If basis is empty, just return zero vector of length l.
    """

    if len(basis) == 0:
        # Base case, return zero vector of length l
        yield l * [ 0 ]
    else:
        # Take first pair (vector, order) from given list
        basis_vector, order = basis[0]
        # Iterate up until order
        for i in range(order):
            # Iterate through all linear combinations of remaining basis
            # elements
            for vector in _enumerate_from_basis(basis[1:], l, N):
                # Return the linear combination of the multiple of this basis
                # vector and the linear combination of the other basis vectors
                yield [ (i * b + v) % N for b, v in zip(basis_vector, vector) ]
