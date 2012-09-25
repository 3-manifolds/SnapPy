import matrix
from ptolemyObstructionClass import PtolemyObstructionClass
from ptolemyVariety import PtolemyVariety

def get_ptolemy_obstruction_classes(manifold):

    """
    Generates a list of obstruction cocycles representing each class in
    H^2(M,bd M; Z/2) suitable as argument for get_ptolemy_variety.
    The first element in the list is always the trivial obstruction class.

    See Definition 1.7 of
    Garoufalidis, Thurston, Zickert
    The Complex Volume of SL(n,C)-Representations of 3-Manifolds
    http://arxiv.org/abs/1111.2828

    s_f_t takes values +/-1 and is the value of evaluating the cocycle on
    face f of tetrahedron t.

    === Examples ===

    Get the obstruction classes for 4_1:
    
    >>> from snappy import Manifold
    >>> M = Manifold("4_1")
    >>> c = get_ptolemy_obstruction_classes(M)

    There are two such clases for 4_1:
    
    >>> len(c)
    2
    
    Print the non-trivial obstruction class:

    >>> c[1]
    PtolemyObstructionClass(s_0_0 + 1, s_1_0 - 1, s_2_0 - 1, s_3_0 + 1, s_0_0 - s_0_1, s_1_0 - s_3_1, s_2_0 - s_2_1, s_3_0 - s_1_1)

    Construct Ptolemy variety for non-trivial obstruction class:
    
    >>> p = get_ptolemy_variety(M, N = 2, obstruction_class = c[1])

    Short cut for the above code:
    
    >>> p = get_ptolemy_variety(M, N = 2, obstruction_class = 1)

    Obstruction class only makes sense for even N:
    
    >>> p = get_ptolemy_variety(M, N = 3, obstruction_class = c[1])
    Traceback (most recent call last):
        ...
    AssertionError: PtolemyObstructionClass only makes sense for even N


    Hence, we get only one variety if we ask for all obstruction classes:
    
    >>> len(get_ptolemy_variety(M, N = 3, obstruction_class = 'all'))
    1
    """

    # compute boundary maps for cellular chain complex
    chain_d3, dummy_rows, dummy_columns = (
        manifold._ptolemy_equations_boundary_map_3())
    chain_d2, dummy_rows, explain_columns = (
        manifold._ptolemy_equations_boundary_map_2())

    # get which faces are identified to face classes
    identified_face_classes = (
        manifold._ptolemy_equations_identified_face_classes())

    # transpose to compute cohomology
    cochain_d2 = matrix.matrix_transpose(chain_d3)
    cochain_d1 = matrix.matrix_transpose(chain_d2)

    # two consecutive maps in a chain complex should give zero
    assert matrix.is_matrix_zero(
        matrix.matrix_mult(cochain_d2, cochain_d1))

    # Change the basis of the groups in the cellular cochain complex by
    # the matrices basechangeN
    # the cellular cochain maps with respect to the new basis will be stored
    # in transformed_d2 and transformed_d1 and will have at most one
    # non-entry per row, respectively, per column.
    basechange3, basechange2, basechange1, transformed_d2, transformed_d1 = (
        matrix.simultaneous_smith_normal_form(cochain_d2, cochain_d1))

    # Perform consistency checks
    matrix.test_simultaneous_smith_normal_form(
        cochain_d2, cochain_d1, 
        basechange3, basechange2, basechange1,
        transformed_d2, transformed_d1)

    # Take the matrices modulo 2 because we want Z/2 coefficients
    transformed_d2 = matrix.matrix_modulo(transformed_d2, 2)
    transformed_d1 = matrix.matrix_modulo(transformed_d1, 2)

    # Perform consistency check
    assert (
        matrix.num_cols(transformed_d2) ==
        matrix.num_rows(transformed_d1))

    # We want to find a representative for each cohomology class
    
    # For this we find a subset of basis vectors such that they span a
    # subspace in C^2 such that each cohomology class has a unique
    # representative.

    # We find all these basis vectors by requiring that it is in the
    # kernel of d2 and that it is not in the image of d1.

    # In other words, we can always change a representative of a cohomology
    # class to take zero values on those entries of a vector which are hit
    # by the image of d1.

    is_generator_of_H2 = [ matrix.col_is_zero(transformed_d2, i) and 
                           matrix.row_is_zero(transformed_d1, i) 
                           for i in range(
                                   matrix.num_rows(transformed_d1))]

    # Compute the subspace of C^2 that is spanned by all the above basis
    # vectors.

    H2_elements_in_new_basis = _enumerate_all_binary_vectors(
        is_generator_of_H2)

    # And package it into obstruction class object we can return

    def construct_obstruction_class(index, H2_element_in_new_basis):

        # change back to the old basis
        H2_element = matrix.vector_modulo(
            matrix.matrix_mult_vector(
                basechange2,
                H2_element_in_new_basis),
            2)

        # convert to python list
        H2_element = [x for x in H2_element]

        # make sure it actually is in the kernel

        assert matrix.is_vector_zero(
            matrix.vector_modulo(
                matrix.matrix_mult_vector(cochain_d2, H2_element),
                2))

        # package it into a proper python object
        return PtolemyObstructionClass(manifold, index,
                                       H2_element, explain_columns,
                                       identified_face_classes)

    return [construct_obstruction_class(index, H2_element_in_new_basis)
            for index, H2_element_in_new_basis
            in enumerate(H2_elements_in_new_basis)]

def get_ptolemy_variety(manifold, N, obstruction_class = None,
                        simplify = True):

    """
    Generates Ptolemy variety as described in
    (1) Garoufalidis, Thurston, Zickert
    The Complex Volume of SL(n,C)-Representations of 3-Manifolds
    http://arxiv.org/abs/1111.2828
    
    (2) Garoufalidis, Goerner, Zickert:
    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
    http://arxiv.org/abs/1207.6711

    The variety can be exported to magma or sage and solved there. The
    solutions can be processed to compute invariants. See below.

    === Arguments ===

    N --- which SL(N,C) we want the variety.
    
    obstruction_class --- class from Definiton 1.7 of (1).
    None for trivial class or a value returned from get_ptolemy_obstruction_classes.
    Short cuts: obstruction_class = 'all' returns a list of Ptolemy varieties
    for each obstruction. For easier iteration, can set obstruction_class to 
    an integer.

    simplify --- boolean to indicate whether to simplify the equations which
    significantly reduces the number of variables.
    Simplifying means that several identified Ptolemy coordinates x = y = z = ...
    are eliminated instead of adding relations x - y = 0, y - z = 0, ...

    === Examples for 4_1 ===
    
    >>> from snappy import Manifold
    >>> M = Manifold("4_1")

    Get the varieties for all obstruction classes at once (use
    help(varieties[0]) for more information):
    
    >>> varieties = get_ptolemy_variety(M, N = 2, obstruction_class = "all")

    Print the variety as an ideal (sage object) for the non-trivial class:

    >>> varieties[1].ideal    #doctest: +SKIP                                                                        
    Ideal (c_0101_0^2 - c_0101_0 + 1, -c_0101_0^2 + c_0101_0 - 1, t*c_0101_0 - 1) of Multivariate Polynomial \
Ring in t, c_0101_0 over Rational Field                                                                       
    (skip doctest because example only works in sage and not plain python)

    >>> for eqn in varieties[1].equations:
    ...     print "    ", eqn
         1 - c_0101_0 + c_0101_0^2
         - 1 + c_0101_0 - c_0101_0^2

    Generate a magma file to compute Groebner basis for N = 3:
    
    >>> p = get_ptolemy_variety(M, N = 3)
    >>> print p.to_magma()          #doctest: +ELLIPSIS
    P<t, c_0012_1, c_0102_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0> := PolynomialRing(RationalField(), 7);
    I := ideal<P |
    c_0102_0 - c_0102_0 * c_1011_0 + c_1101_0,
        ...

    === If you have a magma installation ===

    Call p.compute_solutions() to automatically call magma on the above output
    and produce exact solutions!!!
    
    >>> try:
    ...     sols = p.compute_solutions(verbose)
    ... except:
    ...     sols = None     # magma failed, use precomputed output instead

    === If you do not have a magma installation ===

    Load a precomputed example from magma which is provided with the package:
    
    >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma, triangulation_from_magma
    >>> print _magma_output_for_4_1__sl3      #doctest: +ELLIPSIS
    <BLANKLINE>
    ==TRIANGULATION=BEGINS==
    % Triangulation
    4_1
    geometric_solution  2.02988321
    oriented_manifold
    ...

    Recover the original trigangulation:
    >>> M = triangulation_from_magma(_magma_output_for_4_1__sl3)
    >>> M.is_isometric_to(Manifold("4_1"))
    True

    Parse the file and produce solutions:

    >>> if sols is None:    # calling magma failed, so use precomputed example
    ...     sols = solutions_from_magma(_magma_output_for_4_1__sl3)

    === Continue here whether you have or do not have magma ===

    Pick the first solution of the three different solutions (up to Galois
    conjugates):
    
    >>> len(sols)
    3
    >>> solution = sols[0]

    Read the exact value for c_1020_0 (help(solution) for more information
    on how to compute cross ratios, volumes and other invariants):

    >>> solution['c_1020_0']
    Mod(1/2*x - 1, x^2 - x + 2)

    Example of simplified vs non-simplified variety:

    >>> simplified = get_ptolemy_variety(M, N = 4, obstruction_class = 1)
    >>> full = get_ptolemy_variety(M, N = 4, obstruction_class = 1, simplify = False)
    >>> len(simplified.variables), len(full.variables)
    (17, 70)
    >>> len(simplified.equations), len(full.equations)
    (20, 79)
    """
    
    if obstruction_class == 'all' and N % 2 == 1:
        return [ PtolemyVariety(manifold, N, simplify = simplify) ]

    if not (obstruction_class is None or 
            isinstance(obstruction_class, PtolemyObstructionClass)):
        
        obstruction_classes = get_ptolemy_obstruction_classes(manifold)

        if obstruction_class == 'all':
            return [PtolemyVariety(manifold, N, obstruction_class,
                                   simplify = simplify)
                    for obstruction_class in obstruction_classes]
        
        try:
            obstruction_class = obstruction_classes[obstruction_class]
        except:
            raise Exception("Bad index for obstruction class")

    return PtolemyVariety(manifold, N, obstruction_class,
                          simplify = simplify)

def _enumerate_all_binary_vectors(mask):
    if len(mask) == 0:
        yield [ ]
    else:
        for i in range(2):
            if mask[0] or i == 0:
                for j in _enumerate_all_binary_vectors(mask[1:]):
                    yield [ i ] + j
