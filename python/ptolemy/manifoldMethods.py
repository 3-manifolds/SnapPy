from __future__ import print_function

from . import matrix
from . import homology
from .ptolemyObstructionClass import PtolemyObstructionClass
from .ptolemyGeneralizedObstructionClass import PtolemyGeneralizedObstructionClass
from .ptolemyVariety import PtolemyVariety
from .utilities import MethodMappingList

def _gcd(s, t):
    if t == 0:
        return s
    return _gcd(t, s % t)

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
    AssertionError: PtolemyObstructionClass only makes sense for even N, try PtolemyGeneralizedObstructionClass


    When specifying N = 3, it automatically uses generalized obstruction class.
    
    >>> len(get_ptolemy_variety(M, N = 3, obstruction_class = 'all'))
    2
    """

    # Compute the obstruction classes
    H2_elements, explain_columns = get_obstruction_classes(manifold, 2)

    # get which faces are identified to face classes
    identified_face_classes = (
        manifold._ptolemy_equations_identified_face_classes())

    # package into PtolemyObstructionClass objects
    return [PtolemyObstructionClass(manifold, index,
                                    H2_element, explain_columns,
                                    identified_face_classes)
            for index, H2_element
            in enumerate(H2_elements)]

def get_generalized_ptolemy_obstruction_classes(manifold, N):

    """
    See SnapPy.pyx for documentation

    >>> from snappy import Manifold
    >>> M = Manifold("4_1")
    >>> get_generalized_ptolemy_obstruction_classes(M, 2)
    [PtolemyGeneralizedObstructionClass([0, 0, 0, 0]), PtolemyGeneralizedObstructionClass([1, 0, 0, 1])]
    >>> get_generalized_ptolemy_obstruction_classes(M, 3)
    [PtolemyGeneralizedObstructionClass([0, 0, 0, 0]), PtolemyGeneralizedObstructionClass([2, 0, 0, 1])]
    >>> get_generalized_ptolemy_obstruction_classes(M, 4)
    [PtolemyGeneralizedObstructionClass([0, 0, 0, 0]), PtolemyGeneralizedObstructionClass([3, 0, 0, 1]), PtolemyGeneralizedObstructionClass([2, 0, 0, 2])]
    >>> get_generalized_ptolemy_obstruction_classes(M, 5)
    [PtolemyGeneralizedObstructionClass([0, 0, 0, 0]), PtolemyGeneralizedObstructionClass([4, 0, 0, 1])]
  
    >>> M = Manifold("m202")
    >>> len(get_generalized_ptolemy_obstruction_classes(M, 2))
    4
    >>> len(get_generalized_ptolemy_obstruction_classes(M, 3))
    5
    >>> len(get_generalized_ptolemy_obstruction_classes(M, 4))
    10
    
    >>> M = Manifold("m207")
    >>> len(get_generalized_ptolemy_obstruction_classes(M, 2))
    2
    >>> len(get_generalized_ptolemy_obstruction_classes(M, 3))
    14
    >>> len(get_generalized_ptolemy_obstruction_classes(M, 4))
    3
    """

    # Compute the obstruction classes
    H2_elements, explain_columns = get_obstruction_classes(manifold, N)

    filtered_H2_elements = []
    units = [x for x in range(N) if _gcd(x, N) == 1]

    already_seen = set()

    for H2_element in H2_elements:
        if not tuple(H2_element) in already_seen:
            filtered_H2_elements.append(H2_element)
            for u in units:
                already_seen.add(
                    tuple([(x * u) % N for x in H2_element]))

    return [PtolemyGeneralizedObstructionClass(H2_element, index = index,
                                               N = N, manifold = manifold)
            for index, H2_element
            in enumerate(filtered_H2_elements)]
    

def get_obstruction_classes(manifold, N):

    # compute boundary maps for cellular chain complex
    chain_d3, dummy_rows, dummy_columns = (
        manifold._ptolemy_equations_boundary_map_3())
    chain_d2, dummy_rows, explain_columns = (
        manifold._ptolemy_equations_boundary_map_2())

    # transpose to compute cohomology
    cochain_d2 = matrix.matrix_transpose(chain_d3)
    cochain_d1 = matrix.matrix_transpose(chain_d2)

    return (homology.homology_representatives(cochain_d2, cochain_d1, N),
            explain_columns)

class PtolemyVarietyList(list):
    def retrieve_decomposition(self, *args, **kwargs):
        return MethodMappingList(
            [ p.retrieve_decomposition(*args, **kwargs)
              for p in self ])
    def compute_decomposition(self, *args, **kwargs):
        return MethodMappingList(
            [ p.compute_decomposition(*args, **kwargs)
              for p in self ])
    def retrieve_solutions(self, *args, **kwargs):
        return MethodMappingList(
            [ p.retrieve_solutions(*args, **kwargs)
              for p in self ])
    def compute_solutions(self, *args, **kwargs):
        return MethodMappingList(
            [ p.compute_solutions(*args, **kwargs)
              for p in self ])
        

def get_ptolemy_variety(manifold, N, obstruction_class = None,
                        simplify = True, eliminate_fixed_ptolemys = False):

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

    obstruction_class --- class from Definition 1.7 of (1).
    None for trivial class or a value returned from get_ptolemy_obstruction_classes.
    Short cuts: obstruction_class = 'all' returns a list of Ptolemy varieties
    for each obstruction. For easier iteration, can set obstruction_class to 
    an integer.

    simplify --- boolean to indicate whether to simplify the equations which
    significantly reduces the number of variables.
    Simplifying means that several identified Ptolemy coordinates x = y = z = ...
    are eliminated instead of adding relations x - y = 0, y - z = 0, ...
    Defaults to True.

    eliminate_fixed_ptolemys --- boolean to indicate whether to eliminate
    the Ptolemy coordinates that are set to 1 for fixing the decoration.
    Even though this simplifies the resulting representation, setting it to
    True can cause magma to run longer when finding a Groebner basis.
    Defaults to False.

    === Examples for 4_1 ===
    
    >>> from snappy import Manifold
    >>> M = Manifold("4_1")

    Get the varieties for all obstruction classes at once (use
    help(varieties[0]) for more information):
    
    >>> varieties = get_ptolemy_variety(M, N = 2, obstruction_class = "all")

    Print the variety as an ideal (sage object) for the non-trivial class:

    >>> varieties[1].ideal    #doctest: +SKIP
    Ideal (-c_0011_0^2 + c_0011_0*c_0101_0 + c_0101_0^2, -c_0011_0^2 - c_0011_0*c_0101_0 + c_0101_0^2, c_0011_0 - 1) of Multivariate Polynomial Ring in c_0011_0, c_0101_0 over Rational Field                                                       
    (skip doctest because example only works in sage and not plain python)

    >>> for eqn in varieties[1].equations:
    ...     print("    ", eqn)
         - c_0011_0 * c_0101_0 + c_0011_0^2 + c_0101_0^2
         c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
         - 1 + c_0011_0

    Generate a magma input to compute Groebner basis for N = 3:
    
    >>> p = get_ptolemy_variety(M, N = 3)
    >>> s = p.to_magma()

    The beginning of the magma input

    >>> print(s.strip())       #doctest: +ELLIPSIS
    // Setting up the Polynomial ring and ideal
    <BLANKLINE>
    R<c_0012_0, c_0012_1, c_0102_0, c_0111_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0> := PolynomialRing(RationalField(), 8, "grevlex");
    MyIdeal := ideal<R |
              c_0012_0 * c_1101_0 + c_0102_0 * c_0111_0 - c_0102_0 * c_1011_0,
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
    >>> print(_magma_output_for_4_1__sl3)      #doctest: +ELLIPSIS
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
    Mod(-1/2*x - 3/2, x^2 + 3*x + 4)

    Example of simplified vs non-simplified variety:

    >>> simplified = get_ptolemy_variety(M, N = 4, obstruction_class = 1)
    >>> full = get_ptolemy_variety(M, N = 4, obstruction_class = 1, simplify = False)
    >>> len(simplified.variables), len(full.variables)
    (21, 63)
    >>> len(simplified.equations), len(full.equations)
    (24, 72)

    === ONLY DOCTESTS, NOT PART OF DOCUMENTATION ===

    >>> varieties = get_ptolemy_variety(M, N = 2, obstruction_class = "all", eliminate_fixed_ptolemys = True)
    
    >>> for eqn in varieties[1].equations:
    ...     print("    ", eqn)
         1 - c_0101_0 + c_0101_0^2
         - 1 + c_0101_0 - c_0101_0^2

    >>> p = get_ptolemy_variety(M, N = 3, eliminate_fixed_ptolemys = True)
    >>> s = p.to_magma()
    >>> print(s.split('ring and ideal')[1].strip())          #doctest: +ELLIPSIS
    R<c_0012_1, c_0102_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0> := PolynomialRing(RationalField(), 6, "grevlex");
    MyIdeal := ideal<R |
              c_0102_0 - c_0102_0 * c_1011_0 + c_1101_0,
        ...

    """
    
    if hasattr(manifold, 'cusp_info'):
        if False in manifold.cusp_info('is_complete'):
            raise Exception("Dehn fillings not supported by Ptolemy variety")
    
    # Typing M.ptolemy_variety(N = 3, 1) into sage makes 
    # N have type sage.rings.integer.Integer
    N = int(N)

    # If we are explicitly given an obstruction class or it is None
    # just directly call PtolemyVariety
    if ( obstruction_class is None or
         isinstance(obstruction_class, PtolemyObstructionClass) or
         isinstance(obstruction_class, PtolemyGeneralizedObstructionClass)):
        return PtolemyVariety(
            manifold, N, obstruction_class,
            simplify = simplify,
            eliminate_fixed_ptolemys = eliminate_fixed_ptolemys)

    # Will we return a list of obstruction classes versus a list
    list_obstruction_classes = False

    if obstruction_class == 'all_original':
        # Revert to old behavior of listing obstruction classes:
        # Use Z/2-obstruction class for even N
        if N % 2 == 0:
            obstruction_classes = get_ptolemy_obstruction_classes(
                manifold)
        # And no obstruction class for odd N
        else:
            obstruction_classes = [ None ]
        list_obstruction_classes = True
    elif obstruction_class == 'all_generalized':
        # Use the generalized Z/n obstruction class
        obstruction_classes = get_generalized_ptolemy_obstruction_classes(
            manifold, N)
        list_obstruction_classes = True
    else:
        # New mode: 
        if N == 2:
            # N = 2 uses Z/2-obstruction class (so that we can compute
            # complex volume)
            obstruction_classes = get_ptolemy_obstruction_classes(
                manifold)
        else:
            # N > 2 uses generalized obstruction class
            obstruction_classes = get_generalized_ptolemy_obstruction_classes(
                manifold, N)

        # List all obstruction classes
        if obstruction_class == 'all':
            list_obstruction_classes = True
            
    # Give a list of all obstruction classes
    if list_obstruction_classes:
        return PtolemyVarietyList([
            PtolemyVariety(
                manifold, N, obstruction_class,
                simplify = simplify,
                eliminate_fixed_ptolemys = eliminate_fixed_ptolemys)
            for obstruction_class in obstruction_classes])
    
    # Otherwise try to interpret obstruction_class as an index
    try:
        obstruction_class = obstruction_classes[int(obstruction_class)]
    except:
        raise Exception("Bad index for obstruction class")

    return PtolemyVariety(manifold, N, obstruction_class,
                          simplify = simplify,
                          eliminate_fixed_ptolemys = eliminate_fixed_ptolemys)

